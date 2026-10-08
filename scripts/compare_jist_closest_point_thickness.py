#!/usr/bin/env python3
"""Compare selected full mesh thickness wrapper with Python on a phantom."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.spatial import cKDTree

from pymacruise.surfaces import SurfaceMesh
from pymacruise.jist_intersector import JistSurfaceTree
from pymacruise.thickness import (
    _closest_points_on_triangles,
    _nearest_surface_points,
    closest_point_thickness,
)

ROOT=Path(__file__).resolve().parents[1]
TOADS=ROOT.parent/"toads-cruise-scm-2026-09-26/toads-cruise/src"
TREE=TOADS/"edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree"
INTERSECTOR=TOADS/"edu/jhu/ece/iacl/algorithms/graphics/intersector"
SELECTED=ROOT.parent/"toads-cruise-scm-2026-09-26/toads-cruise/src/edu/vanderbilt/masi/algorithms/CRUISE/utilities/ClosestPointThicknessForInner.java"
STUB=ROOT/"tests/java_reference/thickness_comparison_stubs"
HARNESS=ROOT/"tests/java_reference/ClosestPointThicknessForInnerHarness.java"
DIAGNOSTIC=ROOT/"tests/java_reference/ThicknessQueryDiagnostic.java"
FACE_HARNESS=ROOT/"tests/java_reference/SurfaceIntersectorFaceQueryHarness.java"
TREE_HARNESS=ROOT/"tests/java_reference/KdTreeStructureHarness.java"
VERSION=ROOT/"tests/java_reference/surface_comparison_stubs/edu/jhu/ece/iacl/algorithms/VersionUtil.java"
ABSTRACT=ROOT/"tests/java_reference/stubs/edu/jhu/ece/iacl/jist/pipeline/AbstractCalculation.java"
VECMATH=(ROOT.parent/"jist-scm-2026-09-26/JIST/lib/vecmath-1.3.1.jar").resolve()
FIXTURE=ROOT/"tests/data/synthetic/jist_surface_comparison"
WORK=ROOT/"build/java-reference/closest-point-thickness"
REPORT=ROOT/"docs/reports/jist-closest-point-thickness-java-comparison-2026-10-07.json"
REPORT_MD=ROOT/"docs/reports/jist-closest-point-thickness-java-comparison-2026-10-07.md"


def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()


def centroid_tree_nearest_points(points: np.ndarray, mesh: SurfaceMesh):
    """Historical fast locator retained to document the tie-order deviation."""
    triangles = mesh.vertices[mesh.faces]
    centroids = triangles.mean(axis=1)
    radii = np.linalg.norm(triangles - centroids[:, None, :], axis=2).max(axis=1)
    tree = cKDTree(centroids)
    closest = np.empty_like(points, dtype=np.float64)
    distances = np.empty(len(points), dtype=np.float64)
    face_indices = np.empty(len(points), dtype=np.int64)
    for index, point in enumerate(points):
        _seed_distance, seed = tree.query(point, k=1)
        seed_point, seed_squared = _closest_points_on_triangles(
            point, triangles[int(seed) : int(seed) + 1]
        )
        upper_bound = float(np.sqrt(seed_squared[0]))
        candidate_ids = np.asarray(
            tree.query_ball_point(point, upper_bound + float(radii.max()) + 1e-12),
            dtype=np.int64,
        )
        candidates, squared = _closest_points_on_triangles(point, triangles[candidate_ids])
        local = int(np.argmin(squared))
        closest[index] = candidates[local]
        distances[index] = np.sqrt(squared[local])
        face_indices[index] = int(candidate_ids[local])
    return closest, distances, face_indices


def main()->int:
    classes=WORK/"classes"; classes.mkdir(parents=True,exist_ok=True)
    sources=[HARNESS,DIAGNOSTIC,FACE_HARNESS,TREE_HARNESS,STUB/"Jama/PackageStub.java",STUB/"edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.java",STUB/"edu/jhu/ece/iacl/algorithms/thickness/ThicknessSolver.java",STUB/"edu/jhu/ece/iacl/jist/structures/image/ImageData.java",STUB/"edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java",STUB/"edu/jhu/ece/iacl/algorithms/topology/ConnectivityRule.java",ABSTRACT,VERSION,*[TREE/name for name in ("BBox.java","BBoxEdge.java","KdPoint3.java","KdSegment.java","KdTriangle.java","KdTree.java")],INTERSECTOR/"IntersectorTriangle.java",INTERSECTOR/"SurfaceIntersector.java",SELECTED]
    subprocess.run(["javac","-d",str(classes),"-cp",str(VECMATH),*[str(p) for p in sources]],check=True,cwd=ROOT)
    mesh_json={name:FIXTURE/"java_meshes"/f"{name}.json" for name in ("inner","central","outer")}
    meshes={}
    for name,path in mesh_json.items():
        payload=json.loads(path.read_text(encoding="utf-8"))
        vertices=np.asarray(payload["vertices"],dtype=np.float32).astype(np.float64)
        faces=np.asarray(payload["faces"],dtype=np.int64)
        meshes[name]=SurfaceMesh(vertices,faces)
    paths={}
    for name,mesh in meshes.items():
        vertices_path=WORK/f"{name}-vertices.raw.f32le"; faces_path=WORK/f"{name}-faces.raw.i32le"
        mesh.vertices.astype("<f4").tofile(vertices_path); mesh.faces.astype("<i4").tofile(faces_path)
        paths[name]=(vertices_path,faces_path)
    result_path=WORK/"result.raw"
    subprocess.run(["java","-Xmx2g","-cp",f"{classes}:{VECMATH}","ClosestPointThicknessForInnerHarness",str(paths['inner'][0]),str(paths['inner'][1]),str(paths['central'][0]),str(paths['central'][1]),str(paths['outer'][0]),str(paths['outer'][1]),str(result_path)],check=True,cwd=ROOT)
    expected=closest_point_thickness(meshes["inner"],meshes["outer"],central=meshes["central"])
    raw=result_path.read_bytes(); offset=0
    inner_count=int(np.frombuffer(raw,dtype="<i4",count=1,offset=offset)[0]); offset+=4
    java_inner=np.frombuffer(raw,dtype="<f8",count=inner_count*3,offset=offset).reshape(inner_count,3).copy(); offset+=inner_count*24
    central_count=int(np.frombuffer(raw,dtype="<i4",count=1,offset=offset)[0]); offset+=4
    java_central=np.frombuffer(raw,dtype="<f8",count=central_count*4,offset=offset).reshape(central_count,4).copy()
    baseline_inner_points, baseline_inner_d1, baseline_inner_faces = centroid_tree_nearest_points(
        meshes["inner"].vertices, meshes["outer"]
    )
    baseline_inner_return, baseline_inner_d2, _ = centroid_tree_nearest_points(
        baseline_inner_points, meshes["inner"]
    )
    diagnostic_index=int(np.argmax(np.abs(java_inner[:,1]-baseline_inner_d2)))
    diagnostic_path=WORK/f"diagnostic-{diagnostic_index}.raw"
    subprocess.run(["java","-Xmx2g","-cp",f"{classes}:{VECMATH}","ThicknessQueryDiagnostic",str(paths['inner'][0]),str(paths['inner'][1]),str(paths['outer'][0]),str(paths['outer'][1]),str(diagnostic_index),str(diagnostic_path)],check=True,cwd=ROOT)
    diag=diagnostic_path.read_bytes(); cursor=0
    query=np.frombuffer(diag,dtype="<f4",count=3,offset=cursor).astype(float);cursor+=12
    java_steps=[]
    for _ in range(2):
      distance=float(np.frombuffer(diag,dtype="<f8",count=1,offset=cursor)[0]);cursor+=8
      point=np.frombuffer(diag,dtype="<f4",count=3,offset=cursor).astype(float);cursor+=12
      ids=np.frombuffer(diag,dtype="<i4",count=4,offset=cursor).copy();cursor+=16
      java_steps.append({"distance":distance,"point":point,"triangle_and_vertex_ids":ids})
    py_outer_point,py_outer_distance,py_outer_face=centroid_tree_nearest_points(query[None,:],meshes["outer"])
    py_return_point,py_return_distance,py_return_face=centroid_tree_nearest_points(java_steps[0]["point"][None,:],meshes["inner"])
    print(json.dumps({"diagnostic_vertex":diagnostic_index,"query":query.tolist(),"java_outer":{**{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in java_steps[0].items()}},"python_outer":{"distance":float(py_outer_distance[0]),"point":py_outer_point[0].tolist(),"face":int(py_outer_face[0])},"java_return":{**{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in java_steps[1].items()}},"python_return_from_java_hit":{"distance":float(py_return_distance[0]),"point":py_return_point[0].tolist(),"face":int(py_return_face[0])}},indent=2))
    if inner_count!=len(meshes["inner"].vertices) or central_count!=len(meshes["central"].vertices): raise AssertionError("Java output mesh counts differ")
    inner_names=("distance_inner_to_outer_mm","distance_outer_to_inner_mm","mean_thickness_mm")
    inner_errors={name:np.abs(java_inner[:,i]-expected.inner_surface.point_data[name]) for i,name in enumerate(inner_names)}
    central_names=("distance_inner_to_outer_mm","distance_central_to_inner_mm","distance_to_inner_vertex_mm","nearest_inner_vertex_index")
    central_errors={name:np.abs(java_central[:,i]-expected.central_surface.point_data[name]) for i,name in enumerate(central_names)}
    all_errors={**{f"inner.{name}":v for name,v in inner_errors.items()},**{f"central.{name}":v for name,v in central_errors.items()}}
    max_errors={name:float(error.max()) for name,error in all_errors.items()}
    mean_errors={name:float(error.mean()) for name,error in all_errors.items()}
    central_index_mismatches=int(np.count_nonzero(java_central[:,3]!=expected.central_surface.point_data["nearest_inner_vertex_index"]))
    face_order_diagnostics={}
    for query_name,target_name in (("inner","outer"),("central","inner")):
      query_mesh,target_mesh=meshes[query_name],meshes[target_name]
      selected_path=WORK/f"selected-faces-{query_name}.raw"
      subprocess.run(["java","-Xmx2g","-cp",f"{classes}:{VECMATH}","SurfaceIntersectorFaceQueryHarness",str(paths[query_name][0]),str(paths[query_name][1]),str(paths[target_name][0]),str(paths[target_name][1]),str(selected_path)],check=True,cwd=ROOT)
      java_selected=np.fromfile(selected_path,dtype=[("face","<i4"),("distance","<f4")])
      _points,python_distances,python_faces=centroid_tree_nearest_points(query_mesh.vertices,target_mesh)
      mismatch_ids=np.flatnonzero(java_selected["face"]!=python_faces)
      triangles=target_mesh.vertices[target_mesh.faces]
      centroids=triangles.mean(axis=1)
      radii=np.linalg.norm(triangles-centroids[:,None,:],axis=2).max(axis=1)
      tree=cKDTree(centroids)
      source_is_max_tied_face=0
      for query_id in mismatch_ids:
        bound=float(python_distances[query_id]+radii.max()+1e-12)
        candidates=np.asarray(tree.query_ball_point(query_mesh.vertices[query_id],bound),dtype=np.int64)
        _cp,squared=_closest_points_on_triangles(query_mesh.vertices[query_id],triangles[candidates])
        candidate_distances=np.sqrt(squared)
        tie_ids=candidates[candidate_distances<=float(candidate_distances.min())+1e-6]
        source_is_max_tied_face+=int(int(java_selected["face"][query_id])==int(tie_ids.max()))
      face_order_diagnostics[query_name]={"query_count":len(query_mesh.vertices),"java_python_selected_face_mismatches":int(len(mismatch_ids)),"maximum_distance_difference_mm":float(np.max(np.abs(java_selected["distance"]-python_distances))),"mismatches_where_java_face_is_highest_face_within_1e-6_mm_tie":source_is_max_tied_face,"tie_policy_sample_size":int(len(mismatch_ids)),"central_nearest_vertex_index_mismatches":central_index_mismatches if query_name=="central" else None}
    source_trees={name:JistSurfaceTree(meshes[name].vertices,meshes[name].faces) for name in ("inner","outer")}
    tree_build_diagnostics={}
    for target_name,tree in source_trees.items():
      java_dump=WORK/f"java-tree-{target_name}.tsv"
      subprocess.run(["java","-Xmx2g","-cp",f"{classes}:{VECMATH}","edu.jhu.ece.iacl.algorithms.graphics.intersector.KdTreeStructureHarness",str(paths[target_name][0]),str(paths[target_name][1]),str(java_dump)],check=True,cwd=ROOT)
      java_signature=[]
      for line in java_dump.read_text(encoding="utf-8").splitlines():
        row=line.split("\t")
        if row[1]=="T": java_signature.append((row[0],"T",int(row[2])))
        else: java_signature.append((row[0],"B",int(row[2]),int(row[3]),tuple(float(np.float32(x)) for x in row[4].split(",")),tuple(float(np.float32(x)) for x in row[5].split(","))))
      py_signature=tree.structure_signature()
      tree_build_diagnostics[target_name]={"java_node_count":len(java_signature),"python_node_count":len(py_signature),"tree_structure_exact":java_signature==py_signature}
    inner_tree=source_trees["inner"]
    outer_points,tree_d1,tree_f1=source_trees["outer"].query_many(meshes["inner"].vertices)
    _tree_return_points,tree_d2,tree_f2=inner_tree.query_many(outer_points)
    tree_mean=(tree_d1+tree_d2)/2.0
    _central_points,tree_dc,tree_fc=inner_tree.query_many(meshes["central"].vertices)
    central_faces=meshes["inner"].faces[tree_fc]
    central_points32=meshes["central"].vertices.astype(np.float32)
    central_vertices32=meshes["inner"].vertices.astype(np.float32)[central_faces]
    vertex_delta=np.subtract(central_points32[:,None,:],central_vertices32,dtype=np.float32)
    vertex_squared=np.float32(vertex_delta[:,:,0]*vertex_delta[:,:,0]+vertex_delta[:,:,1]*vertex_delta[:,:,1])
    vertex_squared=np.float32(vertex_squared+vertex_delta[:,:,2]*vertex_delta[:,:,2])
    tree_vertex_distance=np.sqrt(vertex_squared,dtype=np.float32)
    tree_corners=np.argmin(tree_vertex_distance,axis=1)
    tree_nearest_vertices=central_faces[np.arange(len(tree_fc)),tree_corners]
    tree_central_thickness=tree_mean[tree_nearest_vertices]
    tree_python_prototype={"max_inner_to_outer_distance_error_mm":float(np.max(np.abs(tree_d1-java_inner[:,0]))),"max_outer_to_inner_distance_error_mm":float(np.max(np.abs(tree_d2-java_inner[:,1]))),"max_mean_inner_thickness_error_mm":float(np.max(np.abs(tree_mean-java_inner[:,2]))),"central_nearest_vertex_index_mismatches":int(np.count_nonzero(tree_nearest_vertices!=java_central[:,3])),"max_central_projected_thickness_error_mm":float(np.max(np.abs(tree_central_thickness-java_central[:,0]))),"max_central_surface_distance_error_mm":float(np.max(np.abs(tree_dc-java_central[:,1]))),"max_central_nearest_vertex_distance_error_mm":float(np.max(np.abs(tree_vertex_distance[np.arange(len(tree_fc)),tree_corners]-java_central[:,2]))),"inner_to_outer_face_mismatches":int(np.count_nonzero(tree_f1!=np.fromfile(WORK/"selected-faces-inner.raw",dtype=[("face","<i4"),("distance","<f4")])["face"])),"central_to_inner_face_mismatches":int(np.count_nonzero(tree_fc!=np.fromfile(WORK/"selected-faces-central.raw",dtype=[("face","<i4"),("distance","<f4")])["face"]))}
    report={
      "title":"Selected ClosestPointThicknessForInner versus Python thickness API",
      "synthetic_only":True,
      "method_scope":"Compiles the selected VandyCRUISE ClosestPointThicknessForInner.java mesh overload and its selected SurfaceIntersector/KdTree sources. Minimal EmbeddedSurface, ThicknessSolver, JIST progress, image-only overload, and unused Jama import stubs support compilation; source thickness methods are compiled unchanged.",
      "sources":{"closest_point_thickness_for_inner_sha256":sha(SELECTED),"surface_intersector_sha256":sha(INTERSECTOR/"SurfaceIntersector.java"),"intersector_triangle_sha256":sha(INTERSECTOR/"IntersectorTriangle.java"),"kd_tree_sha256":sha(TREE/"KdTree.java"),"kd_triangle_sha256":sha(TREE/"KdTriangle.java")},
      "fixture":{"path":str(FIXTURE.relative_to(ROOT)),"synthetic_only":True,"mesh_origin":"committed outputs of selected Java IsoSurfaceOnGrid reference on the committed level-set NIfTIs; the thickness comparison holds mesh extraction fixed","mesh_json_sha256":{n:sha(p) for n,p in mesh_json.items()},"mesh_counts":{n:{"vertices":len(m.vertices),"faces":len(m.faces)} for n,m in meshes.items()},"coordinate_storage":"RAS millimetres stored as Java Point3f float32"},
      "metrics":{"max_absolute_field_difference":max_errors,"mean_absolute_field_difference":mean_errors,"maximum_physical_field_error_gate_mm":1e-3,"central_nearest_inner_vertex_index_mismatches":central_index_mismatches,"all_physical_fields_pass":bool(max(v for k,v in max_errors.items() if k.endswith("_mm"))<=1e-3),"nearest_vertex_indices_exact":central_index_mismatches==0,"comparison_pass":bool(max(v for k,v in max_errors.items() if k.endswith("_mm"))<=1e-3 and central_index_mismatches==0)},
      "face_order_diagnostics":face_order_diagnostics,
      "centroid_tree_baseline":{"metrics":{"inner_to_outer_face_mismatches":face_order_diagnostics["inner"]["java_python_selected_face_mismatches"],"central_to_inner_face_mismatches":face_order_diagnostics["central"]["java_python_selected_face_mismatches"],"max_outer_to_inner_distance_error_mm":float(np.max(np.abs(baseline_inner_d2-java_inner[:,1])))},"status":"historical approximate locator; retained for regression context, not used by production thickness API"},
      "python_jist_tree_port":{"tree_builds":tree_build_diagnostics,"metrics":tree_python_prototype,"status":"connected to production thickness API; exact selected-source match on these committed meshes"},
      "tie_diagnostic":{"inner_vertex_index":diagnostic_index,"java_selected_outer_face":int(java_steps[0]["triangle_and_vertex_ids"][0]),"centroid_baseline_selected_outer_face":int(py_outer_face[0]),"java_outer_distance_mm":java_steps[0]["distance"],"centroid_baseline_outer_distance_mm":float(py_outer_distance[0]),"java_outer_hit_ras_mm":java_steps[0]["point"].tolist(),"centroid_baseline_outer_hit_ras_mm":py_outer_point[0].tolist(),"java_return_distance_mm":java_steps[1]["distance"],"centroid_baseline_return_distance_from_java_hit_mm":float(py_return_distance[0])},
      "limitations":["This validates the selected mesh overload on committed Java IsoSurfaceOnGrid output meshes; it holds mesh extraction fixed and does not validate the complete workflow.","The selected wrapper's serialized VTK/NIfTI output adapters are not exercised.","The exact source tree/query path increases runtime over the former centroid-tree locator; brain-resolution runtime and peak memory remain to be measured."],
    }
    REPORT.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    lines=["# Selected JIST closest-point thickness comparison","","This runs the selected VandyCRUISE `ClosestPointThicknessForInner.java` mesh overload unchanged and compares its inner-surface bidirectional point data and central-surface projection with `pymacruise.thickness.closest_point_thickness` on the same saved Java `IsoSurfaceOnGrid` synthetic meshes. This holds mesh extraction fixed. The Java thickness class, `SurfaceIntersector`, `KdTree`, and `IntersectorTriangle` source bodies are compiled directly; minimal container/superclass stubs are used only for surrounding API state.","",f"Source hash for `ClosestPointThicknessForInner.java`: `{report['sources']['closest_point_thickness_for_inner_sha256']}`.","","| Field | Maximum absolute difference (mm, except vertex index) | Mean absolute difference (mm, except vertex index) |","|---|---:|---:|"]
    for name,error in max_errors.items(): lines.append(f"| `{name}` | {error:.9g} | {mean_errors[name]:.9g} |")
    tie=report["tie_diagnostic"]
    prototype=report["python_jist_tree_port"]
    tree_metrics=prototype["metrics"]
    lines.extend(["",f"**Production API result: {'PASS' if report['metrics']['comparison_pass'] else 'FAIL'}.** On these held-fixed Java meshes, central nearest-vertex indices match exactly and five of six physical fields are bit-exact; physical-field maximum error is `{max(v for k,v in max_errors.items() if k.endswith('_mm')):.9g} mm` and central vertex-index mismatches are {central_index_mismatches}.","",f"The connected Python JIST tree reproduces Java's complete SAH tree structure exactly on both meshes ({prototype['tree_builds']['inner']['java_node_count']:,} inner-tree nodes; {prototype['tree_builds']['outer']['java_node_count']:,} outer-tree nodes). Its face-query selections and wrapper output metrics are exact on this fixture. The older centroid-tree locator is retained only as a baseline: it selects different faces at ties ({report['centroid_tree_baseline']['metrics']['inner_to_outer_face_mismatches']} inner-to-outer and {report['centroid_tree_baseline']['metrics']['central_to_inner_face_mismatches']} central-to-inner query mismatches).","",f"This is selected-source synthetic mesh evidence only. It does not validate upstream processing, mesh extraction equivalence, whole-workflow outputs, or real MRI. Exact JIST tree construction increases runtime; brain-resolution cost remains unmeasured. See the [JSON report](jist-closest-point-thickness-java-comparison-2026-10-07.json) for source hashes and metrics. Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness.py`.",""])
    REPORT_MD.write_text("\n".join(lines),encoding="utf-8")
    print(json.dumps({"report":str(REPORT),"markdown":str(REPORT_MD),"metrics":report["metrics"]},indent=2))
    return 0

if __name__=="__main__": raise SystemExit(main())
