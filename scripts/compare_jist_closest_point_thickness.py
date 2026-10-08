#!/usr/bin/env python3
"""Compare selected full mesh thickness wrapper with Python on a phantom."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface
from pymacruise.thickness import closest_point_thickness

ROOT=Path(__file__).resolve().parents[1]
TOADS=ROOT.parent/"toads-cruise-scm-2026-09-26/toads-cruise/src"
TREE=TOADS/"edu/jhu/ece/iacl/algorithms/graphics/locator/kdtree"
INTERSECTOR=TOADS/"edu/jhu/ece/iacl/algorithms/graphics/intersector"
SELECTED=ROOT.parent/"toads-cruise-scm-2026-09-26/toads-cruise/src/edu/vanderbilt/masi/algorithms/CRUISE/utilities/ClosestPointThicknessForInner.java"
STUB=ROOT/"tests/java_reference/thickness_comparison_stubs"
HARNESS=ROOT/"tests/java_reference/ClosestPointThicknessForInnerHarness.java"
DIAGNOSTIC=ROOT/"tests/java_reference/ThicknessQueryDiagnostic.java"
VERSION=ROOT/"tests/java_reference/surface_comparison_stubs/edu/jhu/ece/iacl/algorithms/VersionUtil.java"
ABSTRACT=ROOT/"tests/java_reference/stubs/edu/jhu/ece/iacl/jist/pipeline/AbstractCalculation.java"
VECMATH=(ROOT.parent/"jist-scm-2026-09-26/JIST/lib/vecmath-1.3.1.jar").resolve()
FIXTURE=ROOT/"tests/data/synthetic/jist_surface_comparison"
WORK=ROOT/"build/java-reference/closest-point-thickness"
REPORT=ROOT/"docs/reports/jist-closest-point-thickness-java-comparison-2026-10-07.json"
REPORT_MD=ROOT/"docs/reports/jist-closest-point-thickness-java-comparison-2026-10-07.md"


def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()


def main()->int:
    classes=WORK/"classes"; classes.mkdir(parents=True,exist_ok=True)
    sources=[HARNESS,DIAGNOSTIC,STUB/"Jama/PackageStub.java",STUB/"edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.java",STUB/"edu/jhu/ece/iacl/algorithms/thickness/ThicknessSolver.java",STUB/"edu/jhu/ece/iacl/jist/structures/image/ImageData.java",STUB/"edu/jhu/ece/iacl/algorithms/graphics/isosurf/IsoSurfaceOnGrid.java",STUB/"edu/jhu/ece/iacl/algorithms/topology/ConnectivityRule.java",ABSTRACT,VERSION,*[TREE/name for name in ("BBox.java","BBoxEdge.java","KdPoint3.java","KdSegment.java","KdTriangle.java","KdTree.java")],INTERSECTOR/"IntersectorTriangle.java",INTERSECTOR/"SurfaceIntersector.java",SELECTED]
    subprocess.run(["javac","-d",str(classes),"-cp",str(VECMATH),*[str(p) for p in sources]],check=True,cwd=ROOT)
    images={name:nib.load(str(FIXTURE/f"{name}.nii.gz")) for name in ("inner","central","outer")}
    meshes={name:extract_levelset_surface(image) for name,image in images.items()}
    meshes={name:SurfaceMesh(mesh.vertices.astype(np.float32).astype(np.float64),mesh.faces) for name,mesh in meshes.items()}
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
    diagnostic_path=WORK/"diagnostic-1124.raw"
    subprocess.run(["java","-Xmx2g","-cp",f"{classes}:{VECMATH}","ThicknessQueryDiagnostic",str(paths['inner'][0]),str(paths['inner'][1]),str(paths['outer'][0]),str(paths['outer'][1]),"1124",str(diagnostic_path)],check=True,cwd=ROOT)
    diag=diagnostic_path.read_bytes(); cursor=0
    query=np.frombuffer(diag,dtype="<f4",count=3,offset=cursor).astype(float);cursor+=12
    java_steps=[]
    for _ in range(2):
      distance=float(np.frombuffer(diag,dtype="<f8",count=1,offset=cursor)[0]);cursor+=8
      point=np.frombuffer(diag,dtype="<f4",count=3,offset=cursor).astype(float);cursor+=12
      ids=np.frombuffer(diag,dtype="<i4",count=4,offset=cursor).copy();cursor+=16
      java_steps.append({"distance":distance,"point":point,"triangle_and_vertex_ids":ids})
    py_outer_point,py_outer_distance,py_outer_face=__import__("pymacruise.thickness",fromlist=["_nearest_surface_points"])._nearest_surface_points(query[None,:],meshes["outer"])
    py_return_point,py_return_distance,py_return_face=__import__("pymacruise.thickness",fromlist=["_nearest_surface_points"])._nearest_surface_points(java_steps[0]["point"][None,:],meshes["inner"])
    print(json.dumps({"diagnostic_vertex":1124,"query":query.tolist(),"java_outer":{**{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in java_steps[0].items()}},"python_outer":{"distance":float(py_outer_distance[0]),"point":py_outer_point[0].tolist(),"face":int(py_outer_face[0])},"java_return":{**{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in java_steps[1].items()}},"python_return_from_java_hit":{"distance":float(py_return_distance[0]),"point":py_return_point[0].tolist(),"face":int(py_return_face[0])}},indent=2))
    if inner_count!=len(meshes["inner"].vertices) or central_count!=len(meshes["central"].vertices): raise AssertionError("Java output mesh counts differ")
    inner_names=("distance_inner_to_outer_mm","distance_outer_to_inner_mm","mean_thickness_mm")
    inner_errors={name:np.abs(java_inner[:,i]-expected.inner_surface.point_data[name]) for i,name in enumerate(inner_names)}
    central_names=("distance_inner_to_outer_mm","distance_central_to_inner_mm","distance_to_inner_vertex_mm","nearest_inner_vertex_index")
    central_errors={name:np.abs(java_central[:,i]-expected.central_surface.point_data[name]) for i,name in enumerate(central_names)}
    all_errors={**{f"inner.{name}":v for name,v in inner_errors.items()},**{f"central.{name}":v for name,v in central_errors.items()}}
    max_errors={name:float(error.max()) for name,error in all_errors.items()}
    mean_errors={name:float(error.mean()) for name,error in all_errors.items()}
    central_index_mismatches=int(np.count_nonzero(java_central[:,3]!=expected.central_surface.point_data["nearest_inner_vertex_index"]))
    report={
      "title":"Selected ClosestPointThicknessForInner versus Python thickness API",
      "synthetic_only":True,
      "method_scope":"Compiles the selected VandyCRUISE ClosestPointThicknessForInner.java mesh overload and its selected SurfaceIntersector/KdTree sources. Minimal EmbeddedSurface, ThicknessSolver, JIST progress, image-only overload, and unused Jama import stubs support compilation; source thickness methods are compiled unchanged.",
      "sources":{"closest_point_thickness_for_inner_sha256":sha(SELECTED),"surface_intersector_sha256":sha(INTERSECTOR/"SurfaceIntersector.java"),"intersector_triangle_sha256":sha(INTERSECTOR/"IntersectorTriangle.java"),"kd_tree_sha256":sha(TREE/"KdTree.java"),"kd_triangle_sha256":sha(TREE/"KdTriangle.java")},
      "fixture":{"path":str(FIXTURE.relative_to(ROOT)),"synthetic_only":True,"mesh_counts":{n:{"vertices":len(m.vertices),"faces":len(m.faces)} for n,m in meshes.items()},"coordinate_storage":"RAS millimetres rounded to Java Point3f float32"},
      "metrics":{"max_absolute_field_difference":max_errors,"mean_absolute_field_difference":mean_errors,"maximum_physical_field_error_gate_mm":1e-3,"central_nearest_inner_vertex_index_mismatches":central_index_mismatches,"all_physical_fields_pass":bool(max(v for k,v in max_errors.items() if k.endswith("_mm"))<=1e-3),"nearest_vertex_indices_exact":central_index_mismatches==0,"comparison_pass":bool(max(v for k,v in max_errors.items() if k.endswith("_mm"))<=1e-3 and central_index_mismatches==0)},
      "tie_diagnostic":{"inner_vertex_index":1124,"java_selected_outer_face":int(java_steps[0]["triangle_and_vertex_ids"][0]),"python_selected_outer_face":int(py_outer_face[0]),"java_outer_distance_mm":java_steps[0]["distance"],"python_outer_distance_mm":float(py_outer_distance[0]),"java_outer_hit_ras_mm":java_steps[0]["point"].tolist(),"python_outer_hit_ras_mm":py_outer_point[0].tolist(),"java_return_distance_mm":java_steps[1]["distance"],"python_return_distance_from_java_hit_mm":float(py_return_distance[0])},
      "limitations":["This validates the selected mesh overload on Lewiner-generated input meshes, not JIST-generated meshes; the latter are a separately measured substitution.","The selected wrapper's serialized VTK/NIfTI output adapters and full pipeline are not exercised.","The overall comparison fails: Java and Python choose different equidistant outer faces at an exact folded-surface tie, producing a 0.0457 mm return-distance difference. The central projection also has 21 nearest-triangle/vertex selection differences and up to 0.307 mm projected-thickness difference. Java uses float32 triangle arithmetic plus source KD traversal ordering; Python uses float64 geometry and centroid-tree ordering. This report records the unresolved source-ordering deviation; it is not hidden by the looser surface-position tolerances."],
    }
    REPORT.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    lines=["# Selected JIST closest-point thickness comparison","","This runs the selected VandyCRUISE `ClosestPointThicknessForInner.java` mesh overload unchanged and compares its inner-surface bidirectional point data and central-surface projection with `pymacruise.thickness.closest_point_thickness` on the same float32-stored synthetic meshes. The Java class, `SurfaceIntersector`, `KdTree`, and `IntersectorTriangle` source bodies are compiled directly; minimal container/superclass stubs are used only for surrounding API state.","",f"Source hash for `ClosestPointThicknessForInner.java`: `{report['sources']['closest_point_thickness_for_inner_sha256']}`.","","| Field | Maximum absolute difference (mm, except vertex index) | Mean absolute difference (mm, except vertex index) |","|---|---:|---:|"]
    for name,error in max_errors.items(): lines.append(f"| `{name}` | {error:.9g} | {mean_errors[name]:.9g} |")
    tie=report["tie_diagnostic"]
    lines.extend(["",f"**Result: FAIL.** The wrapper exceeds the focused `0.001 mm` maximum-error gate: the return-distance error peaks at `{max_errors['inner.distance_outer_to_inner_mm']:.6f} mm` and central projected-thickness error at `{max_errors['central.distance_inner_to_outer_mm']:.6f} mm`; {central_index_mismatches} central nearest-inner-vertex indices differ.","",f"At inner vertex {tie['inner_vertex_index']}, Java selects outer face {tie['java_selected_outer_face']} and Python face {tie['python_selected_outer_face']}. Both are equidistant (`{tie['java_outer_distance_mm']:.9f} mm` in Java; `{tie['python_outer_distance_mm']:.9f} mm` in Python), but project to opposite sides of the folded outer mesh. Returning from Java's hit gives Python `{tie['python_return_distance_from_java_hit_mm']:.9f} mm`, matching Java `{tie['java_return_distance_mm']:.9f} mm`; Python's own selected hit instead leads to the reported return-distance discrepancy. This isolates the source of the largest inner error to tie ordering, not triangle-distance arithmetic.","","The current production operator remains a Python geometric implementation and does not reproduce Java's float32 triangle arithmetic and KD traversal ordering at exact ties. The report does not establish workflow equivalence or real-MRI behavior. It uses Lewiner meshes to isolate thickness behavior; surface topology remains a separate P37 deviation. The compact output with source hashes and mesh counts is in the [JSON report](jist-closest-point-thickness-java-comparison-2026-10-07.json). Reproduce with `./.venv/bin/python scripts/compare_jist_closest_point_thickness.py`.",""])
    REPORT_MD.write_text("\n".join(lines),encoding="utf-8")
    print(json.dumps({"report":str(REPORT),"markdown":str(REPORT_MD),"metrics":report["metrics"]},indent=2))
    return 0

if __name__=="__main__": raise SystemExit(main())
