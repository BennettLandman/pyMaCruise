# GIFTI surfaces and surface information

PyMaCRUISE reads and writes GIFTI triangular surfaces (`.surf.gii`, or generic
`.gii`) through NiBabel, which is already a runtime dependency. Geometry and
named per-vertex arrays can live in one file. Dense measurements can also be
exported separately as `.shape.gii`, `.func.gii`, or `.label.gii` files.

## Commands

Extract directly from a level set:

```sh
pymacruise surface-extract --level-set inner.nii.gz --output inner.surf.gii
```

Convert an existing BASE or PLUS workflow output, retaining its vertex data:

```sh
pymacruise surface-convert --input inner_surface_thickness.vtk \
  --output inner.surf.gii --data-output thickness.shape.gii \
  --provenance conversion.json
pymacruise surface-info --input inner.surf.gii --output surface-info.json
```

`surface-info` prints JSON with vertex and triangle counts, RAS-mm bounds,
surface metadata, reference-space code, vertex field names/shapes/intents,
field metadata, and label names. For an older untagged VTK input, specify
`--input-frame RAS` or `--input-frame LPS` on conversion or inspection.

Attach separate measurements before conversion:

```sh
pymacruise surface-convert --input cortex.surf.gii \
  --point-data thickness.shape.gii --point-data regions.label.gii \
  --output cortex-with-data.surf.gii
```

Each array must have one row per vertex **in the same vertex order**. Matching
counts alone cannot establish correspondence. Duplicate names are rejected.
Sparse `NODE_INDEX` data and files with multiple geometries are unsupported.
Unnamed arrays receive `array_<index>` names in original file order.

`thickness-surface` accepts either VTK or GIFTI on each input/output:

```sh
pymacruise thickness-surface --inner inner.surf.gii --outer outer.surf.gii \
  --inner-output inner-thickness.surf.gii
```

`base-surface-thickness --surface-format gii` writes its four surface products
as `.gii`; the default remains `vtk`. Complete restartable BASE/PLUS runs
retain their existing VTK output contracts. Use `surface-convert` to export
their surfaces afterward.

## Python API

```python
from pymacruise import (
    read_surface, write_surface, load_gifti_point_data, write_gifti_point_data,
)

mesh = read_surface("inner_surface_thickness.vtk")
write_surface(mesh, "inner.surf.gii")
write_gifti_point_data(mesh, "thickness.shape.gii", names=["mean_thickness_mm"])
mesh = load_gifti_point_data(read_surface("inner.surf.gii"), "regions.label.gii")
```

`read_gifti` and `write_gifti` are also public. `SurfaceMesh.point_data` holds
finite scalar or multicomponent vertex arrays. Its `metadata`,
`point_data_metadata`, `point_data_intents`, `label_table`, and
`gifti_dataspace` fields retain GIFTI surface/measurement information. File and
pointset metadata are merged on import (pointset entries take precedence).
Measurement file metadata fills entries absent from the geometry's metadata.
GIFTI label arrays retain their integer intent and label names/RGBA table.
New scalar measurements default to `NIFTI_INTENT_SHAPE`; multicomponent arrays
default to `NIFTI_INTENT_NONE` unless an intent is supplied.

## Coordinates, precision, and validation scope

All meshes use RAS+ millimetres. A GIFTI pointset's coordinate transform is
applied once on import; exports store an identity transform in the resulting
reference space. Reflections reverse triangle winding. Nonidentity transforms
with an unknown target space are rejected. An unspecified space remains
unknown, without asserting scanner or template registration. Per-vertex
arrays stay attached by index; spatial vector/tensor components are not
reoriented by these geometry transforms.

Exports use standard float32 coordinates/measurements and int32 triangles or
labels with compressed base64 encoding. Import also supports ASCII and
uncompressed base64 through NiBabel. Float32 export can round binary64 geometry
and measurements, potentially affecting near-tied closest-point queries. Keep
VTK originals for comparisons requiring binary64 geometry. VTK conversion
retains numeric vertex arrays but its existing writer does not carry GIFTI
metadata, intents, or label tables.

Tests cover geometry and measurement round trips, label tables, encoding
variants, transforms and reference spaces, malformed data, conversion,
inspection, extraction, and surface thickness. These are engineering and I/O
checks; they do not establish whole-workflow Java equivalence or clinical
accuracy.

Format API reference: [NiBabel GIFTI documentation](https://nipy.org/nibabel/reference/nibabel.gifti.html).
