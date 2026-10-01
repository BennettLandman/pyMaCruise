# TOADS startup relation state

This fixture checks the selected startup sequence
`computeRelations` → `computeBoundaryDistances` → `computeClosestLabels` on a
synthetic 9×8×7 categorical phantom and deterministic fixed memberships. The
initial class-wise boundary-distance map is all zero, matching the
`MultipleToadSegmentation` constructor before the startup calls. The reference
records the source-selected relation tables/maps, boundary distances, closest
label slots, and pass counts. It contains no imaging data.

Generate the selected-method reference with
`scripts/generate_java_toads_initial_relation_state_reference.sh`. It verifies
the `MultipleToadSegmentation.java` SHA-256 before compiling small Java
reference excerpts for the constituent methods. The fixture validates only
this initial preparation state; class variance estimation, atlas shape
sampling, membership updating, and iterative TOADS behavior remain separate.
