# LesionTOADS OUTMASK normalization integration

Three selected-Java one-voxel cases compose direct/prior OUTMASK mixing and
outlier-map updates with the shared denominator, class membership
normalization, repeated outlier-map division, and local energy accumulation.
Cases cover direct and prior relations, zero/unit initial outlier weight,
and two sequential OUTMASK classes sharing one outlier map.

The class/outmask values are already inverted inputs. Regenerate the fixture
with `scripts/generate_java_lesion_toads_outmask_normalization_reference.sh`.
