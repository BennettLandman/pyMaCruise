# LesionTOADS OUTMASK membership reference

Twelve selected-Java scalar cases capture the OUTMASK assignment after the
class and outmask energy terms have already been inverted. Direct relations
blend the class and outmask memberships using the existing outlier
membership, then update the outlier map. Prior relations additionally apply
the source shape factor to the mixed class value and to the outlier value
before multiplying by the outmask membership.

The cases cover outlier weights at 0, 0.25, 0.5, 0.75, and 1; direct and
prior relation paths; shape factors below, at, and above saturation; and
float32-sensitive inputs. Denominator normalization is a separate operation.
Regenerate with `scripts/generate_java_lesion_toads_outmask_membership_reference.sh`.
