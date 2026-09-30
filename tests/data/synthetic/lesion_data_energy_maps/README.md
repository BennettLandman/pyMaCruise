# P39zzzzm volume-level LesionTOADS data-energy fixture

`images_ras.nii.gz` is a 5×6×7×3 synthetic multi-channel volume with true RAS
axis codes and 1.5×2.0×2.5 mm voxel spacing. Intensities are deterministic
binary fractions so the Java and Python fixtures receive identical float32
samples. It was created by `scripts/generate_lesion_data_energy_maps_fixture.py`.

`java_reference.json` contains raw float32 bit patterns for three class-order
cases: OBJECT→LESION→BLACKHOLE, OBJECT→BLACKHOLE→LESION, and
LESION→OBJECT→BLACKHOLE→LESION. The Java excerpt is
`tests/java_reference/LesionDataEnergyMapsReferenceHarness.java`; regenerate
the reference with `scripts/generate_java_lesion_data_energy_maps_reference.sh`.
The selected source is VandyCRUISE `LesionToadSegmentation.java` (SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`), method
`computeMemberships`; wrapper SHA-256 is
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`.

The reference covers only unfielded OBJECT, LESION, and BLACKHOLE raw
residual-energy maps. It excludes Sulcal-CSF special handling, field and
PowerTable weighting, FLAIR weighting, smoothing, priors, inversion,
membership updates, and objective accumulation. Energy maps are zero on the
one-voxel faces skipped by the Java interior scan. LESION and BLACKHOLE share
the lesion-energy slot, so the later class in source order overwrites that
slot; BLACKHOLE alone writes the black-hole slot.
