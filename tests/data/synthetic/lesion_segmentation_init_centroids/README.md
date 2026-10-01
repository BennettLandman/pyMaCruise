# LesionTOADS initial centroid handoff

This fixed 2-channel x 4-class float32 matrix represents the output of the
preprocessing centroid initializer. The method under test copies it into the
segmentation engine's current and prior centroid matrices. These are scalar
class statistics, not a spatial image, so an RAS NIfTI fixture would misstate
their meaning; the fixture stores exact float32 bit patterns in JSON instead.

The source harness pins `LesionToadSegmentation.initCentroids(float[][])`.
The active wrapper calls it at `AlgorithmLesionToads.runAlgorithm` after
constructing the segmentation object and before setting lesion-specific
centroids. This test does not validate centroid estimation or segmentation.
