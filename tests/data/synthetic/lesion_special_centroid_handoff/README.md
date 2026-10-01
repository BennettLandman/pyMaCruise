# LesionTOADS special centroid setters

This fixture contains exact float32 bit patterns for the lesion and black-hole
centroid vectors supplied by the wrapper after `initCentroids`. The values are
small scalar statistics, not spatial image data, so they are stored as JSON.

The Java harness verifies both setter assignments preserve the input array
references. It also records wrapper call order: lesion centroid first,
black-hole centroid second. The test does not validate how preprocessing
estimates either vector.
