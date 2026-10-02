# TOADS GM-derived crop phantom

Asymmetric 13×12×11 RAS synthetic maps on a 0.8×1.2×2.0 mm grid. The GM crop mask has an irregular positive region spanning x=3..9, y=2..8, z=3..7 with one internal hole; WM-fill stores coordinate-coded float32 values, and WM-mask stores a distinct binary pattern. With edge width 2, source-style inclusive bounds are x=1..11, y=0..10, z=1..9, giving output dimensions 11×11×9. This fixture checks one shared GM-derived crop and target voxel order.
