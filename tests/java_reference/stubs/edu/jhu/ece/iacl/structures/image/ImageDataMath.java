package edu.jhu.ece.iacl.structures.image;

import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

public final class ImageDataMath {
    private ImageDataMath() {}

    public static ImageDataFloat vectorMagnitude4D(ImageDataFloat vector) {
        int x = vector.getRows(), y = vector.getCols(), z = vector.getSlices();
        float[][][][] values = vector.toArray4d();
        ImageDataFloat magnitude = new ImageDataFloat(x, y, z);
        float[][][] output = magnitude.toArray3d();
        for (int i = 0; i < x; i++)
            for (int j = 0; j < y; j++)
                for (int k = 0; k < z; k++) {
                    double sum = 0;
                    for (int c = 0; c < vector.getComponents(); c++) sum += Math.pow(values[i][j][k][c], 2.0);
                    output[i][j][k] = (float) Math.sqrt(sum);
                }
        return magnitude;
    }

    public static void addFloatImage(ImageDataFloat target, ImageDataFloat increment) {
        float[][][] a = target.toArray3d(), b = increment.toArray3d();
        for (int i = 0; i < target.getRows(); i++)
            for (int j = 0; j < target.getCols(); j++)
                for (int k = 0; k < target.getSlices(); k++) a[i][j][k] += b[i][j][k];
    }
}
