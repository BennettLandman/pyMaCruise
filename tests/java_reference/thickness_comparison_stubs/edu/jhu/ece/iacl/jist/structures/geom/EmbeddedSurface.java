package edu.jhu.ece.iacl.jist.structures.geom;

import javax.vecmath.Point3f;

/** Minimal coordinate/index/point-data shell for selected thickness source. */
public final class EmbeddedSurface {
    private final Point3f[] points;
    private final int[] indices;
    private double[][] vertexData;
    private String name = "synthetic";

    public EmbeddedSurface(Point3f[] points, int[] indices) {
        this.points = points;
        this.indices = indices;
    }

    public int getIndexCount() { return indices.length; }
    public int getCoordinateIndex(int index) { return indices[index]; }
    public Point3f getVertex(int index) { return points[index]; }
    public int getVertexCount() { return points.length; }
    public void setVertexData(double[][] values) { vertexData = values; }
    public double getVertexDataAtOffset(int vertex, int offset) { return vertexData[vertex][offset]; }
    public double[] getVertexData(int vertex) { return vertexData[vertex]; }
    public void setName(String value) { name = value; }
    public String getName() { return name; }
}
