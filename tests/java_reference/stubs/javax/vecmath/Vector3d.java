package javax.vecmath;

public class Vector3d {
    public double x, y, z;
    public Vector3d() {}
    public Vector3d(double x, double y, double z) { this.x=x; this.y=y; this.z=z; }
    public Vector3d(Vector3d other) { this(other.x, other.y, other.z); }
    public void scale(double value) { x*=value; y*=value; z*=value; }
    public void add(Vector3d other) { x+=other.x; y+=other.y; z+=other.z; }
    public double dot(Vector3d other) { return x*other.x + y*other.y + z*other.z; }
    public double length() { return Math.sqrt(x*x + y*y + z*z); }
    public void normalize() { double value=length(); if(value>0) scale(1.0/value); }
}
