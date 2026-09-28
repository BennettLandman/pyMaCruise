package javax.vecmath;

public class Vector3f {
    public float x, y, z;
    public Vector3f() {}
    public Vector3f(float x, float y, float z) { this.x = x; this.y = y; this.z = z; }
    public void normalize() {
        float length = (float) Math.sqrt(x * x + y * y + z * z);
        if (length != 0) { x /= length; y /= length; z /= length; }
    }
}
