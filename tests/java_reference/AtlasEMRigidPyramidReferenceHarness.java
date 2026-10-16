import edu.jhmi.rad.medic.libraries.ImageFunctions;

/** Selected DemonToadDeformableAtlas.registerShapesPyramid sampling excerpt. */
public final class AtlasEMRigidPyramidReferenceHarness {
    private static final int X=8,Y=8,Z=8,C=3;

    private static float[][][][] field(int salt) {
        float[][][][] out=new float[X][Y][Z][C];
        for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) for(int c=0;c<C;c++)
            out[x][y][z][c]=(float)((x*17+y*11+z*7+c*13+salt)%97)/32.0f;
        return out;
    }

    private static void bits(StringBuilder out,float[][][][] values) {
        out.append('[');boolean first=true;
        for(int x=0;x<values.length;x++) for(int y=0;y<values[x].length;y++) for(int z=0;z<values[x][y].length;z++) for(int c=0;c<values[x][y][z].length;c++) {
            if(!first)out.append(',');first=false;out.append(Integer.toUnsignedLong(Float.floatToIntBits(values[x][y][z][c])));
        }
        out.append(']');
    }

    public static void main(String[] args) {
        float[][][][] mem=field(3),shape=field(19);
        StringBuilder out=new StringBuilder("{\"demon_atlas_source_sha256\":\"96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a\",\"image_functions_sha256\":\"be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892\",\"levels\":[");
        boolean firstLevel=true;
        for(int factor:new int[]{4,2}) {
            if(!firstLevel)out.append(',');firstLevel=false;
            float[][][][] halfMem=ImageFunctions.subsample(mem,X,Y,Z,C,factor);
            float[][][][] halfShape=ImageFunctions.subsample(shape,X,Y,Z,C,factor);
            out.append("{\"scale\":").append(factor).append(",\"dims\":[").append(halfMem.length).append(',').append(halfMem[0].length).append(',').append(halfMem[0][0].length).append("],\"membership_bits\":");bits(out,halfMem);
            out.append(",\"shape_bits\":");bits(out,halfShape);out.append('}');
        }
        System.out.println(out.append("]}").toString());
    }
}
