/*
 * Test-only transcription of the supplied ImageFunctions.subsample(float[][][][])
 * method. Its nested block-offset loops, factor-cubed float divisor, and
 * truncation to complete blocks mirror the VandyCRUISE source.
 */
public final class ToadsPyramidSubsampleHarness {
    private static final int NX=9, NY=8, NZ=7, NC=3;

    private static float value(int x,int y,int z,int k) {
        return ((11*x+5*y+3*z+7*k)%29)/28.0f;
    }

    private static float[][][][] subsample(int factor) {
        int nsx=NX/factor,nsy=NY/factor,nsz=NZ/factor;
        float[][][][] output=new float[nsx][nsy][nsz][NC];
        float scale=factor*factor*factor;
        for(int x=0;x<nsx;x++) for(int y=0;y<nsy;y++) for(int z=0;z<nsz;z++)
            for(int k=0;k<NC;k++) {
                output[x][y][z][k]=0.0f;
                for(int i=0;i<factor;i++) for(int j=0;j<factor;j++) for(int l=0;l<factor;l++)
                    output[x][y][z][k]+=value(x*factor+i,y*factor+j,z*factor+l,k);
                output[x][y][z][k]/=scale;
            }
        return output;
    }

    private static void append(int factor,StringBuilder out) {
        float[][][][] values=subsample(factor);
        if(out.length()>1) out.append(',');
        out.append('"').append(factor).append("\":{\"shape\":[")
            .append(values.length).append(',').append(values[0].length).append(',')
            .append(values[0][0].length).append(',').append(NC).append("],\"values\":[");
        boolean first=true;
        for(int x=0;x<values.length;x++) for(int y=0;y<values[0].length;y++)
            for(int z=0;z<values[0][0].length;z++) for(int k=0;k<NC;k++) {
                if(!first) out.append(',');
                out.append(Float.toString(values[x][y][z][k]));
                first=false;
            }
        out.append("]}");
    }

    public static void main(String[] args) {
        StringBuilder out=new StringBuilder("{");
        append(2,out);
        append(4,out);
        System.out.println(out.append('}'));
    }
}
