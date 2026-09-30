/*
 * Selected VandyCRUISE LesionToadPreprocess.initialMemberships() excerpt.
 * Source SHA-256: 72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08.
 * It preserves class/channel accumulation, 1e-20/1e20 handling, repeated
 * first-channel centroid grouping, normalization, and the padded active crop.
 * Preprocessing and atlas registration are intentionally outside the excerpt.
 */
public final class LesionInitialMembershipsReferenceHarness {
    static final int NX=7, NY=8, NZ=9, NC=2, K=4, PAD=1;
    static final int X0=2, XN=7, Y0=3, YN=9, Z0=4, ZN=11;
    static final float ZERO=1.0e-20f, INF=1.0e20f;

    static float image(int c,int x,int y,int z) {
        if(x<PAD || x>=NX-PAD || y<PAD || y>=NY-PAD || z<PAD || z>=NZ-PAD) return 0.0f;
        int fold=(x+2*y+3*z+c)%5-2;
        return (7*x-5*y+3*z+11*c+fold)/64.0f;
    }
    static float[][] centroid() {
        return new float[][] {
            {image(0,1,1,1), image(0,1,1,1), -0.25f, 0.125f},
            {image(1,1,1,1), 0.375f, 0.4f, -0.375f}
        };
    }
    static int bits(float v) { return Float.floatToIntBits(v); }
    static void appendBits(StringBuilder b,int[] values) {
        b.append('['); for(int i=0;i<values.length;i++){if(i>0)b.append(','); b.append(values[i]);} b.append(']');
    }
    static String run() {
        float[][] means=centroid();
        int[] output=new int[NX*NY*NZ*K];
        for(int sx=X0;sx<XN;sx++) for(int sy=Y0;sy<YN;sy++) for(int sz=Z0;sz<ZN;sz++) {
            int x=sx-X0+PAD, y=sy-Y0+PAD, z=sz-Z0+PAD;
            float den=0.0f;
            float[] nums=new float[K], grouping=new float[K];
            for(int k=0;k<K;k++) {
                float num=0.0f;
                for(int c=0;c<NC;c++) {
                    float d=image(c,x,y,z)-means[c][k];
                    num += d*d;
                }
                if(num>ZERO) num=1.0f/num; else num=INF;
                nums[k]=num;
                den += num;
                grouping[k]=0.0f;
                for(int m=0;m<K;m++) if(means[0][m]==means[0][k]) grouping[k]++;
            }
            if(den>0.0f) for(int k=0;k<K;k++) {
                output[((x*NY+y)*NZ+z)*K+k]=bits(grouping[k]*nums[k]/den);
            }
        }
        int[] centroidBits=new int[NC*K];
        for(int c=0;c<NC;c++) for(int k=0;k<K;k++) centroidBits[c*K+k]=bits(means[c][k]);
        StringBuilder b=new StringBuilder("{\"source_sha256\":\"72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08\",\"shape\":["+NX+","+NY+","+NZ+"],\"channels\":"+NC+",\"classes\":"+K+",\"pad_border\":"+PAD+",\"crop_bounds\":[["+X0+","+XN+"],["+Y0+","+YN+"],["+Z0+","+ZN+"]],\"centroid_bits\":");
        appendBits(b,centroidBits); b.append(",\"membership_bits\":"); appendBits(b,output); return b.append('}').toString();
    }
    public static void main(String[] args) { System.out.println(run()); }
}
