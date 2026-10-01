/* Selected VandyCRUISE LesionToadPreprocess.initialMemberships() excerpt.
 * Source SHA-256: 72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08.
 * Isolates zero, nonzero sub-threshold, and above-threshold squared distances.
 */
public final class LesionInitialMembershipSmallDistanceReferenceHarness {
    static final int NX=4, NY=4, NZ=4, NC=2, K=3, PAD=1;
    static final int X0=5, XN=7, Y0=8, YN=10, Z0=2, ZN=4;
    static final float ZERO=1.0e-20f, INF=1.0e20f;
    static float image(int c,int x,int y,int z) {
        if(x<PAD || x>=NX-PAD || y<PAD || y>=NY-PAD || z<PAD || z>=NZ-PAD) return 0.0f;
        int pattern=(x+2*y+3*z)%3;
        if(pattern==1) return c<2?1.0e-12f:0.0f;
        if(pattern==2 && c==0) return 5.0e-11f;
        return 0.0f;
    }
    static float[][] centroids() {
        return new float[][] {{0.0f,1.0e-12f,1.0e-9f},{0.0f,1.0e-12f,0.0f}};
    }
    static int bits(float value) { return Float.floatToRawIntBits(value); }
    static String run() {
        float[][] means=centroids();
        int[] output=new int[NX*NY*NZ*K];
        int[] centroidBits=new int[NC*K];
        for(int c=0;c<NC;c++) for(int k=0;k<K;k++) centroidBits[c*K+k]=bits(means[c][k]);
        for(int sx=X0;sx<XN;sx++) for(int sy=Y0;sy<YN;sy++) for(int sz=Z0;sz<ZN;sz++) {
            int x=sx-X0+PAD,y=sy-Y0+PAD,z=sz-Z0+PAD;
            float den=0.0f;
            float[] nums=new float[K],grouping=new float[K];
            for(int k=0;k<K;k++) {
                float num=0.0f;
                for(int c=0;c<NC;c++) { float d=image(c,x,y,z)-means[c][k]; num+=d*d; }
                if(num>ZERO) num=1.0f/num; else num=INF;
                nums[k]=num; den+=num; grouping[k]=0.0f;
                for(int m=0;m<K;m++) if(means[0][m]==means[0][k]) grouping[k]++;
            }
            if(den>0.0f) for(int k=0;k<K;k++) output[((x*NY+y)*NZ+z)*K+k]=bits(grouping[k]*nums[k]/den);
        }
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08\",\"shape\":["+NX+","+NY+","+NZ+"],\"centroid_bits\":[");
        for(int i=0;i<centroidBits.length;i++){if(i>0)out.append(',');out.append(centroidBits[i]);}
        out.append("],\"membership_bits\":[");
        for(int i=0;i<output.length;i++){if(i>0)out.append(',');out.append(output[i]);}
        return out.append("]}").toString();
    }
    public static void main(String[] args) { System.out.println(run()); }
}
