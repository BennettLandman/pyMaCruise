/*
 * Selected wrapper call-order excerpt from VandyCRUISE AlgorithmLesionToads
 * and LesionToadAtlas. AlgorithmLesionToads.java SHA-256:
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 * LesionToadAtlas.java SHA-256:
 * 718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e.
 * ParametricTransform.java SHA-256:
 * 22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac.
 * The harness records init -> preprocessing transform seed -> same-mode update
 * -> pyramid order. P39o separately validates the pyramid optimizer numerics.
 */
public final class LesionInitialAtlasAlignmentReferenceHarness {
    static final int NX=33, NY=33, NZ=33, LEVELS=4, INIT_ALIGN_ITER=50;
    static final int NC=2, CLASSES=3, PAD=1, X0=5, Y0=7, Z0=9, XN=36, YN=38, ZN=40;
    static final float COUPLING=1.0f;
    static final float[] PREPROCESS_TRANSFORM={0.008f,-0.012f,0.006f,0.08f,-0.06f,0.05f};
    static final String[] CALLS={
        "preprocess.initialMemberships",
        "atlas.initShapeRegistration(memberships,\"rigid\",initAlignIter,levels,1.0f)",
        "atlas.setTransform(preprocess.getTransform())",
        "atlas.updateShapeRegistration(memberships,\"rigid\",initAlignIter,levels,1.0f)",
        "atlas.registerShapesPyramid"
    };
    static int[] bits(float[] values) {
        int[] result=new int[values.length];
        for(int i=0;i<values.length;i++) result[i]=Float.floatToIntBits(values[i]);
        return result;
    }
    static void append(StringBuilder out,int[] values) {
        out.append('[');
        for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(values[i]);}
        out.append(']');
    }
    static String jsonQuoted(String value) { return "\""+value.replace("\\","\\\\").replace("\"","\\\"")+"\""; }
    static float image(int c,int x,int y,int z) {
        if(x<PAD || x>=NX-PAD || y<PAD || y>=NY-PAD || z<PAD || z>=NZ-PAD)return 0.0f;
        int fold=(x+2*y+3*z+c)%7-3;
        return (5*x-3*y+7*z+9*c+fold)/512.0f;
    }
    static float[][] centroids() {
        int center=NX/2;
        return new float[][]{{image(0,center,center,center),image(0,center,center,center),-0.25f},
            {image(1,center,center,center),0.35f,-0.15f}};
    }
    static int[] membershipBits() {
        float[][] means=centroids();
        int[] values=new int[NX*NY*NZ*CLASSES];
        for(int sx=X0;sx<XN;sx++)for(int sy=Y0;sy<YN;sy++)for(int sz=Z0;sz<ZN;sz++) {
            int x=sx-X0+PAD,y=sy-Y0+PAD,z=sz-Z0+PAD;
            float den=0.0f;float[] numerator=new float[CLASSES],grouping=new float[CLASSES];
            for(int k=0;k<CLASSES;k++) {
                float num=0.0f;
                for(int c=0;c<NC;c++){float d=image(c,x,y,z)-means[c][k];num+=d*d;}
                if(num>1.0e-20f)num=1.0f/num;else num=1.0e20f;
                numerator[k]=num;den+=num;grouping[k]=0.0f;
                for(int m=0;m<CLASSES;m++)if(means[0][m]==means[0][k])grouping[k]++;
            }
            if(den>0.0f)for(int k=0;k<CLASSES;k++) {
                values[((x*NY+y)*NZ+z)*CLASSES+k]=Float.floatToIntBits(grouping[k]*numerator[k]/den);
            }
        }
        return values;
    }
    public static void main(String[] args) {
        float[] transform=new float[6];
        float[] imageCenter={NX/2.0f,NY/2.0f,NZ/2.0f};
        int subsample=2;
        // initShapeRegistration(memberships, "rigid", INIT_ALIGN_ITER, LEVELS, COUPLING)
        // setTransform(preprocess.getTransform())
        System.arraycopy(PREPROCESS_TRANSFORM,0,transform,0,6);
        // updateShapeRegistration with the same "rigid" type calls
        // changeTransformType(rigid, rigid), which returns the same parameters.
        int[] scales=new int[LEVELS-1];
        int scale=1;
        for(int i=1;i<LEVELS;i++)scale*=2;
        for(int l=LEVELS;l>1;l--){scales[LEVELS-l]=scale;scale/=2;}
        StringBuilder out=new StringBuilder("{\"source\":\"selected LesionTOADS initial membership and atlas alignment call-order excerpt\",\"algorithm_source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"atlas_source_sha256\":\"718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e\",\"transform_source_sha256\":\"22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac\",\"dimensions\":["+NX+","+NY+","+NZ+"],\"levels\":"+LEVELS+",\"init_align_iterations\":"+INIT_ALIGN_ITER+",\"transform_mode\":\"rigid\",\"coupling_factor_bits\":"+Float.floatToIntBits(COUPLING)+",\"subsample\":"+subsample+",\"preprocessing_transform_bits\":");
        append(out,bits(PREPROCESS_TRANSFORM));
        out.append(",\"post_update_transform_bits\":");append(out,bits(transform));
        out.append(",\"image_center_bits\":");append(out,bits(imageCenter));
        out.append(",\"pyramid_scales\":");append(out,scales);
        out.append(",\"membership_bits\":");append(out,membershipBits());
        out.append(",\"call_order\":[");
        for(int i=0;i<CALLS.length;i++){if(i>0)out.append(',');out.append(jsonQuoted(CALLS[i]));}
        System.out.println(out.append("]}").toString());
    }
}
