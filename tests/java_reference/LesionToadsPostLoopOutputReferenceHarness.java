/**
 * Selected in-memory output tail of AlgorithmLesionToads.runAlgorithm.
 * Source: VandyCRUISE AlgorithmLesionToads.java
 * SHA-256: 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb
 * The harness intentionally excludes cleanup, JIST image creation, UI and I/O.
 */
public final class LesionToadsPostLoopOutputReferenceHarness {
    private static final int NX=2, NY=1, NZ=1, NC=3;
    private static final byte[][][] BASE={{ {26} }, { {26} }};
    private static final short[][][] SUBTYPE={{ {1} }, { {0} }};
    private static final byte[] TEMPLATE={10,20,26};
    private static final int[] SEGMENTATION={0,2};
    private static final String[] NAMES={"SulcalCSF","CorticalGM","Other"};
    private static final float[][][][] MEM={
        {{{0.8f,0.1f,0.1f}}}, {{{0.1f,0.2f,0.7f}}}
    };
    private static final float[][][] LESION={{{0.05f}},{{0.10f}}};
    private static final float[][][] BLACK_HOLE={{{0.02f}},{{0.03f}}};

    private static String quote(String value) { return "\""+value+"\""; }
    private static String ints(byte[][][] a) {
        StringBuilder s=new StringBuilder("["); int n=0;
        for (int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            if(n++>0)s.append(','); s.append(a[x][y][z]);
        }
        return s.append(']').toString();
    }
    private static String shorts(short[][][] a) {
        StringBuilder s=new StringBuilder("["); int n=0;
        for (int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            if(n++>0)s.append(','); s.append(a[x][y][z]);
        }
        return s.append(']').toString();
    }
    private static String floats(float[][][][] a, int channels) {
        StringBuilder s=new StringBuilder("["); int n=0;
        for (int k=0;k<channels;k++) for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            if(n++>0)s.append(','); s.append(Float.floatToRawIntBits(a[k][x][y][z]));
        }
        return s.append(']').toString();
    }
    private static boolean isCsf(String name) {
        return name.equals("SulcalCSF") || name.equals("Sulcal-CSF");
    }
    private static boolean isCruiseGm(String name) {
        return name.equals("CorticalGM") || name.equals("Cerebrum-GM") || name.equals("CerebralGM");
    }
    private static boolean isDuraGm(String name) {
        return isCruiseGm(name) || name.equals("Cerebellum-GM") || name.equals("CerebellarGM");
    }
    private static boolean excludedCruise(String name) {
        return isCsf(name) || isCruiseGm(name) || name.startsWith("Cerebell")
            || name.equals("Brainstem") || name.equals("Background");
    }
    private static boolean excludedDura(String name) {
        return isCsf(name) || isDuraGm(name) || name.equals("Background");
    }
    private static float[][][][] routedChannels(String mode) {
        float[][][][] out=new float[4][NX][NY][NZ];
        for(int x=0;x<NX;x++) {
            if(mode.equals("cruise inputs")) {
                for(int k=0;k<NC;k++) {
                    if(isCsf(NAMES[k]))out[0][x][0][0]=MEM[x][0][0][k];
                    if(isCruiseGm(NAMES[k]))out[1][x][0][0]=MEM[x][0][0][k];
                }
                out[2][x][0][0]=LESION[x][0][0]+BLACK_HOLE[x][0][0];
                for(int k=0;k<NC;k++)if(!excludedCruise(NAMES[k])) {
                    out[2][x][0][0]+=MEM[x][0][0][k];
                    if(SEGMENTATION[x]==k)out[3][x][0][0]=1.0f;
                }
            } else {
                out[2][x][0][0]=LESION[x][0][0]+BLACK_HOLE[x][0][0];
                for(int k=0;k<NC;k++) {
                    if(isCsf(NAMES[k]))out[0][x][0][0]=MEM[x][0][0][k];
                    if(!excludedDura(NAMES[k])) {
                        out[2][x][0][0]+=MEM[x][0][0][k];
                        if(SEGMENTATION[x]==k)out[3][x][0][0]=1.0f;
                    }
                    if(isDuraGm(NAMES[k])) {
                        out[1][x][0][0]+=MEM[x][0][0][k];
                        String selected=NAMES[SEGMENTATION[x]];
                        if(selected.equals("Cerebellum-GM") || selected.equals("CerebellarGM)"))
                            out[3][x][0][0]=1.0f;
                    }
                }
            }
        }
        return out;
    }
    private static void addProduct(StringBuilder out, String name, String values, boolean[] first) {
        if(!first[0])out.append(','); first[0]=false;
        out.append("{\"name\":").append(quote(name)).append(",\"values\":").append(values).append('}');
    }
    private static String order(String mode, boolean fromMembership, boolean include) {
        StringBuilder s=new StringBuilder("[\"cleanUp\",\"exportClassification\",\"exportLesionClassification\"");
        if(fromMembership)s.append(",\"updateClassificationFromMemberships\",\"getClassification\"");
        if(mode.equals("class+member"))s.append(",\"exportMemberships\",\"getLesionMemberships\"");
        else if(mode.equals("cruise inputs"))s.append(",\"generateCruiseInputs\"");
        else if(mode.equals("dura removal inputs"))s.append(",\"generateDuraRemovalOutputs\"");
        if(include)s.append(",\"includeLesionnClassification\"");
        s.append(']'); return s.toString();
    }
    public static void main(String[] args) {
        String mode=args.length>0?args[0]:"class+member";
        boolean fromMembership=args.length>1 && Boolean.parseBoolean(args[1]);
        boolean include=args.length>2 && Boolean.parseBoolean(args[2]);
        byte[][][] classification=new byte[NX][NY][NZ];
        short[][][] subtype=new short[NX][NY][NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            classification[x][y][z]=BASE[x][y][z]; subtype[x][y][z]=SUBTYPE[x][y][z];
        }
        byte[][][] classificationMem=null;
        if(fromMembership) {
            classificationMem=new byte[NX][NY][NZ];
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                int best=0; float max=0.0f;
                for(int k=0;k<NC;k++) if(MEM[x][y][z][k]>max) {best=k;max=MEM[x][y][z][k];}
                classificationMem[x][y][z]=max>0.0f?TEMPLATE[best]:0;
            }
        }
        float[][][][] stack=null;
        if(mode.equals("class+member")) {
            stack=new float[NC+1][NX][NY][NZ];
            for(int k=0;k<NC;k++) for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) stack[k][x][y][z]=MEM[x][y][z][k];
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) stack[NC][x][y][z]=LESION[x][y][z];
        }
        if(include) for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            classification[x][y][z]-=15*subtype[x][y][z];
            if(classification[x][y][z]==26) classification[x][y][z]=25;
        }
        StringBuilder products=new StringBuilder("["); boolean[] first={true};
        if(stack!=null)addProduct(products,"memberships_plus_lesion",floats(stack,NC+1),first);
        addProduct(products,"classification",ints(classification),first);
        if(classificationMem!=null)addProduct(products,"classification_from_membership",ints(classificationMem),first);
        addProduct(products,"lesion_subtype",shorts(subtype),first);
        if(mode.equals("cruise inputs"))addProduct(products,"cruise_inputs",floats(routedChannels(mode),4),first);
        if(mode.equals("dura removal inputs"))addProduct(products,"dura_removal_inputs",floats(routedChannels(mode),4),first);
        products.append(']');
        System.out.print("{\"order\":"+order(mode,fromMembership,include)+",\"classification\":"+ints(classification)+",\"subtype\":"+shorts(subtype)+",\"classification_mem\":"+(classificationMem==null?"null":ints(classificationMem))+",\"stack_bits\":"+(stack==null?"null":floats(stack,NC+1))+",\"products\":"+products+"}");
    }
}
