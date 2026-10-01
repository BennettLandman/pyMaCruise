import java.io.DataInputStream;
import java.io.FileInputStream;

/**
 * Selected P39ar terminal-state plus AlgorithmLesionToads dura-output dispatch.
 * It accepts the already branch-referenced terminal arrays and reproduces the
 * post-loop selected wrapper calls/products; it is not a full-workflow harness.
 */
public final class LesionPostRegistrationDuraDispatchReferenceHarness {
    private static final String[] NAMES={
        "Background","Sulcal-CSF","Ventricles","Cerebellum-GM","Cerebrum-GM",
        "Caudate","Thalamus","Putamen","Brainstem","Cerebellum-WM",
        "OuterVentricular-WM","InterVentricular-WM","Frontal-WM"
    };
    private static int nx,ny,nz,nc,lesionClass;
    private static boolean outputFromMembership,includeLesion;
    private static byte[] template,segmentation,classification;
    private static short[] subtype;
    private static float[] memberships,lesion,blackHole;

    private static int voxel(int x,int y,int z) { return (x*ny+y)*nz+z; }
    private static int membershipOffset(int voxel,int k) { return voxel*nc+k; }
    private static boolean isCsf(String name) { return name.equals("SulcalCSF") || name.equals("Sulcal-CSF"); }
    private static boolean isGm(String name) {
        return name.equals("CorticalGM") || name.equals("Cerebrum-GM")
            || name.equals("CerebralGM") || name.equals("Cerebellum-GM")
            || name.equals("CerebellarGM");
    }
    private static boolean excluded(String name) {
        return isCsf(name) || isGm(name) || name.equals("Background");
    }
    private static byte[] membershipClassification() {
        byte[] out=new byte[nx*ny*nz];
        for(int v=0;v<out.length;v++) {
            int best=0; float max=0.0f;
            for(int k=0;k<nc;k++) if(memberships[membershipOffset(v,k)]>max) {
                best=k; max=memberships[membershipOffset(v,k)];
            }
            float lm=lesion[v];
            out[v]=lm>max?template[lesionClass]:(max>0.0f?template[best]:0);
        }
        return out;
    }
    private static float[][] duraChannels() {
        int n=nx*ny*nz;
        float[][] out=new float[4][n];
        for(int k=0;k<nc;k++) {
            if(isCsf(NAMES[k])) for(int v=0;v<n;v++) out[0][v]=memberships[membershipOffset(v,k)];
        }
        for(int v=0;v<n;v++) out[2][v]=lesion[v]+blackHole[v];
        for(int k=0;k<nc;k++) {
            if(!excluded(NAMES[k])) for(int v=0;v<n;v++) {
                out[2][v]+=memberships[membershipOffset(v,k)];
                if((segmentation[v]&0xff)==k)out[3][v]=1.0f;
            }
        }
        for(int k=0;k<nc;k++) if(isGm(NAMES[k])) {
            for(int v=0;v<n;v++)out[1][v]+=memberships[membershipOffset(v,k)];
            for(int v=0;v<n;v++) {
                String selected=NAMES[segmentation[v]&0xff];
                if(selected.equals("Cerebellum-GM") || selected.equals("CerebellarGM)"))out[3][v]=1.0f;
            }
        }
        return out;
    }
    private static void appendInts(StringBuilder s, byte[] values) {
        s.append('['); for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]);} s.append(']');
    }
    private static void appendShorts(StringBuilder s, short[] values) {
        s.append('['); for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]);} s.append(']');
    }
    private static void appendFloatBits(StringBuilder s, float[] values) {
        s.append('['); for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.floatToRawIntBits(values[i]));} s.append(']');
    }
    private static void product(StringBuilder s, String name, String values, boolean[] first) {
        if(!first[0])s.append(','); first[0]=false;
        s.append("{\"name\":\"").append(name).append("\",\"values\":").append(values).append('}');
    }
    private static void read(String path) throws Exception {
        try(DataInputStream in=new DataInputStream(new FileInputStream(path))) {
            nx=in.readInt();ny=in.readInt();nz=in.readInt();nc=in.readInt();lesionClass=in.readInt();
            outputFromMembership=in.readInt()!=0;includeLesion=in.readInt()!=0;
            int n=nx*ny*nz;
            template=new byte[nc]; in.readFully(template);
            segmentation=new byte[n]; in.readFully(segmentation);
            classification=new byte[n]; in.readFully(classification);
            subtype=new short[n]; for(int i=0;i<n;i++)subtype[i]=in.readShort();
            memberships=new float[n*nc]; for(int i=0;i<memberships.length;i++)memberships[i]=in.readFloat();
            lesion=new float[n]; for(int i=0;i<n;i++)lesion[i]=in.readFloat();
            blackHole=new float[n]; for(int i=0;i<n;i++)blackHole[i]=in.readFloat();
        }
    }
    public static void main(String[] args) throws Exception {
        read(args[0]);
        byte[] classificationMem=outputFromMembership?membershipClassification():null;
        float[][] channels=duraChannels();
        if(includeLesion) for(int i=0;i<classification.length;i++) {
            classification[i]=(byte)(classification[i]-15*subtype[i]);
            if(classification[i]==26)classification[i]=25;
        }
        StringBuilder products=new StringBuilder("["); boolean[] first={true}; StringBuilder values=new StringBuilder();
        appendInts(values,classification); product(products,"classification",values.toString(),first);
        if(classificationMem!=null) { values.setLength(0);appendInts(values,classificationMem);product(products,"classification_from_membership",values.toString(),first); }
        values.setLength(0);appendShorts(values,subtype);product(products,"lesion_subtype",values.toString(),first);
        String[] channelNames={"csf","gray_matter","white_matter_fill","white_matter_mask"};
        for(int c=0;c<4;c++){values.setLength(0);appendFloatBits(values,channels[c]);product(products,channelNames[c],values.toString(),first);}
        products.append(']');
        System.out.print("{\"products\":"+products+"}");
    }
}
