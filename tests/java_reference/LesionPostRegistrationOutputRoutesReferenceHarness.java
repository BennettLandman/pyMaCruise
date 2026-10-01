import java.io.DataInputStream;
import java.io.FileInputStream;

/**
 * Selected P39ar terminal-state composition with the class+member and CRUISE
 * output routes from VandyCRUISE AlgorithmLesionToads.runAlgorithm.
 * This consumes referenced terminal arrays; it does not execute a full
 * LesionTOADS workflow.
 */
public final class LesionPostRegistrationOutputRoutesReferenceHarness {
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
    private static int mo(int v,int k) { return v*nc+k; }
    private static boolean csf(String name) { return name.equals("SulcalCSF")||name.equals("Sulcal-CSF"); }
    private static boolean cruiseGm(String name) {
        return name.equals("CorticalGM")||name.equals("Cerebrum-GM")||name.equals("CerebralGM");
    }
    private static boolean excludedCruise(String name) {
        return csf(name)||cruiseGm(name)||name.startsWith("Cerebell")
            ||name.equals("Brainstem")||name.equals("Background");
    }
    private static byte[] membershipClassification() {
        byte[] out=new byte[nx*ny*nz];
        for(int v=0;v<out.length;v++) {
            int best=0; float max=0.0f;
            for(int k=0;k<nc;k++) if(memberships[mo(v,k)]>max) { best=k; max=memberships[mo(v,k)]; }
            float lm=lesion[v];
            out[v]=lm>max?template[lesionClass]:(max>0.0f?template[best]:0);
        }
        return out;
    }
    private static float[] cruiseChannels() {
        int n=nx*ny*nz;
        float[] out=new float[4*n];
        for(int v=0;v<n;v++) {
            for(int k=0;k<nc;k++) {
                String name=NAMES[k];
                if(csf(name))out[v]=memberships[mo(v,k)];
                if(cruiseGm(name))out[n+v]=memberships[mo(v,k)];
            }
            out[2*n+v]=lesion[v]+blackHole[v];
            for(int k=0;k<nc;k++) if(!excludedCruise(NAMES[k])) {
                out[2*n+v]+=memberships[mo(v,k)];
                if((segmentation[v]&0xff)==k)out[3*n+v]=1.0f;
            }
        }
        return out;
    }
    private static void read(String path) throws Exception {
        try(DataInputStream in=new DataInputStream(new FileInputStream(path))) {
            nx=in.readInt();ny=in.readInt();nz=in.readInt();nc=in.readInt();lesionClass=in.readInt();
            outputFromMembership=in.readInt()!=0;includeLesion=in.readInt()!=0;
            int n=nx*ny*nz;
            template=new byte[nc];in.readFully(template);
            segmentation=new byte[n];in.readFully(segmentation);
            classification=new byte[n];in.readFully(classification);
            subtype=new short[n];for(int i=0;i<n;i++)subtype[i]=in.readShort();
            memberships=new float[n*nc];for(int i=0;i<memberships.length;i++)memberships[i]=in.readFloat();
            lesion=new float[n];for(int i=0;i<n;i++)lesion[i]=in.readFloat();
            blackHole=new float[n];for(int i=0;i<n;i++)blackHole[i]=in.readFloat();
        }
        if(nc!=NAMES.length)throw new IllegalArgumentException("unexpected selected atlas channel count");
    }
    private static String ints(byte[] values) {
        StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]);}return s.append(']').toString();
    }
    private static String shorts(short[] values) {
        StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(values[i]);}return s.append(']').toString();
    }
    private static String bits(float[] values) {
        StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.floatToRawIntBits(values[i]));}return s.append(']').toString();
    }
    private static void add(StringBuilder s,String name,String values,boolean[] first) {
        if(!first[0])s.append(',');first[0]=false;
        s.append("{\"name\":\"").append(name).append("\",\"values\":").append(values).append('}');
    }
    public static void main(String[] args) throws Exception {
        String mode=args[1];read(args[0]);
        if(includeLesion)for(int i=0;i<classification.length;i++) {
            classification[i]=(byte)(classification[i]-15*subtype[i]);
            if(classification[i]==26)classification[i]=25;
        }
        byte[] membershipClass=outputFromMembership?membershipClassification():null;
        StringBuilder order=new StringBuilder("[\"cleanUp\",\"exportClassification\",\"exportLesionClassification\"");
        if(outputFromMembership)order.append(",\"updateClassificationFromMemberships\",\"getClassification\"");
        float[] stack=null,cruise=null;
        if(mode.equals("class+member")) {
            order.append(",\"exportMemberships\",\"getLesionMemberships\"");
            int n=nx*ny*nz;stack=new float[(nc+1)*n];
            for(int k=0;k<nc;k++)for(int v=0;v<n;v++)stack[k*n+v]=memberships[mo(v,k)];
            System.arraycopy(lesion,0,stack,nc*n,n);
        } else if(mode.equals("cruise inputs")) {
            order.append(",\"generateCruiseInputs\"");cruise=cruiseChannels();
        } else throw new IllegalArgumentException("unsupported route: "+mode);
        if(includeLesion)order.append(",\"includeLesionnClassification\"");
        order.append(']');
        StringBuilder products=new StringBuilder("[");boolean[] first={true};
        if(stack!=null)add(products,"memberships_plus_lesion",bits(stack),first);
        add(products,"classification",ints(classification),first);
        if(membershipClass!=null)add(products,"classification_from_membership",ints(membershipClass),first);
        add(products,"lesion_subtype",shorts(subtype),first);
        if(cruise!=null)add(products,"cruise_inputs",bits(cruise),first);
        products.append(']');
        System.out.print("{\"order\":"+order+",\"products\":"+products+"}");
    }
}
