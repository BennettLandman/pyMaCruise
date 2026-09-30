import java.util.ArrayList;
import java.util.List;

/** Selected-method excerpt of LesionToadSegmentation.propagateThinning. */
public final class LesionToadsThinningTransitionReferenceHarness {
    private static final int C = 3;
    private static final String DIRECT_ACCEPT = "212221122120200120112112212";
    private static final String DIRECT_REJECT = "221001200120102002101100002";
    private static final String FALLBACK_ACCEPT = "222222101122202022120200001";

    private static final class Scenario {
        final String name;
        final short[][][] seg = new short[3][3][3];
        final short[][][] skeleton = new short[3][3][3];
        final float[] mem = {0.125f, 0.25f, 0.0625f};
        final float[][][] lesion = new float[3][3][3];
        final float[][][] blackHole = new float[3][3][3];
        final float[][][][] shape = new float[3][3][3][C];
        final boolean[] available = {true, true, true};
        final float[] volume = {12.5f, 17.25f, 9.0f};
        final float[] minimum = {2.0f, 2.0f, 2.0f};
        final int k;
        final boolean useMax;

        Scenario(String name, String encoded, int k, boolean useMax) {
            this.name=name; this.k=k; this.useMax=useMax;
            for (int n=0;n<27;n++) {
                int x=n/9, y=(n/3)%3, z=n%3;
                seg[x][y][z]=(short)(encoded.charAt(n)-'0');
                skeleton[x][y][z]=seg[x][y][z];
                lesion[x][y][z]=0.125f;
                blackHole[x][y][z]=0.0625f;
                for (int m=0;m<C;m++) shape[x][y][z][m]=0.0f;
            }
            skeleton[1][1][1]=(short)k;
            if (name.startsWith("direct") || name.startsWith("fallback")) {
                // A positive shape prior is the fallback criterion, without
                // a relation-matrix lookup in this selected Java class.
                for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++)
                    if (i*i+j*j+l*l==1) shape[1+i][1+j][1+l][1]=0.375f;
            }
            if (name.equals("critical_volume")) volume[k]=minimum[k];
            if (name.equals("best_already")) mem[0]=0.75f;
            if (name.equals("lesion_volume_max")) {
                mem[1]=0.125f; mem[0]=0.0f; mem[2]=0.0f;
                lesion[1][1][1]=0.5f; blackHole[1][1][1]=0.375f;
            }
        }
    }

    private static float lesionScore(float ordinary, float lesion, float black, boolean useMax) {
        if (useMax) {
            float inner=lesion>black?lesion:black;
            return ordinary>inner?ordinary:inner;
        }
        return (ordinary+lesion)+black;
    }

    private static int best(Scenario s) {
        int id=0; float max=0.0f;
        for (int m=0;m<C;m++) {
            float score=(m==1||m==2)
                ? lesionScore(s.mem[m],s.lesion[1][1][1],s.blackHole[1][1][1],s.useMax)
                : s.mem[m];
            if (score>max) { id=m; max=s.mem[m]; } // source's retained ordinary-membership quirk
        }
        return id;
    }

    private static boolean topology(short[][][] labels, int candidate) {
        // The same selected LUT-backed predicate already checked against the
        // method excerpt in ToadsHomeomorphicChangeReferenceHarness.
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(labels, (short)1, (short)1,
            (short)1, (short)candidate);
    }

    private static String bits(float value) { return Integer.toString(Float.floatToRawIntBits(value)); }
    private static void emit(Scenario s, boolean first) {
        if (!first) System.out.print(",");
        int id=best(s), candidate=-1; boolean change=false; Boolean accepted=null;
        if (id!=s.k) {
            float prior=0.0f; int next=-1; boolean direct=false;
            for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if(i*i+j*j+l*l==1) {
                short neighbor=s.seg[1+i][1+j][1+l];
                if (neighbor==id) direct=true;
                else if (neighbor!=s.k && s.shape[1+i][1+j][1+l][id]>prior) {
                    next=neighbor; prior=s.shape[1+i][1+j][1+l][id];
                }
            }
            if (direct) {candidate=id; change=true;}
            else if (next!=-1) {candidate=next; change=true;}
            if (change) {accepted=topology(s.seg,candidate); change=accepted;}
        }
        boolean processed=s.seg[1][1][1]==s.k && s.volume[s.k]>s.minimum[s.k];
        short label=s.seg[1][1][1], skel=s.skeleton[1][1][1];
        float[] volumes=s.volume.clone(); boolean[] available=s.available.clone();
        float order=0.375f;
        if (processed) {
            available[s.k]=false; order=-0.625f;
            float decrement=s.k==1
                ? (s.useMax?lesionScore(s.mem[s.k],s.lesion[1][1][1],s.blackHole[1][1][1],true)
                    :(s.mem[s.k]+s.lesion[1][1][1])+s.blackHole[1][1][1])
                : s.mem[s.k];
            volumes[s.k]-=decrement; skel=-1;
            if(change) label=(short)candidate;
        }
        System.out.printf("{\"name\":\"%s\",\"best\":%d,\"candidate\":%d,\"topology\":%s,\"processed\":%s,\"seg\":%d,\"skeleton\":%d,\"volume_bits\":[%s,%s,%s],\"available\":[%s,%s,%s],\"ordering_bits\":%s}",
            s.name,id,candidate,accepted==null?"null":accepted.toString(),processed,label,skel,
            bits(volumes[0]),bits(volumes[1]),bits(volumes[2]),available[0],available[1],available[2],bits(order));
    }

    public static void main(String[] args) {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        List<Scenario> cases=new ArrayList<>();
        cases.add(new Scenario("direct_accept",DIRECT_ACCEPT,0,false));
        cases.add(new Scenario("direct_reject",DIRECT_REJECT,0,false));
        cases.add(new Scenario("fallback_accept",FALLBACK_ACCEPT,0,false));
        cases.add(new Scenario("no_connected_target","000000000000000000000000000",0,false));
        cases.add(new Scenario("best_already","000000000000000000000000000",0,false));
        cases.add(new Scenario("critical_volume","000000000000000000000000000",0,false));
        cases.add(new Scenario("lesion_volume_max","111111111111111111111111111",1,true));
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.propagateThinning selected event excerpt\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"source_lines\":\"2124-2198\",\"cases\":[");
        for (int i=0;i<cases.size();i++) emit(cases.get(i),i==0);
        System.out.println("]}");
    }
}
