import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import edu.jhu.ece.iacl.algorithms.tgdm.DynamicPressureTGDM;
import edu.jhu.ece.iacl.algorithms.tgdm.ProfileTGDM;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import gov.nih.mipav.model.structures.ModelImage;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Executes the selected Dura stage sequence on attached-dura and clean-cortex
 * synthetic fields. The sequence and output loop follow
 * MedicAlgorithmRemoveDura09.execute(); the two ProcessingAlgorithm topology
 * calls use the selected AlgorithmTopologyCorrection implementation.
 */
public final class DuraAttachedNoDuraReferenceHarness {
    private static final int NX = Integer.getInteger("dura.nx", 17);
    private static final int NY = Integer.getInteger("dura.ny", 17);
    private static final int NZ = Integer.getInteger("dura.nz", 17);
    private static final double SCALE = Math.min(NX, Math.min(NY, NZ)) / 17.0;
    private static final float MEMBERSHIP = Float.parseFloat(System.getProperty("dura.membership", "0.4"));
    private static final float MASK_TRANSLATION = 2.0f;
    private static final float LEVEL_SET_TRANSLATION = 0.0f;
    private static final int ITERATIONS = Integer.getInteger("dura.iterations", 5);

    static final class Probe extends DynamicPressureTGDM {
        final List<Integer> steps = new ArrayList<>();
        final List<String> passHashes = new ArrayList<>();
        final List<float[][][]> passLevelSets = new ArrayList<>();
        int passSteps;

        Probe() { super(ConnectivityRule.CONNECT_18_6); }

        @Override public void incrementCompletedUnits() {
            super.incrementCompletedUnits();
            passSteps++;
        }

        @Override protected void endNarrowBand(ProfileTGDM profile) {
            super.endNarrowBand(profile);
            steps.add(passSteps);
            passHashes.add(sha256(implicitLevelSet));
            passLevelSets.add(copy(implicitLevelSet));
            passSteps = 0;
        }
    }

    static final class CaseResult {
        String name;
        float[][][] t1, wm, gm, paint, wmTc, union, gmTc;
        float[][][] inner, outer, innerDistance, outerDistance, levelSet;
        float[][][] stripped, wmOut, gmOut, maskLevelSet;
        byte[][][] mask;
        Probe evolution;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("usage: OUTPUT_JSON");
        CaseResult attached = run("attached_dura", true);
        CaseResult clean = run("no_dura", false);
        try (BufferedWriter out = new BufferedWriter(new FileWriter(args[0]))) {
            out.write("{\"source\":\"MedicAlgorithmRemoveDura09 selected execute sequence\",");
            out.write("\"medic_algorithm_remove_dura09_sha256\":\"1ed180daf4a4334c995ad310828ec91cfd5593413eec17ab15404a421d7805a4\",");
            out.write("\"algorithm_topology_correction_sha256\":\"9dd7cd0db8c83e47a21059fe5acc52c7737cb0001222d6019619bf62526159a8\",");
            out.write("\"generic_tgdm_sha256\":\"70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1\",");
            out.write("\"dynamic_pressure_tgdm_sha256\":\"53b872a5828539bb51d70db0c048d1fd8162f6b7a6b63176266ef7e2e796b6fc\",");
            out.write("\"cortex_profile_sha256\":\"f85912d621198de297e65d6b9205505b3ac2de72e67e6840f3968c58a850f505\",");
            out.write("\"distance_field_sha256\":\"92294be2db74351446c3e5898fc35b58497f89ccb38772d4d565885f3f0c35b7\",");
            out.write("\"dimensions_xyz\":["+NX+","+NY+","+NZ+"],\"scale\":"+SCALE+",\"membership_level\":"+MEMBERSHIP+",");
            out.write("\"mask_translation\":"+MASK_TRANSLATION+",\"level_set_translation\":"+LEVEL_SET_TRANSLATION+",\"outer_iterations\":"+ITERATIONS+",");
            out.write("\"cases\":[");
            emitCase(out, attached); out.write(','); emitCase(out, clean);
            out.write("]}");
        }
    }

    private static CaseResult run(String name, boolean attachedDura) {
        CaseResult result = new CaseResult();
        result.name = name;
        result.t1 = new float[NX][NY][NZ];
        result.wm = new float[NX][NY][NZ];
        result.gm = new float[NX][NY][NZ];
        result.paint = new float[NX][NY][NZ];
        BitSet paintBits = new BitSet(NX * NY * NZ);
        for (int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            double cx=(NX-1.0)/2.0, cy=(NY-1.0)/2.0, cz=(NZ-1.0)/2.0;
            double dx=x-cx, dy=y-cy, dz=z-cz;
            double radius=Math.sqrt(dx*dx+dy*dy+dz*dz);
            float wmValue=(float)(255.0*Math.max(0.0, Math.min(1.0, (4.4*SCALE-radius)/(1.7*SCALE))));
            float gmValue=(float)(255.0*Math.max(0.0, Math.min(1.0, 1.0-Math.abs(radius-4.25*SCALE)/(1.8*SCALE))));
            if (attachedDura && Math.round(y)==Math.round(cy) && Math.round(z)==Math.round(cz)
                    && x>=cx+4.0*SCALE && x<=cx+7.0*SCALE) gmValue=190.0f;
            if (attachedDura && x>=cx+6.0*SCALE && x<=cx+7.0*SCALE
                    && y>=cy-SCALE && y<=cy+SCALE && z>=cz-SCALE && z<=cz+SCALE) gmValue=230.0f;
            result.wm[x][y][z]=wmValue;
            result.gm[x][y][z]=gmValue;
            result.t1[x][y][z]=(float)(80.0+2.0*x+3.0*y+5.0*z);
            if (radius<=2.9*SCALE) {
                result.paint[x][y][z]=1.0f;
                paintBits.set(flat(x,y,z));
            }
        }

        return execute(result, paintBits);
    }

    static CaseResult runFromInputs(String name, float[][][] t1, float[][][] wm, float[][][] gm, float[][][] paint) {
        CaseResult result = new CaseResult();
        result.name = name; result.t1 = t1; result.wm = wm; result.gm = gm; result.paint = paint;
        BitSet paintBits = new BitSet(NX * NY * NZ);
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
            if(paint[x][y][z] > 0.0f) paintBits.set(flat(x,y,z));
        return execute(result, paintBits);
    }

    private static CaseResult execute(CaseResult result, BitSet paintBits) {
        result.wmTc=topology(result.wm, paintBits, "wm");
        result.union=new float[NX][NY][NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
            result.union[x][y][z]=result.wmTc[x][y][z]+result.gm[x][y][z];
        result.gmTc=topology(result.union, paintBits, "gm");
        // The wrapper captures these arrays before handing GM to the TGDM.
        result.gmTc=copy(result.gmTc);

        result.inner=new float[NX][NY][NZ];
        result.outer=new float[NX][NY][NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            result.inner[x][y][z]=0.5f-result.wmTc[x][y][z]/255.0f;
            result.outer[x][y][z]=MEMBERSHIP-result.gmTc[x][y][z]/255.0f;
        }
        DistanceField distance = new DistanceField();
        result.innerDistance=distance.solve(new ImageDataFloat(result.inner),15.0).toArray3d();
        result.outerDistance=distance.solve(new ImageDataFloat(result.outer),10.0).toArray3d();
        Probe tgdm=new Probe();
        tgdm.setReferenceLevelSet(new ImageDataFloat(result.innerDistance));
        ImageDataFloat evolved=tgdm.solveBrainMaskSurface(new ImageDataFloat(result.outerDistance),null,
            new ImageDataFloat(copy(result.gmTc)),null,null,0.5f,0.0f,1.0f,0.0f,0.5f,0.5f,0.0f,
            MEMBERSHIP,ITERATIONS);
        result.levelSet=evolved.toArray3d();
        result.evolution=tgdm;
        result.stripped=copy(result.t1);
        // MedicAlgorithmRemoveDura09's output loop masks the original input
        // WM/GM volumes in place; the topology-corrected copies are used for
        // level-set construction and are not the corrected membership outputs.
        result.wmOut=copy(result.wm);
        result.gmOut=copy(result.gm);
        result.maskLevelSet=copy(result.levelSet);
        result.mask=new byte[NX][NY][NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
            if(result.levelSet[x][y][z]>=MASK_TRANSLATION) {
                result.wmOut[x][y][z]=0.0f;
                result.gmOut[x][y][z]=0.0f;
                result.stripped[x][y][z]=0.0f;
                result.mask[x][y][z]=0;
            } else result.mask[x][y][z]=1;
            result.maskLevelSet[x][y][z]-=LEVEL_SET_TRANSLATION;
        }
        return result;
    }

    private static float[][][] topology(float[][][] values, BitSet paint, String name) {
        ModelImage source=new ModelImage(new int[]{NX,NY,NZ},name+"_input");
        ModelImage destination=new ModelImage(new int[]{NX,NY,NZ},name+"_output");
        source.setData(flatten(values));
        AlgorithmTopologyCorrection algorithm=new AlgorithmTopologyCorrection(destination,source,
            "scalar_image",0.05f,3,6,18,"paint_mask",paint,0.0f,1.0f,true,"object->background");
        algorithm.runAlgorithm();
        if(!algorithm.isCompleted()) throw new IllegalStateException("topology call failed: "+name);
        return inflate(destination.getData());
    }

    private static int flat(int x,int y,int z) { return x+NX*y+NX*NY*z; }
    private static float[] flatten(float[][][] v) {
        float[] a=new float[NX*NY*NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) a[flat(x,y,z)]=v[x][y][z];
        return a;
    }
    private static float[][][] inflate(float[] a) {
        float[][][] v=new float[NX][NY][NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) v[x][y][z]=a[flat(x,y,z)];
        return v;
    }
    private static float[][][] copy(float[][][] v) {
        float[][][] c=new float[NX][NY][NZ];
        for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) System.arraycopy(v[x][y],0,c[x][y],0,NZ);
        return c;
    }

    private static void emitCase(BufferedWriter out, CaseResult c) throws Exception {
        out.write("{\"name\":\""+c.name+"\",\"trace_steps\":[");
        for(int i=0;i<c.evolution.steps.size();i++){if(i>0)out.write(',');out.write(Integer.toString(c.evolution.steps.get(i)));}
        out.write("],\"trace_levelset_sha256\":[");
        for(int i=0;i<c.evolution.passHashes.size();i++){if(i>0)out.write(',');out.write("\""+c.evolution.passHashes.get(i)+"\"");}
        out.write("],\"trace_levelsets\":[");
        for(int i=0;i<c.evolution.passLevelSets.size();i++){if(i>0)out.write(',');out.write('{');field(out,"pass_"+i,c.evolution.passLevelSets.get(i));out.write('}');}
        out.write("],\"arrays\":{");
        field(out,"t1",c.t1);out.write(',');field(out,"wm",c.wm);out.write(',');field(out,"gm",c.gm);out.write(',');field(out,"paint",c.paint);out.write(',');
        field(out,"wm_tc",c.wmTc);out.write(',');field(out,"union",c.union);out.write(',');field(out,"gm_tc",c.gmTc);out.write(',');
        field(out,"inner_distance",c.innerDistance);out.write(',');field(out,"outer_distance",c.outerDistance);out.write(',');
        field(out,"level_set",c.levelSet);out.write(',');field(out,"stripped_t1",c.stripped);out.write(',');
        field(out,"corrected_wm",c.wmOut);out.write(',');field(out,"corrected_gm",c.gmOut);out.write(',');
        field(out,"brain_mask_level_set",c.maskLevelSet);out.write(',');byteField(out,"brain_mask",c.mask);
        out.write("}}");
    }
    private static void field(BufferedWriter out,String name,float[][][] v) throws Exception {
        out.write("\""+name+"\":[");boolean first=true;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            if(!first)out.write(',');first=false;out.write(Float.toString(v[x][y][z]));
        }
        out.write(']');
    }
    private static void byteField(BufferedWriter out,String name,byte[][][] v) throws Exception {
        out.write("\""+name+"\":[");boolean first=true;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            if(!first)out.write(',');first=false;out.write(Byte.toString(v[x][y][z]));
        }
        out.write(']');
    }
    private static String sha256(float[][][] values) {
        try {
            ByteBuffer bytes=ByteBuffer.allocate(NX*NY*NZ*4).order(ByteOrder.LITTLE_ENDIAN);
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)bytes.putFloat(values[x][y][z]);
            byte[] digest=MessageDigest.getInstance("SHA-256").digest(bytes.array());
            StringBuilder hex=new StringBuilder();for(byte b:digest)hex.append(String.format("%02x",b));
            return hex.toString();
        } catch(Exception e) { throw new IllegalStateException(e); }
    }
}
