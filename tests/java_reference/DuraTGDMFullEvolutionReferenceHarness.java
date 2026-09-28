import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import edu.jhu.ece.iacl.algorithms.tgdm.DynamicPressureTGDM;
import edu.jhu.ece.iacl.algorithms.tgdm.ProfileTGDM;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;

/** Executes the selected VandyCRUISE Dura TGDM path on a tiny synthetic field. */
public final class DuraTGDMFullEvolutionReferenceHarness {
    private static final int NX = Integer.getInteger("dura.nx", 11);
    private static final int NY = Integer.getInteger("dura.ny", 10);
    private static final int NZ = Integer.getInteger("dura.nz", 9);
    private static final boolean FOLDED = Boolean.getBoolean("dura.folded");

    private static final class Probe extends DynamicPressureTGDM {
        final List<float[][][]> passes = new ArrayList<>();
        final List<float[][][]> preReinitialization = new ArrayList<>();
        final List<Integer> steps = new ArrayList<>();
        final List<Integer> minePasses = new ArrayList<>();
        int passSteps = 0;
        int mineHits = 0;
        Probe() { super(ConnectivityRule.CONNECT_18_6); }

        float[][][] initialPointUpdates(float[][][] phi, float[][][] gm) {
            rows=NX; cols=NY; slices=NZ;
            implicitLevelSet=phi;
            regionForces=new float[NX][NY][NZ];
            float denominator=0.4f*255.0f;
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++)
                regionForces[x][y][z]=(gm[x][y][z]-denominator)/denominator;
            setAlpha(0.5f); setAlignmentForce(1.0f); curvatureBiasForce=0.0f;
            timeStep=0.1;
            ProfileTGDM profile=new edu.jhu.ece.iacl.algorithms.tgdm.DynamicPressureCortexProfile(
                    -0.5,0.0,0.5,5);
            float[][][] updated=new float[NX][NY][NZ];
            for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) {
                NBPoint point=new NBPoint(); point.x=x; point.y=y; point.z=z;
                updated[x][y][z]=(float)updatePoint(point,profile);
            }
            return updated;
        }

        @Override public void incrementCompletedUnits() {
            super.incrementCompletedUnits();
            passSteps++;
        }

        @Override protected void updateNarrowBand(ProfileTGDM profile) {
            super.updateNarrowBand(profile);
            if (mineHit) mineHits++;
        }

        @Override protected void endNarrowBand(ProfileTGDM profile) {
            float[][][] before = new float[NX][NY][NZ];
            for (int x=0;x<NX;x++) for (int y=0;y<NY;y++)
                System.arraycopy(implicitLevelSet[x][y], 0, before[x][y], 0, NZ);
            preReinitialization.add(before);
            super.endNarrowBand(profile);
            float[][][] copy = new float[NX][NY][NZ];
            for (int x=0;x<NX;x++) for (int y=0;y<NY;y++)
                System.arraycopy(implicitLevelSet[x][y], 0, copy[x][y], 0, NZ);
            passes.add(copy);
            steps.add(passSteps);
            minePasses.add(mineHits);
            passSteps = 0;
            mineHits = 0;
        }
    }

    private static float[][][] field(int which) {
        float[][][] values = new float[NX][NY][NZ];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            if (FOLDED) {
                double cx=(NX-1.0)/2.0, cy=(NY-1.0)/2.0, cz=(NZ-1.0)/2.0;
                double dx=x-cx, dy=y-cy, dz=z-cz;
                double radius=Math.sqrt(dx*dx+dy*dy+dz*dz);
                double theta=Math.atan2(dy,dx);
                double polar=(radius == 0.0) ? 0.0 : Math.acos(dz/radius);
                double folds=1.15*Math.sin(7.0*theta)*Math.sin(4.0*polar)
                        +0.45*Math.cos(5.0*theta+0.7)*Math.sin(3.0*polar);
                if (which == 0) values[x][y][z]=(float)(radius-(8.0+folds));
                else if (which == 1) values[x][y][z]=(float)(radius-(4.8+0.52*folds));
                else values[x][y][z]=(float)(118.0+42.0*Math.sin(0.24*x+0.16*y)
                        +35.0*Math.cos(0.19*y-0.21*z));
                continue;
            }
            double dx=x-5.0, dy=y-4.5, dz=z-4.0;
            if (which == 0) {
                values[x][y][z]=(float)(Math.sqrt(dx*dx+dy*dy+dz*dz)-3.2);
            } else if (which == 1) {
                double rx=x-5.25, ry=y-4.25, rz=z-4.0;
                values[x][y][z]=(float)(Math.sqrt(rx*rx+ry*ry+rz*rz)-2.6);
            } else {
                values[x][y][z]=(float)(102.0 + 31.0*Math.sin(0.38*x+0.20*y)
                        + 27.0*Math.cos(0.31*y-0.25*z));
            }
        }
        return values;
    }

    private static void emitVolume(StringBuilder out, float[][][] volume) {
        out.append('[');
        boolean first=true;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            if(!first) out.append(','); first=false;
            out.append(String.format("%.9g", volume[x][y][z]));
        }
        out.append(']');
    }

    public static void main(String[] args) {
        PrintStream stdout=System.out;
        System.setOut(new PrintStream(new ByteArrayOutputStream()));
        float[][][] outer=field(0), inner=field(1), gm=field(2);
        Probe tgdm=new Probe();
        tgdm.setReferenceLevelSet(new ImageDataFloat(inner));
        float[][][] initialUpdates=tgdm.initialPointUpdates(outer,gm);
        ImageDataFloat result=tgdm.solveBrainMaskSurface(new ImageDataFloat(outer),null,
                new ImageDataFloat(gm),null,null,0.5f,0.0f,1.0f,0.0f,0.5f,0.5f,
                0.0f,0.4f,5);
        System.setOut(stdout);

        StringBuilder json=new StringBuilder();
        json.append("{\"source\":\"VandyCRUISE GenericTGDM + DynamicPressureTGDM + DynamicPressureCortexProfile; local DistanceFieldyh seam\",");
        json.append("\"phantom\":\""+(FOLDED ? "folded_cortex" : "nested_surfaces")+"\",");
        json.append("\"generic_sha256\":\"70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1\",");
        json.append("\"dims\":["+NX+","+NY+","+NZ+"],\"membership_level\":0.4,\"alpha\":0.5,");
        json.append("\"alignment_force\":1.0,\"curvature_force\":0.5,\"pressure_force\":0.5,");
        json.append("\"outer_iterations\":5,\"pass_steps\":[");
        for(int i=0;i<tgdm.steps.size();i++){if(i>0)json.append(',');json.append(tgdm.steps.get(i));}
        json.append("],\"mine_hit_steps\":[");
        for(int i=0;i<tgdm.minePasses.size();i++){if(i>0)json.append(',');json.append(tgdm.minePasses.get(i));}
        json.append("],\"pre_reinitialization_passes\":[");
        for(int i=0;i<tgdm.preReinitialization.size();i++){if(i>0)json.append(',');emitVolume(json,tgdm.preReinitialization.get(i));}
        json.append("],\"passes\":[");
        for(int i=0;i<tgdm.passes.size();i++){if(i>0)json.append(',');emitVolume(json,tgdm.passes.get(i));}
        json.append("],\"initial_point_update\":"); emitVolume(json,initialUpdates);
        json.append(",\"final\":"); emitVolume(json,result.toArray3d());
        json.append('}');
        System.out.println(json);
    }
}
