/** Fixed-input OBJECT+LESION excerpt from selected computeExactMemberships. */
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;
import java.util.Locale;

public final class ToadsLesionMembershipPassReferenceHarness {
    private static final String SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final int NX=5, NY=4, NZ=5, C=2, K=2, NIN=16;
    private static final float[][] MOD={{0.75f,0.35f},{1.125f,0.625f}};
    private static final float[][] VAR={{0.8f,1.1f},{1.2f,0.7f}};
    private static final float[] IMGVAR={0.8f,1.25f};
    private static final float[][] CENTERS={{0.25f,0.125f},{-0.125f,-0.25f}}; // modality, OBJECT/LESION
    private static final float[] LESION_CENTER={0.125f,-0.25f};
    private static final float[] SMOOTH_PRIOR={0.6f,0.9f}; // OBJECT, LESION
    private static final float[] ATLAS_PRIOR={0.21f,0.34f};
    private static final float[][] SIM={{0.0f,0.4f},{0.35f,0.0f}};
    private static final int[][] OFF={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

    private static float[][][][] inputs() {
        float[][][][] d=new float[NX][NY][NZ][NIN];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            int code=x+2*y+3*z;
            for(int c=0;c<C;c++) {
                d[x][y][z][c]=((7*x+5*y+3*z+4*c)%17-8)/9.0f;
                d[x][y][z][2+c]=0.82f+((x+2*y+z+c)%5)*0.09f;
                d[x][y][z][4+c]=((2*x+3*y+5*z+c)%9)/31.0f;
            }
            if(x==2&&y==1&&z==2){d[x][y][z][0]=0.125f;d[x][y][z][1]=-0.25f;}
            else if(x==3&&y==1&&z==2){d[x][y][z][0]=0.25f;d[x][y][z][1]=-0.125f;}
            for(int k=0;k<K;k++) {
                d[x][y][z][6+k]=((5*x+7*y+11*z+3*k)%31)/30.0f;
                d[x][y][z][9+k]=0.1f+((x+3*y+2*z+5*k)%11)*0.04f;
            }
            d[x][y][z][8]=((3*x+5*y+7*z)%13)/16.0f;
            d[x][y][z][11]=(code%3==1)?1.0f:0.0f;
            int relation=code%4;
            d[x][y][z][12]=relation==0?0.0f:-1.0f;
            d[x][y][z][13]=relation==1?1.0f:-1.0f;
            d[x][y][z][14]=(relation==1||relation==2)?1.0f:-1.0f;
            d[x][y][z][15]=relation==0?0.0f:-1.0f; // OBJECT appears DIRECT and NEIGHBOR; DIRECT wins.
        }
        return d;
    }

    private static String digest(float[][][][] d) throws Exception {
        MessageDigest md=MessageDigest.getInstance("SHA-256");
        ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<NIN;c++) {
            b.clear();b.putInt(Float.floatToRawIntBits(d[x][y][z][c]));md.update(b.array());
        }
        StringBuilder out=new StringBuilder();for(byte v:md.digest())out.append(String.format(Locale.ROOT,"%02x",v&255));return out.toString();
    }
    private static int lesionIndex(boolean reverse){return reverse?0:1;}
    private static int objectIndex(boolean reverse){return reverse?1:0;}
    private static float[][][][] permuteClasses(float[][][][] d, boolean reverse) {
        float[][][][] a=new float[NX][NY][NZ][K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++)
            a[x][y][z][k]=d[x][y][z][6+(reverse?1-k:k)];
        return a;
    }
    private static float[][][][] permuteShapes(float[][][][] d, boolean reverse) {
        float[][][][] a=new float[NX][NY][NZ][K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<K;k++)
            a[x][y][z][k]=d[x][y][z][9+(reverse?1-k:k)];
        return a;
    }
    private static int[][][][] closest(float[][][][] d, boolean reverse, boolean allDirect) {
        int[][][][] a=new int[NX][NY][NZ][4];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            for(int s=0;s<4;s++)a[x][y][z][s]=(int)d[x][y][z][12+s];
            if(reverse)for(int s=0;s<4;s++)if(a[x][y][z][s]>=0)a[x][y][z][s]=1-a[x][y][z][s];
            if(allDirect){a[x][y][z][0]=0;a[x][y][z][1]=1;a[x][y][z][2]=-1;a[x][y][z][3]=-1;}
        }
        return a;
    }
    private static float factor(int rel,float shape,float priorFactor) {
        if(rel==0)return 1.0f;
        float scaled=(2.0f*priorFactor)*shape;
        if(rel==1)return Math.max(0.25f,Math.min(scaled,1.0f));
        return Math.min(1.0f,scaled);
    }
    private static int relation(int[][][][] closest,int x,int y,int z,int k) {
        if(closest[x][y][z][0]==k||closest[x][y][z][1]==k)return 0;
        if(closest[x][y][z][2]==k||closest[x][y][z][3]==k)return 1;
        return 2;
    }
    private static float[][][][] run(float[][][][] d, boolean reverse, boolean useField, boolean useBias,
            float smoothing,float priorCoeff,boolean allDirect,Result result) {
        int li=lesionIndex(reverse),oi=objectIndex(reverse);
        int[] kind={reverse?1:0,reverse?0:1}; // source class order: 0 OBJECT, 1 LESION, reversed as requested
        float[][][][] mem=permuteClasses(d,reverse),shapes=permuteShapes(d,reverse);
        int[][][][] near=closest(d,reverse,allDirect);
        float[][][][] classTerms=new float[NX][NY][NZ][K],lesionTerms=new float[NX][NY][NZ][K];
        float[][][][] energy=new float[NX][NY][NZ][K];
        float[][][] lesion=new float[NX][NY][NZ],denMap=new float[NX][NY][NZ];
        float[][][][] field=new float[NX][NY][NZ][C],bias=new float[NX][NY][NZ][C];
        float[][] centers=new float[C][K];float[][] varWeight=new float[C][K];float[] smoothPrior=new float[K],atlasPrior=new float[K];float[][] similarity=new float[K][K];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            lesion[x][y][z]=d[x][y][z][8];
            for(int c=0;c<C;c++){field[x][y][z][c]=d[x][y][z][2+c];bias[x][y][z][c]=d[x][y][z][4+c];}
        }
        for(int c=0;c<C;c++)for(int k=0;k<K;k++) {
            int logical=reverse?1-k:k;
            centers[c][k]=CENTERS[c][logical];varWeight[c][k]=VAR[c][logical];
        }
        for(int k=0;k<K;k++) {
            int logical=reverse?1-k:k;smoothPrior[k]=SMOOTH_PRIOR[logical];atlasPrior[k]=ATLAS_PRIOR[logical];
            for(int m=0;m<K;m++)similarity[k][m]=SIM[logical][reverse?1-m:m];
        }
        float[][][][] priorShapes=priorCoeff>0?shapes:new float[NX][NY][NZ][K];
        float[][][][] priorNum=new float[NX][NY][NZ][K],priorObjective=new float[NX][NY][NZ][K];
        if(priorCoeff>0)for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++)for(int k=0;k<K;k++) {
            float inv=0.0f;for(int m=0;m<K;m++)if(m!=k){float term=atlasPrior[m]*similarity[k][m];term=term*priorShapes[x][y][z][m];term=term*priorShapes[x][y][z][m];inv=inv+term;}
            priorNum[x][y][z][k]=priorCoeff*inv;priorObjective[x][y][z][k]=(0.5f*priorCoeff)*inv;
        }
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
            int seg=((int)d[x][y][z][11]==1)?li:oi;
            for(int k=0;k<K;k++) {
                float num=0.0f,numL=0.0f;
                for(int c=0;c<C;c++) {
                    float corrected=useField?field[x][y][z][c]*d[x][y][z][c]:d[x][y][z][c];
                    if(useBias)corrected=corrected-bias[x][y][z][c];
                    float residual=corrected-centers[c][k];
                    float weight=MOD[c][0]*varWeight[c][k];weight=weight*IMGVAR[c];
                    float contribution=weight*residual;contribution=contribution*residual;num=num+contribution;
                    if(kind[k]==1&&seg==k) {
                        float lesionResidual=corrected-LESION_CENTER[c];
                        float lesionWeight=MOD[c][1]*varWeight[c][k];lesionWeight=lesionWeight*IMGVAR[c];
                        float lesionContribution=lesionWeight*lesionResidual;lesionContribution=lesionContribution*lesionResidual;numL=numL+lesionContribution;
                    }
                }
                classTerms[x][y][z][k]=num;if(kind[k]==1)lesionTerms[x][y][z][k]=numL;
            }
        }
        float total=0.0f;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++) {
            float den=0.0f;float[] localObjective=new float[K];int seg=((int)d[x][y][z][11]==1)?li:oi;
            for(int k=0;k<K;k++) {
                int logical=reverse?1-k:k;boolean assigned=kind[k]==1&&seg==k;
                float num=classTerms[x][y][z][k],numL=lesionTerms[x][y][z][k];
                float e=assigned?(1.0f-lesion[x][y][z])*num+lesion[x][y][z]*numL:num;
                if(smoothing>0.0f) {
                    float ngb=0.0f;for(int[] off:OFF)for(int m=0;m<K;m++)if(m!=k){float v=mem[x+off[0]][y+off[1]][z+off[2]][m];ngb=ngb+v*v;}
                    num=num+smoothPrior[k]*smoothing*ngb;
                    if(kind[k]==1&&seg==k)numL=numL+smoothPrior[k]*smoothing*ngb;
                    e=e+0.5f*smoothPrior[k]*smoothing*ngb;
                }
                if(priorCoeff>0.0f){num=num+priorNum[x][y][z][k];if(assigned)numL=numL+priorNum[x][y][z][k];e=e+priorObjective[x][y][z][k];}
                energy[x][y][z][k]=e;localObjective[k]=e;
                float rec=num>1e-10f?1.0f/num:1e10f,recL=assigned?(numL>1e-10f?1.0f/numL:1e10f):0.0f;
                int rel=relation(near,x,y,z,k);float f=factor(rel,priorShapes[x][y][z][k],0.57f);
                if(kind[k]==1) {
                    if(assigned){float mixed=(1.0f-lesion[x][y][z])*rec+lesion[x][y][z]*recL;mem[x][y][z][k]=f*mixed;lesion[x][y][z]=(f*lesion[x][y][z])*recL;}
                    else {mem[x][y][z][k]=(rel==0||rel==1)?rec:f*rec;lesion[x][y][z]=0.0f;}
                } else mem[x][y][z][k]=f*rec;
                den=den+mem[x][y][z][k];
            }
            denMap[x][y][z]=den;
            if(den>0.0f)for(int k=0;k<K;k++) {
                mem[x][y][z][k]=mem[x][y][z][k]/den;
                if(kind[k]==1&&seg==k)lesion[x][y][z]=lesion[x][y][z]/den;
                total=total+localObjective[k]*mem[x][y][z][k]*mem[x][y][z][k];
            }
        }
        result.classTerms=classTerms;result.lesionTerms=lesionTerms;result.energy=energy;result.memberships=mem;result.lesion=lesion;result.den=denMap;result.objective=total;return mem;
    }
    private static final class Result {float[][][][] classTerms,lesionTerms,energy,memberships;float[][][] lesion,den;float objective;}
    private static int[] bits(float[][][][] a){int[] b=new int[NX*NY*NZ*a[0][0][0].length];int q=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(float v:a[x][y][z])b[q++]=Float.floatToRawIntBits(v);return b;}
    private static int[] bits(float[][][] a){int[] b=new int[NX*NY*NZ];int q=0;for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)b[q++]=Float.floatToRawIntBits(a[x][y][z]);return b;}
    private static String json(int[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}return s.append(']').toString();}
    private static void emit(StringBuilder out,String name,Result r,boolean comma){if(comma)out.append(',');out.append("{\"name\":\"").append(name).append("\",\"class_data_terms_bits\":").append(json(bits(r.classTerms))).append(",\"lesion_data_terms_bits\":").append(json(bits(r.lesionTerms))).append(",\"class_energy_bits\":").append(json(bits(r.energy))).append(",\"memberships_bits\":").append(json(bits(r.memberships))).append(",\"lesion_membership_bits\":").append(json(bits(r.lesion))).append(",\"denominator_bits\":").append(json(bits(r.den))).append(",\"objective_bits\":").append(Float.floatToRawIntBits(r.objective)).append('}');}
    public static void main(String[] args) throws Exception {
        float[][][][] d=inputs();StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"input_sha256\":\"").append(digest(d)).append("\",\"shape\":[5,4,5],\"channels\":16,\"cases\":[");
        Result r=new Result();run(d,false,true,true,0.18f,0.31f,false,r);emit(out,"object_lesion_field_bias_smooth_prior",r,false);
        r=new Result();run(d,true,false,false,0.0f,0.0f,true,r);emit(out,"lesion_object_direct_no_field_no_bias",r,true);
        r=new Result();run(d,false,true,false,0.14f,0.0f,true,r);emit(out,"object_lesion_field_smoothing_only",r,true);
        r=new Result();run(d,true,false,true,0.0f,0.22f,false,r);emit(out,"lesion_object_bias_prior_only",r,true);
        out.append("]}");System.out.print(out.toString());
    }
}
