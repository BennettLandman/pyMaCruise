import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;

/**
 * Selected VandyCRUISE wrapper composition for P39zzzzs.
 * Source anchors: AlgorithmLesionToads.runAlgorithm() statistics calls then
 * computeMemberships()/Eratio, and LesionToadSegmentation's selected methods.
 * Source SHA-256 values are pinned in the companion fixture JSON.
 */
public final class LesionStatisticsMembershipBridgeReferenceHarness {
    static final int NX=6, NY=5, NZ=4, NC=2, K=4;
    static final float SMOOTH=0.23f, PRIOR_COEF=0.31f, PRIOR_SCALE=0.67f;
    static final float OUTLIER_DISTANCE=2.75f, PREVIOUS_ENERGY=3.25f;
    static final float ZERO=1.0e-10f, INF=1.0e10f, EPS=1.0e-30f;
    static final float[] SMOOTH_PRIOR={0.30f,0.45f,0.60f,0.75f};
    static final float[] ATLAS_PRIOR={0.15f,0.25f,0.35f,0.45f};
    static final boolean[] PRIOR_REL={false,true,false,true};
    static final int[][] DELTA={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    static float[][][][] image=new float[NX][NY][NZ][NC], field=new float[NX][NY][NZ][NC], mem=new float[NX][NY][NZ][K], shape=new float[NX][NY][NZ][K];
    static int[][][] seg=new int[NX][NY][NZ];
    static float[][][] lesion=new float[NX][NY][NZ], bh=new float[NX][NY][NZ], outlier=new float[NX][NY][NZ];
    static float[][] initialCentroid={{0.22f,0.34f,0.46f,0.58f},{0.30f,0.42f,0.54f,0.66f}};

    static int bits(float v){return Float.floatToRawIntBits(v);}
    static float sq(float v){return v*v;}
    static void makeData(){
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            int fold=(x*11+y*7+z*3+x*y+2*y*z)%23;
            for(int c=0;c<NC;c++){
                image[x][y][z][c]=0.15f+((fold+5*c)%23)*0.025f;
                field[x][y][z][c]=0.82f+((x*2+y*3+z*4+c)%5)*0.04f;
            }
            seg[x][y][z]=(x+2*y+z+((x*y+2*z)%3))%K;
            for(int k=0;k<K;k++){
                mem[x][y][z][k]=0.08f+((fold+3*k)%9)*0.045f;
                shape[x][y][z][k]=0.45f+((x*3+y*5+z*7+k*2)%11)*0.035f;
            }
            lesion[x][y][z]=0.04f+((x+2*y+z)%8)*0.06f;
            bh[x][y][z]=0.03f+((2*x+y+3*z)%8)*0.065f;
            outlier[x][y][z]=((x+3*y+z)%5)*0.12f;
        }
    }
    static final class Stats { float[][] centroids,std; float[] lesionCentroid,bhCentroid,lesionStd,bhStd,weights; }
    static Stats refresh(boolean useField){
        Stats s=new Stats(); s.centroids=new float[NC][K];s.std=new float[NC][K];s.lesionCentroid=new float[NC];s.bhCentroid=new float[NC];s.lesionStd=new float[NC];s.bhStd=new float[NC];s.weights=new float[NC];
        for(int c=0;c<NC;c++)for(int k=0;k<K;k++)s.centroids[c][k]=initialCentroid[c][k];
        float[] oldL={0.72f,0.83f},oldB={0.44f,0.55f};
        // computeCentroids: channel, class, x/y/z; direct auxiliary modes are active.
        for(int c=0;c<NC;c++){
            float[] n=new float[K],d=new float[K]; float nl=0,dl=0,nb=0,db=0;
            for(int k=0;k<K;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                float m=mem[x][y][z][k], sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                if(k==0){ n[k]+=m*m*sample; d[k]+=m*m; }
                else if(k==1||k==2){
                    float lm=lesion[x][y][z],bm=bh[x][y][z],winner=Math.max(m,Math.max(lm,bm));
                    if(winner==m){n[k]+=m*m*sample;d[k]+=m*m;}
                    else if(winner==lm){nl+=lm*lm*sample;dl+=lm*lm;}
                    else {nb+=bm*bm*sample;db+=bm*bm;}
                }else {n[k]+=m*m*sample;d[k]+=m*m;}
            }
            for(int k=0;k<K;k++)s.centroids[c][k]=d[k]>0?n[k]/d[k]:0.0f;
            s.lesionCentroid[c]=dl>0?nl/dl:0.0f;s.bhCentroid[c]=db>0?nb/db:0.0f;
        }
        // computeVariances consumes the refreshed means; preserve source k/x/y/z order.
        for(int c=0;c<NC;c++){
            float ln=0,ld=0,bn=0,bd=0;
            for(int k=0;k<K;k++){
                float n=0,d=0;
                for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                    float m=mem[x][y][z][k],lm=lesion[x][y][z],bm=bh[x][y][z];
                    float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                    if((k==1&&m<lm)||(k==2&&m<bm)){
                        if(lm>=bm){float diff=sample-s.lesionCentroid[c],w=lm*lm+EPS;ln+=w*diff*diff;ld+=w;}
                        else {float diff=sample-s.bhCentroid[c],w=lm*lm+EPS;bn+=w*diff*diff;bd+=w;}
                    }else{float diff=sample-s.centroids[c][k],w=m*m+EPS;n+=w*diff*diff;d+=w;}
                }
                s.std[c][k]=d>0?(float)Math.sqrt((double)(n/d)):0.0f;
            }
            s.lesionStd[c]=ld>0?(float)Math.sqrt((double)(ln/ld)):0.0f;s.bhStd[c]=bd>0?(float)Math.sqrt((double)(bn/bd)):0.0f;
        }
        // computeChannelWeights: class, x/y/z; auxiliary source linear residual terms.
        float weightSum=0;
        for(int c=0;c<NC;c++){
            float n=0,d=0;
            for(int k=0;k<K;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
                float m=mem[x][y][z][k],sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                if(k==0){float diff=sample-s.centroids[c][k],w=m*m+EPS;n+=w*diff*diff;d+=w;}
                else if(k==1||k==2){float lm=lesion[x][y][z],bm=bh[x][y][z],dc=sample-s.centroids[c][k],dl=sample-s.lesionCentroid[c],db=sample-s.bhCentroid[c];
                    n+=(m*m+EPS)*dc+(lm*lm+EPS)*dl+(bm*bm+EPS)*db;d+=(m*m+EPS)+(lm*lm+EPS)+(bm*bm+EPS);
                }
            }
            s.weights[c]=n>0?d/n:1.0e10f;weightSum+=s.weights[c];
        }
        for(int c=0;c<NC;c++)s.weights[c]/=weightSum;
        return s;
    }
    static final class Similarity {float[][] cls;float[] lesion;float[] bh;}
    static float sim(float distance){float r=(1.0f-PRIOR_SCALE)/PRIOR_SCALE;return 1.0f/(1.0f+r*distance);}
    static Similarity similarity(Stats s){
        Similarity o=new Similarity();o.cls=new float[K][K];o.lesion=new float[K];o.bh=new float[]{1,1,1,1};
        for(int k=0;k<K;k++)for(int m=0;m<K;m++)if(k!=m){float d=0;for(int c=0;c<NC;c++){float v=s.centroids[c][k]-s.centroids[c][m];d+=s.weights[c]*v*v/(s.std[c][k]+s.std[c][m])/(s.std[c][k]+s.std[c][m]);}o.cls[k][m]=sim(d);}
        for(int k=0;k<K;k++){float d=0;for(int c=0;c<NC;c++){float v=s.lesionCentroid[c]-s.centroids[c][k];d+=s.weights[c]*v*v/(s.lesionStd[c]+s.std[c][k])/(s.lesionStd[c]+s.std[c][k]);}o.lesion[k]=(k==1||k==2)?0.0f:sim(d);}
        return o;
    }
    static final class Energies {float[][][][] cls;float[][][] lesion,bh,outmask;}
    static Energies energies(Stats s,boolean useField){
        Energies o=new Energies();o.cls=new float[NX][NY][NZ][K];o.lesion=new float[NX][NY][NZ];o.bh=new float[NX][NY][NZ];o.outmask=new float[NX][NY][NZ];
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            for(int k=0;k<3;k++){
                float num=0,nl=0,nb=0;
                for(int c=0;c<NC;c++){
                    float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                    float d=sample-s.centroids[c][k];num+=s.weights[c]*d*d;
                    if(k==1){float dl=sample-s.lesionCentroid[c];nl+=s.weights[c]*dl*dl;}
                    if(k==2){float dl=sample-s.lesionCentroid[c],db=sample-s.bhCentroid[c];nl+=s.weights[c]*dl*dl;nb+=s.weights[c]*db*db;}
                }
                o.cls[x][y][z][k]=num;
                if(k==1||k==2)o.lesion[x][y][z]=nl;
                if(k==2)o.bh[x][y][z]=nb;
            }
            for(int c=0;c<NC;c++)if(image[x][y][z][c]!=0.0f){o.outmask[x][y][z]+=s.weights[c];}
            o.cls[x][y][z][3]=o.outmask[x][y][z];
        }
        return o;
    }
    static float invert(float x){return x>ZERO?1.0f/x:INF;}
    static float shapeFactor(float x){return Math.min(1.0f,2.0f*x);}
    static final class Pass {float[][][][] mem,energy;float[][][] lesion,bh,out,den,lesionEnergy,bhEnergy;float objective;}
    static Pass membership(Stats stats,Similarity sim,Energies raw){
        Pass o=new Pass();o.mem=new float[NX][NY][NZ][K];o.energy=new float[NX][NY][NZ][K];o.lesion=new float[NX][NY][NZ];o.bh=new float[NX][NY][NZ];o.out=new float[NX][NY][NZ];o.den=new float[NX][NY][NZ];o.lesionEnergy=new float[NX][NY][NZ];o.bhEnergy=new float[NX][NY][NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){System.arraycopy(mem[x][y][z],0,o.mem[x][y][z],0,K);o.lesion[x][y][z]=lesion[x][y][z];o.bh[x][y][z]=bh[x][y][z];o.out[x][y][z]=outlier[x][y][z];}
        float total=0;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0;
            for(int k=0;k<K;k++){
                float num=raw.cls[x][y][z][k],numL=0,numB=0,numO=0;
                if(k==1||k==2)numL=raw.lesion[x][y][z];if(k==2)numB=raw.bh[x][y][z];if(k==3)numO=OUTLIER_DISTANCE*stats.weights[0]+OUTLIER_DISTANCE*stats.weights[1];
                float e=k==3?(1-o.out[x][y][z])*num+o.out[x][y][z]*numO:num;
                float ngb=0,ngbL=0,ngbB=0,ngbS=0;
                for(int[] d:DELTA){int xx=x+d[0],yy=y+d[1],zz=z+d[2];for(int m=0;m<K;m++)if(m!=k)ngb+=o.mem[xx][yy][zz][m]*o.mem[xx][yy][zz][m];ngbL+=o.lesion[xx][yy][zz]*o.lesion[xx][yy][zz];ngbB+=o.bh[xx][yy][zz]*o.bh[xx][yy][zz];}
                if(k==1||k==2)for(int face:new int[]{1,0,3,2,5,4}){int[] d=DELTA[face];float v=o.mem[x+d[0]][y+d[1]][z+d[2]][k];ngbS+=v*v;}
                float p=SMOOTH_PRIOR[k];num+=p*SMOOTH*(ngb+ngbL+ngbB);e+=0.5f*p*SMOOTH*(ngb+ngbL+ngbB);
                float le=(k==1||k==2)?raw.lesion[x][y][z]:0,bhe=k==2?raw.bh[x][y][z]:0;
                if(k==1||k==2)numL+=p*SMOOTH*(ngb+ngbS+ngbB);if(k==2)numB+=p*SMOOTH*(ngb+ngbS+ngbL);if(k==3)numO+=p*SMOOTH*(ngb+ngbL);
                if(k==1||k==2)le+=0.5f*p*SMOOTH*(ngb+ngbS+ngbB);if(k==2)bhe+=0.5f*p*SMOOTH*(ngb+ngbS+ngbL);
                float inv=0,invL=0,invB=0;
                for(int m=0;m<K;m++){if(m!=k)inv+=ATLAS_PRIOR[m]*sim.cls[k][m]*shape[x][y][z][m]*shape[x][y][z][m];if(k==1||k==2)invL+=ATLAS_PRIOR[m]*sim.lesion[m]*shape[x][y][z][m]*shape[x][y][z][m];if(k==2)invB+=ATLAS_PRIOR[m]*sim.bh[m]*shape[x][y][z][m]*shape[x][y][z][m];}
                inv+=ATLAS_PRIOR[2]*sim.lesion[k]*shape[x][y][z][2]*shape[x][y][z][2];inv+=ATLAS_PRIOR[2]*sim.bh[k]*shape[x][y][z][2]*shape[x][y][z][2];
                num+=PRIOR_COEF*inv;if(k==1||k==2)numL+=PRIOR_COEF*invL;if(k==2)numB+=PRIOR_COEF*invB;if(k==3)numO+=PRIOR_COEF*inv;
                e+=0.5f*PRIOR_COEF*inv;if(k==1||k==2)le+=0.5f*PRIOR_COEF*invL;if(k==2)bhe+=0.5f*PRIOR_COEF*invB;
                o.energy[x][y][z][k]=e;if(k==1||k==2)o.lesionEnergy[x][y][z]=le;if(k==2)o.bhEnergy[x][y][z]=bhe;
                float classMem=invert(num),lesMem=(k==1||k==2)?invert(numL):0,bhMem=k==2?invert(numB):0;
                if(PRIOR_REL[k]){float f=shapeFactor(shape[x][y][z][k]);classMem*=f;if(k==1||k==2)lesMem*=f;if(k==2)bhMem*=f;}
                if(k==3){float c=invert(num),om=invert(numO),f=PRIOR_REL[k]?shapeFactor(shape[x][y][z][k]):1;classMem=f*((1-o.out[x][y][z])*c+o.out[x][y][z]*om);o.out[x][y][z]=f*o.out[x][y][z]*om;}
                o.mem[x][y][z][k]=classMem;if(k==1||k==2)o.lesion[x][y][z]=lesMem;if(k==2)o.bh[x][y][z]=bhMem;
                if(k==1)den+=o.lesion[x][y][z];else if(k==2)den+=o.lesion[x][y][z]+o.bh[x][y][z];den+=o.mem[x][y][z][k];
            }
            o.den[x][y][z]=den;
            if(den>0){for(int k=0;k<K;k++){o.mem[x][y][z][k]/=den;if(k==1||k==2)o.lesion[x][y][z]/=den;if(k==2)o.bh[x][y][z]/=den;else if(k==3)o.out[x][y][z]/=den;}
                for(int k=0;k<K;k++)total+=o.energy[x][y][z][k]*o.mem[x][y][z][k]*o.mem[x][y][z][k];}
            total+=o.lesionEnergy[x][y][z]*o.lesion[x][y][z]*o.lesion[x][y][z]+o.bhEnergy[x][y][z]*o.bh[x][y][z]*o.bh[x][y][z];
        }
        o.objective=total;return o;
    }
    static float ratio(float old,float value){return Math.abs(2.0f*(old-value)/(old+value));}
    static String flat(float[][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(float[] row:a)for(float v:row){if(!f)b.append(',');f=false;b.append(bits(v));}return b.append(']').toString();}
    static String flat(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}return b.append(']').toString();}
    static String flat(float[][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(float[][] p:a)for(float[] q:p)for(float v:q){if(!f)b.append(',');f=false;b.append(bits(v));}return b.append(']').toString();}
    static String flat(float[][][][] a){StringBuilder b=new StringBuilder("[");boolean f=true;for(float[][][] p:a)for(float[][] q:p)for(float[] r:q)for(float v:r){if(!f)b.append(',');f=false;b.append(bits(v));}return b.append(']').toString();}
    public static void main(String[] args){
        makeData();StringBuilder out=new StringBuilder("{");
        for(int fieldMode=0;fieldMode<2;fieldMode++){
            boolean useField=fieldMode==1;Stats s=refresh(useField);Similarity sim=similarity(s);Energies e=energies(s,useField);Pass p=membership(s,sim,e);float er=ratio(PREVIOUS_ENERGY,p.objective);
            if(fieldMode>0)out.append(',');String key=useField?"with_field":"no_field";
            out.append('"').append(key).append("\":{");
            out.append("\"centroids_bits\":").append(flat(s.centroids)).append(',').append("\"lesion_centroid_bits\":").append(flat(s.lesionCentroid)).append(',').append("\"black_hole_centroid_bits\":").append(flat(s.bhCentroid)).append(',');
            out.append("\"class_std_bits\":").append(flat(s.std)).append(',').append("\"lesion_std_bits\":").append(flat(s.lesionStd)).append(',').append("\"black_hole_std_bits\":").append(flat(s.bhStd)).append(',').append("\"weights_bits\":").append(flat(s.weights)).append(',');
            out.append("\"class_similarity_bits\":").append(flat(sim.cls)).append(',').append("\"lesion_similarity_bits\":").append(flat(sim.lesion)).append(',');
            out.append("\"raw_class_energy_bits\":").append(flat(e.cls)).append(',').append("\"raw_lesion_energy_bits\":").append(flat(e.lesion)).append(',').append("\"raw_black_hole_energy_bits\":").append(flat(e.bh)).append(',').append("\"outmask_energy_bits\":").append(flat(e.outmask)).append(',');
            out.append("\"class_memberships_bits\":").append(flat(p.mem)).append(',').append("\"lesion_membership_bits\":").append(flat(p.lesion)).append(',').append("\"black_hole_membership_bits\":").append(flat(p.bh)).append(',').append("\"outlier_membership_bits\":").append(flat(p.out)).append(',').append("\"denominator_bits\":").append(flat(p.den)).append(',').append("\"objective_energy_bits\":").append(flat(p.energy)).append(',').append("\"lesion_objective_energy_bits\":").append(flat(p.lesionEnergy)).append(',').append("\"black_hole_objective_energy_bits\":").append(flat(p.bhEnergy)).append(',');
            out.append("\"objective_bits\":").append(bits(p.objective)).append(',').append("\"energy_ratio_bits\":").append(bits(er)).append('}');
        }
        out.append("}");System.out.print(out);
    }
}
