import edu.jhmi.rad.medic.utilities.PowerTable;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

/** P18czo+: selected two-channel continuation with supplied fields and nonzero transformed shapes. */
public final class ToadsP18czoNonzeroShapesNextGrowthReferenceHarness extends ToadsDefault11ContinueThinningReferenceHarness {
    private static final int NC=2;
    private static final float SMOOTH=.043f, FACTOR=.64f, COEFFICIENT=.16f, PRIOR_SCALE=.05f;
    private static final float[] MODALITY={.73f,.46f};
    private static final float[] SMOOTH_PRIOR={.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_PRIOR={.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] VARIANCE_WEIGHT={
        {.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f},
        {1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}
    };
    private static final float[][] EXPONENT=new float[NC][11];
    private static final int[] CATEGORY={-1,0,0,1,1,1,1,1,2,2,2};
    private final float[][][][] images2=new float[N][N][N][NC], fields2=new float[N][N][N][NC];
    private final float[][] centers=new float[NC][C], deviations=new float[NC][C];
    private final float[][] similarity=new float[C][C];
    private final float[][][][] energy=new float[N][N][N][C];
    private final float[][][] denominator=new float[N][N][N];
    private final PowerTable[][] powers=new PowerTable[3][NC];
    private final float[] imageVariance=new float[NC], imageStddev=new float[NC];
    private float objective=0.0f;
    private int directBranchCount,neighborBranchCount,shapeOnlyBranchCount,smoothingPriorCount,inverseShapePriorCount;
    private final BinaryHeap4D nextGrowQueue=new BinaryHeap4D(N*N*3,BinaryHeap4D.MINTREE);
    private final boolean[][][][] nextGrowAvailable=new boolean[N][N][N][C];
    private final java.util.List<String> nextGrowthEvents=new java.util.ArrayList<>();
    private float[] nextVolumes;

    private static float secondImage(int x,int y,int z) {
        int label=(x*3+y*5+z*7+(x*y)%C)%C;
        for(int k=0;k<C;k++)if(x==1+k&&y==1+(k*3)%12&&z==1+(k*5)%12)label=k;
        int delta=(x*5+y*11+z*7)%19-9;
        return (11-label)/12.0f+delta/1900.0f;
    }
    private static int[] closest(int x,int y,int z) {
        return new int[]{(x+y+z)%C,(x+2*y+3*z+1)%C,(x+3*y+2*z+4)%C,(2*x+y+z+7)%C};
    }
    private void loadTwoChannels(Path fieldPath) throws IOException {
        byte[] bytes=Files.readAllBytes(fieldPath);
        int count=N*N*N*NC;
        if(bytes.length!=count*Float.BYTES)throw new IllegalArgumentException("expected 14x14x14x2 little-endian float32 fields");
        ByteBuffer input=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
            float f0=input.getFloat(),f1=input.getFloat();
            field[x][y][z][0]=f0;
            fields2[x][y][z][0]=f0;
            fields2[x][y][z][1]=f1;
            images2[x][y][z][0]=images[x][y][z][0];
            images2[x][y][z][1]=secondImage(x,y,z);
        }
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)EXPONENT[c][k]=CATEGORY[k]<0?1.0f:(CATEGORY[k]==0?.8f:1.0f);
    }
    private void loadNonzeroShapes(Path shapePath) throws IOException {
        byte[] bytes=Files.readAllBytes(shapePath);
        int count=N*N*N*C;
        if(bytes.length!=count*Float.BYTES)throw new IllegalArgumentException("expected 14x14x14x11 little-endian float32 transformed shapes");
        ByteBuffer input=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)for(int k=0;k<C;k++)
            shape[x][y][z][k]=input.getFloat();
    }
    private void captureThinningVolumes() throws Exception {
        java.lang.reflect.Field field=ToadsDefault11ContinueThinningReferenceHarness.class.getDeclaredField("volume");
        field.setAccessible(true);nextVolumes=((float[])field.get(this)).clone();
    }
    private void refreshMultichannelStatistics() {
        PowerTable[] tables={new PowerTable(0,4,.000001f,.8f),new PowerTable(0,4,.000001f,1.0f),new PowerTable(0,4,.000001f,1.0f)};
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++) {
            float old=c==0?oldCentroids[k]:(11-k)/12.0f, num=0.0f, den=0.0f;
            for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c], mem=mems[x][y][z][k], weight=mem*mem;
                if(k>0) {
                    float residual=corrected-old,distance=residual*residual,factor=1.0f;
                    if(distance>1e-10f)factor=(float)(tables[CATEGORY[k]].lookup((double)distance,(double)EXPONENT[c][k])/distance);
                    weight=factor*mem;weight*=mem;
                }
                num+=weight*corrected;den+=weight;
            }
            centers[c][k]=den>0.0f?.9f*num/den+.1f*old:old;
        }
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++) {
            float num=0.0f,den=0.0f;
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++)if(seg[x][y][z]==k) {
                float mem=mems[x][y][z][k],weight=mem*mem+1e-30f;
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c],diff=corrected-centers[c][k];
                num+=weight*diff*diff;den+=weight;
            }
            deviations[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
        }
    }
    private void computeImageVariance() {
        float sum=0.0f;
        for(int c=0;c<NC;c++) {
            float num=0.0f,den=0.0f;
            // P18's selected default atlas class zero is MASK. The source
            // computeExactImageVariance switch contributes optimized classes
            // only; the MASK branch contributes neither numerator nor weight.
            for(int k=1;k<C;k++)for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
                float mem=mems[x][y][z][k],weight=mem*mem+1e-30f;
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c];
                float residual=corrected-centers[c][k],diff=residual*residual;
                num+=weight*diff*diff;
                den+=weight;
            }
            imageVariance[c]=num>0.0f?den/num:1.0e10f;
            sum+=imageVariance[c];
        }
        for(int c=0;c<NC;c++) {
            imageStddev[c]=(float)Math.sqrt(imageVariance[c]);
            imageVariance[c]=imageVariance[c]/sum;
        }
        float scale=(1.0f-PRIOR_SCALE)/PRIOR_SCALE;
        for(int k=0;k<C;k++)for(int m=0;m<C;m++)if(k!=m) {
            float distance=0.0f;
            for(int c=0;c<NC;c++) {
                float d=centers[c][k]-centers[c][m],sumStd=deviations[c][k]+deviations[c][m];
                float term=imageVariance[c]*d;term*=d;term/=sumStd;term/=sumStd;distance+=term;
            }
            similarity[k][m]=1.0f/(1.0f+scale*distance);
        }
    }
    private void updateMemberships() {
        for(int type=0;type<3;type++)for(int c=0;c<NC;c++)powers[type][c]=new PowerTable(0,4,.000001f,type==0?.8f:1.0f);
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)for(int k=0;k<C;k++) {
            if(CATEGORY[k]<0) {
                for(int c=0;c<NC;c++){float w=MODALITY[c]*VARIANCE_WEIGHT[c][k];w*=imageVariance[c];float v=images2[x][y][z][c];if(v!=0.0f)energy[x][y][z][k]+=w;}
            } else for(int c=0;c<NC;c++) {
                float corrected=fields2[x][y][z][c]*images2[x][y][z][c],residual=corrected-centers[c][k],d=residual*residual;
                double term=powers[CATEGORY[k]][c].lookup((double)d,(double)EXPONENT[c][k]);
                float w=MODALITY[c]*VARIANCE_WEIGHT[c][k];w*=imageVariance[c];
                energy[x][y][z][k]=(float)(energy[x][y][z][k]+w*term);
            }
        }
        int[][] offsets={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++) {
            int[] labels=closest(x,y,z);float den=0.0f;float[] objectiveAt=new float[C];
            for(int k=0;k<C;k++) {
                boolean direct=labels[0]==k||labels[1]==k,neighbor=labels[2]==k||labels[3]==k;
                if(!direct&&!neighbor&&shape[x][y][z][k]==0.0f){energy[x][y][z][k]=0.0f;mems[x][y][z][k]=0.0f;continue;}
                if(direct)directBranchCount++;else if(neighbor)neighborBranchCount++;else shapeOnlyBranchCount++;
                float num=energy[x][y][z][k],obj=num,neighborSum=0.0f;
                for(int[] d:offsets)for(int m=0;m<C;m++)if(m!=k){float v=mems[x+d[0]][y+d[1]][z+d[2]][m];neighborSum+=v*v;}
                if(neighborSum>0.0f)smoothingPriorCount++;
                float increment=SMOOTH_PRIOR[k]*SMOOTH*neighborSum;num+=increment;obj+=.5f*SMOOTH_PRIOR[k]*SMOOTH*neighborSum;
                float inverse=0.0f;
                for(int m=0;m<C;m++)if(m!=k){float v=ATLAS_PRIOR[m]*similarity[k][m];v*=shape[x][y][z][m];v*=shape[x][y][z][m];inverse+=v;}
                if(inverse>0.0f)inverseShapePriorCount++;
                num+=COEFFICIENT*inverse;obj+=.5f*COEFFICIENT*inverse;
                float factor=shape[x][y][z][k];
                if(direct)factor=1.0f;else if(neighbor)factor=Math.max(.25f,Math.min(2.0f*FACTOR*factor,1.0f));else factor=Math.min(1.0f,2.0f*FACTOR*factor);
                float value=factor*(num>1e-10f?1.0f/num:1e10f);mems[x][y][z][k]=value;den+=value;objectiveAt[k]=obj;
            }
            denominator[x][y][z]=den;
            if(den>0.0f)for(int k=0;k<C;k++){mems[x][y][z][k]/=den;float contribution=objectiveAt[k]*mems[x][y][z][k];objective+=contribution*mems[x][y][z][k];}
        }
    }
    private int bestNextMembership(int x,int y,int z) {
        int best=0;float maximum=0.0f;
        for(int k=0;k<C;k++)if(mems[x][y][z][k]>maximum){best=k;maximum=mems[x][y][z][k];}
        return best;
    }
    private float nextMembershipRatio(int x,int y,int z,int label) {
        int best=bestNextMembership(x,y,z);
        return label==best?1.0f:(mems[x][y][z][best]+1e-10f)/(mems[x][y][z][label]+1e-10f);
    }
    private boolean homeomorphicNextGrow(int x,int y,int z,int label) {
        short[][][] local=new short[N][N][N];
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int l=0;l<N;l++)local[i][j][l]=seg[i][j][l];
        return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(local,(short)x,(short)y,(short)z,(short)label);
    }
    private void propagateNextGrowing() {
        nextGrowQueue.reset();nextGrowQueue.setMinTree();
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)
            for(int k=0;k<C;k++)nextGrowAvailable[x][y][z][k]=skel[x][y][z]==EMPTY;
        int[][] face={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
        for(int x=1;x<N-1;x++)for(int y=1;y<N-1;y++)for(int z=1;z<N-1;z++)for(int k=0;k<C;k++) {
            if(!nextGrowAvailable[x][y][z][k])continue;
            boolean boundary=false;for(int[] d:face)if(skel[x+d[0]][y+d[1]][z+d[2]]==k){boundary=true;break;}
            if(!boundary)continue;
            nextGrowAvailable[x][y][z][k]=false;
            nextGrowQueue.addValue(ordering[x][y][z]+nextMembershipRatio(x,y,z,k),x,y,z,k);
        }
        while(nextGrowQueue.isNotEmpty()) {
            float value=nextGrowQueue.getFirst();int x=nextGrowQueue.getFirstX(),y=nextGrowQueue.getFirstY(),z=nextGrowQueue.getFirstZ(),k=nextGrowQueue.getFirstK();
            nextGrowQueue.removeFirst();int accepted=-1;java.util.List<String> inserted=new java.util.ArrayList<>();
            if(skel[x][y][z]==EMPTY) {
                boolean regular=homeomorphicNextGrow(x,y,z,k);accepted=regular?1:0;
                if(regular) {
                    skel[x][y][z]=(byte)k;for(int m=0;m<C;m++)nextGrowAvailable[x][y][z][m]=false;
                    ordering[x][y][z]=value;nextVolumes[k]+=mems[x][y][z][k];seg[x][y][z]=(byte)k;
                    for(int[] d:face) {
                        int xi=x+d[0],yj=y+d[1],zl=z+d[2];if(!nextGrowAvailable[xi][yj][zl][k])continue;
                        float next=value+nextMembershipRatio(xi,yj,zl,k);nextGrowAvailable[xi][yj][zl][k]=false;
                        nextGrowQueue.addValue(next,xi,yj,zl,k);
                        inserted.add("["+Float.floatToIntBits(next)+","+xi+","+yj+","+zl+","+k+"]");
                    }
                } else nextGrowAvailable[x][y][z][k]=true;
            }
            nextGrowthEvents.add("["+Float.floatToIntBits(value)+","+x+","+y+","+z+","+k+","+accepted+",["+String.join(",",inserted)+"]]");
        }
    }
    private static String hash(float[][][][] a)throws Exception {MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[][][] x:a)for(float[][] y:x)for(float[] z:y)for(float v:z){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hash(float[][][] a)throws Exception {MessageDigest d=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);for(float[][] x:a)for(float[] y:x)for(float v:y){b.clear();b.putInt(Float.floatToIntBits(v));d.update(b.array());}return hex(d.digest());}
    private static String hex(byte[] a){StringBuilder s=new StringBuilder();for(byte b:a)s.append(String.format("%02x",b));return s.toString();}
    private static String labels(byte[][][] a){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(v<10?(char)('0'+v):(char)('a'+v-10));return b.toString();}
    private static String bits(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.floatToIntBits(a[i]));}return b.append(']').toString();}
    public static void main(String[] args) {
        if(args.length!=2)throw new IllegalArgumentException("usage: harness two-channel-fields.bin transformed-shapes.bin");
        try {
            ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
            ToadsP18czoNonzeroShapesNextGrowthReferenceHarness h=new ToadsP18czoNonzeroShapesNextGrowthReferenceHarness();
            h.initialize();h.grow();h.refreshRelationAndInitializeThinning();h.thin();h.captureThinningVolumes();h.loadTwoChannels(Path.of(args[0]));h.refreshMultichannelStatistics();h.computeImageVariance();h.loadNonzeroShapes(Path.of(args[1]));h.updateMemberships();
            String preSeg=labels(h.seg),preSkel=gridLabels(h.skel,true),preOrdering=gridBits(h.ordering),preVolumes=bits(h.nextVolumes);
            h.propagateNextGrowing();
            System.out.print("{\"reference_scope\":\"selected Java formula harness; fixed source-connected grow/thin state, dynamic two-channel exact statistics/image variances, nonzero transformed-shape membership then source-ordered next propagateGrowing; queue/topology continuation only\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"wrapper_source_sha256\":\"b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d\",\"power_table_source_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"shape\":["+N+","+N+","+N+","+C+"],\"source_order\":[\"computeExactCentroids\",\"computeExactVariances\",\"computeExactImageVariance\",\"computeExactMemberships\",\"convergenceCheck.continue\",\"propagateGrowing.initGrowingLabels\",\"propagateGrowing.queue\",\"propagateGrowing.topology\"],\"pre_growth_segmentation\":\""+preSeg+"\",\"pre_growth_skeleton\":\""+preSkel+"\",\"pre_growth_ordering_bits\":"+preOrdering+",\"pre_growth_volume_bits\":"+preVolumes+",\"image_variance_bits\":"+bits(h.imageVariance)+",\"image_stddev_bits\":"+bits(h.imageStddev)+",\"centroids_channel0_bits\":"+bits(h.centers[0])+",\"centroids_channel1_bits\":"+bits(h.centers[1])+",\"deviations_channel0_bits\":"+bits(h.deviations[0])+",\"deviations_channel1_bits\":"+bits(h.deviations[1])+",\"similarity_bits\":"+matrixBits(h.similarity)+",\"shape_sha256\":\""+hash(h.shape)+"\",\"shape_nonzero_count\":"+shapeNonzeroCount(h.shape)+",\"membership_sha256\":\""+hash(h.mems)+"\",\"denominator_sha256\":\""+hash(h.denominator)+"\",\"class_energy_sha256\":\""+hash(h.energy)+"\",\"objective_bits\":"+Float.floatToIntBits(h.objective)+",\"growth_event_count\":"+h.nextGrowthEvents.size()+",\"growth_events\":["+String.join(",",h.nextGrowthEvents)+"],\"post_growth_segmentation\":\""+labels(h.seg)+"\",\"post_growth_skeleton\":\""+gridLabels(h.skel,true)+"\",\"post_growth_ordering_bits\":"+gridBits(h.ordering)+",\"post_growth_volume_bits\":"+bits(h.nextVolumes)+"}");
        } catch(Exception e) { throw new RuntimeException(e); }
    }
    private static String matrixBits(float[][] values){StringBuilder b=new StringBuilder("[");boolean first=true;for(float[] row:values)for(float v:row){if(!first)b.append(',');first=false;b.append(Float.floatToIntBits(v));}return b.append(']').toString();}
    private static int shapeNonzeroCount(float[][][][] values){int count=0;for(float[][][] x:values)for(float[][] y:x)for(float[] z:y)for(float v:z)if(v!=0.0f)count++;return count;}
    private static String gridLabels(byte[][][] a,boolean empty){StringBuilder b=new StringBuilder();for(byte[][] x:a)for(byte[] y:x)for(byte v:y)b.append(empty&&v==EMPTY?'.':(v<10?(char)('0'+v):(char)('a'+v-10)));return b.toString();}
    private static String gridBits(float[][][] a){StringBuilder b=new StringBuilder("[");boolean first=true;for(float[][] x:a)for(float[] y:x)for(float v:y){if(!first)b.append(',');first=false;b.append(Float.floatToIntBits(v));}return b.append(']').toString();}
}
