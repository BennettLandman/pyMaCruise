import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Selected-source call-sequence reference for P18cyn. Relation, boundary,
 * closest-label, exact-membership, and first-thinning operations are executed
 * in that order using the companion method-excerpt harnesses and copied
 * selected-source membership/thinning formulas. This is a fixed-upstream
 * startup slice, not a standalone TOADS run.
 */
public final class ToadsInitialThinningStartupReferenceHarness {
    private static final int C = 11;
    private static final int[][] FACES = {{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};
    private static final String TOADS_SHA = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884";
    private static final String POWER_SHA = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678";
    private static final String LUT_SHA = "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f";
    private static final String HEAP_SHA = "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a";
    private static final float[][] CENTROID = {{.28f,.31f,.34f,.37f,.40f,.43f,.46f,.49f,.52f,.55f,.58f}, {.24f,.27f,.30f,.33f,.36f,.39f,.42f,.45f,.48f,.51f,.54f}};
    private static final float[][] VARIANCE = {{.81f,.93f,1.17f,.76f,.88f,1.02f,1.14f,.69f,.96f,1.21f,.73f}, {1.19f,.72f,.87f,1.08f,.79f,1.25f,.91f,1.13f,.68f,.98f,1.16f}};
    private static final float[] MODALITY = {.73f,.46f}, IMAGE_VARIANCE = {.61f,.34f};
    private static final float[] EXP_CSF = {.78f,.91f}, EXP_GM = {1.16f,1.29f}, EXP_WM = {1.55f,1.73f};
    private static final float[] SMOOTH_PRIOR = {.65f,.71f,.84f,.91f,1.02f,1.13f,.76f,.88f,1.05f,1.17f,.82f};
    private static final float[] ATLAS_PRIOR = {.10f,.08f,.09f,.11f,.12f,.09f,.10f,.11f,.07f,.06f,.07f};
    private static final float[][] SIM = {
        {0,.11f,.14f,.18f,.13f,.16f,.12f,.15f,.09f,.10f,.08f}, {.12f,0,.13f,.17f,.11f,.16f,.14f,.10f,.09f,.12f,.08f},
        {.15f,.12f,0,.13f,.18f,.12f,.16f,.11f,.10f,.07f,.09f}, {.14f,.16f,.12f,0,.13f,.18f,.11f,.14f,.08f,.10f,.12f},
        {.11f,.13f,.17f,.12f,0,.11f,.19f,.14f,.10f,.09f,.07f}, {.16f,.12f,.14f,.17f,.11f,0,.12f,.18f,.09f,.08f,.10f},
        {.13f,.15f,.11f,.14f,.18f,.12f,0,.13f,.08f,.09f,.11f}, {.12f,.10f,.14f,.16f,.13f,.18f,.11f,0,.07f,.09f,.08f},
        {.10f,.09f,.12f,.08f,.11f,.09f,.10f,.07f,0,.13f,.16f}, {.11f,.12f,.08f,.10f,.09f,.07f,.11f,.08f,.14f,0,.12f},
        {.08f,.10f,.09f,.12f,.07f,.11f,.10f,.09f,.17f,.13f,0}
    };
    private static final float SMOOTH=.043f, FACTOR=.64f, COEFFICIENT=.16f;
    private final int nx, ny, nz, relationLimit;
    private final byte[][][] segmentation;
    private final float[][][][] image, field, bias, initial, shapes;
    private float[][][][] memberships;
    private final byte[][][][] closest;
    private final short[][][] relationMap;
    private final boolean[][] relation;
    private final byte[][][][] boundary;
    private float[][][][] energy;
    private float[][][] denominator;
    private float objective;
    private final byte[][][] skeleton;
    private final boolean[][][][] available;
    private final float[] volumes, minimums;
    private float[] initialVolumes;
    private final int[] boundaryCounts;
    private final List<String> processed = new ArrayList<>();
    private final List<String> stageOrder = new ArrayList<>();
    private final BinaryHeap4D queue;
    private final float[][][][] attachedFieldByChannel;
    private final boolean channelFieldsAttached;

    private static byte[] gunzip(String path) throws Exception {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(path))) {
            byte[] buffer = new byte[8192]; int count;
            while ((count = in.read(buffer)) >= 0) out.write(buffer, 0, count);
        }
        return out.toByteArray();
    }
    private static ByteBuffer data(String path) throws Exception {
        return ByteBuffer.wrap(gunzip(path)).order(ByteOrder.LITTLE_ENDIAN);
    }
    private static int[] shape(ByteBuffer b, int expectedDims, int expectedDatatype) {
        int expectedBitpix=expectedDatatype==2?8:32;
        if (b.getInt(0) != 348 || (b.getShort(40) & 0xffff) != expectedDims
                || (b.getShort(70) & 0xffff) != expectedDatatype || b.getShort(72) != expectedBitpix)
            throw new IllegalArgumentException("unexpected NIfTI header");
        int[] shape = new int[expectedDims];
        for (int axis=0; axis<expectedDims; axis++) shape[axis] = b.getShort(42+2*axis);
        return shape;
    }
    private static float[][][][] readFloat4(String path, int channels) throws Exception {
        ByteBuffer b=data(path); int[] d=shape(b,4,16);
        if (d[3] != channels) throw new IllegalArgumentException("unexpected NIfTI channel count");
        int offset=(int)b.getFloat(108); float[][][][] out=new float[d[0]][d[1]][d[2]][channels];
        for(int c=0;c<channels;c++)for(int z=0;z<d[2];z++)for(int y=0;y<d[1];y++)for(int x=0;x<d[0];x++)
            out[x][y][z][c]=b.getFloat(offset+4*(((c*d[2]+z)*d[1]+y)*d[0]+x));
        return out;
    }
    private static byte[][][] readByte3(String path) throws Exception {
        ByteBuffer b=data(path); int[] d=shape(b,3,2); int offset=(int)b.getFloat(108);
        byte[][][] out=new byte[d[0]][d[1]][d[2]];
        for(int z=0;z<d[2];z++)for(int y=0;y<d[1];y++)for(int x=0;x<d[0];x++)
            out[x][y][z]=b.get(offset+(z*d[1]+y)*d[0]+x);
        return out;
    }
    private static String sha(byte[] bytes) throws Exception {
        byte[] hash=MessageDigest.getInstance("SHA-256").digest(bytes); StringBuilder s=new StringBuilder();
        for(byte b:hash)s.append(String.format("%02x",b&255)); return s.toString();
    }
    private static void putInt(MessageDigest d,int v){d.update(ByteBuffer.allocate(4).order(ByteOrder.BIG_ENDIAN).putInt(v).array());}
    private static String floatDigest(float[][][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(float[][][] p:a)for(float[][] r:p)for(float[] v:r)for(float q:v)putInt(d,Float.floatToIntBits(q));return hex(d.digest());}
    private static String floatDigest(float[][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(float[][] p:a)for(float[] r:p)for(float q:r)putInt(d,Float.floatToIntBits(q));return hex(d.digest());}
    private static String byteDigest(byte[][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(byte[][] p:a)for(byte[] r:p)d.update(r);return hex(d.digest());}
    private static String byte4Digest(byte[][][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(byte[][][] p:a)for(byte[][] r:p)for(byte[] v:r)d.update(v);return hex(d.digest());}
    private static String bool4Digest(boolean[][][][] a)throws Exception{StringBuilder s=new StringBuilder();for(boolean[][][] p:a)for(boolean[][] r:p)for(boolean[] v:r)for(boolean q:v)s.append(q?'1':'0');return sha(s.toString().getBytes(StandardCharsets.US_ASCII));}
    private static String short3Digest(short[][][] a)throws Exception{MessageDigest d=MessageDigest.getInstance("SHA-256");for(short[][] p:a)for(short[] r:p)for(short v:r)putInt(d,v);return hex(d.digest());}
    private static String bool2Digest(boolean[][] a)throws Exception{StringBuilder s=new StringBuilder();for(boolean[] r:a)for(boolean v:r)s.append(v?'1':'0');return sha(s.toString().getBytes(StandardCharsets.US_ASCII));}
    private static String hex(byte[] bytes){StringBuilder s=new StringBuilder();for(byte b:bytes)s.append(String.format("%02x",b&255));return s.toString();}
    private static int category(int k){return k<=2?0:k<=7?1:2;}
    private static float exponent(int category,int channel){return category==0?EXP_CSF[channel]:category==1?EXP_GM[channel]:EXP_WM[channel];}

    private ToadsInitialThinningStartupReferenceHarness(String dir) throws Exception {
        segmentation=readByte3(dir+"/segmentation_ras.nii.gz"); nx=segmentation.length;ny=segmentation[0].length;nz=segmentation[0][0].length;
        image=readFloat4(dir+"/images_ras.nii.gz",2);field=readFloat4(dir+"/fields_ras.nii.gz",2);
        bias=readFloat4(dir+"/psf_bias_ras.nii.gz",2);initial=readFloat4(dir+"/initial_memberships_ras.nii.gz",C);
        shapes=readFloat4(dir+"/transformed_shapes_ras.nii.gz",C);relationLimit=2;
        attachedFieldByChannel=new float[2][][][];float[][][][] channelInputs=new float[2][][][];boolean[] useField={false};
        for(int c=0;c<2;c++){
            channelInputs[c]=new float[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)channelInputs[c][x][y][z]=field[x][y][z][c];
            addInhomogeneityCorrection(attachedFieldByChannel,channelInputs[c],c,useField);
            stageOrder.add("addInhomogeneityCorrection[channel="+c+"]");
        }
        channelFieldsAttached=useField[0]&&attachedFieldByChannel[0]==channelInputs[0]&&attachedFieldByChannel[1]==channelInputs[1];
        skeleton=new byte[nx][ny][nz];available=new boolean[nx][ny][nz][C];volumes=new float[C];minimums=new float[C];boundaryCounts=new int[C];
        queue=new BinaryHeap4D(nx*ny*3,BinaryHeap4D.MAXTREE);

        ToadsRelationDiffusionReferenceHarness.Outcome relationOutcome =
            ToadsRelationDiffusionReferenceHarness.computeRelations(segmentation,C,initial,relationLimit);
        stageOrder.add("computeRelations");
        relationMap=relationOutcome.map;
        boolean[][] relationRows=ToadsRelationListReferenceHarness.build(segmentation,C).relations;
        relation=new boolean[C][C];
        for(int k=0;k<C;k++)System.arraycopy(relationRows[k],0,relation[k],0,C);
        ToadsBoundaryDistancesReferenceHarness.Outcome boundaryOutcome =
            ToadsBoundaryDistancesReferenceHarness.compute(segmentation,
                ToadsBoundaryDistancesReferenceHarness.initial(nx,ny,nz,C,false),C,relationLimit);
        stageOrder.add("computeBoundaryDistances");
        boundary=boundaryOutcome.distances;
        closest=ToadsClosestLabelsReferenceHarness.closestLabels(segmentation,relationLimit);
        stageOrder.add("computeClosestLabels");
        stageOrder.add("fixed:computeVariances");
        stageOrder.add("fixed:computeImageVariance+centroidSimilarity+transformedShapes");
        computeMemberships();
        initializeThinning();
        propagateThinning(3.0f);
    }

    private boolean direct(int x,int y,int z,int k){return closest[x][y][z][0]==k||closest[x][y][z][1]==k;}
    private boolean neighbor(int x,int y,int z,int k){return closest[x][y][z][2]==k||closest[x][y][z][3]==k;}
    private static void addInhomogeneityCorrection(float[][][][] fields,float[][][] field,int channel,boolean[] useField){fields[channel]=field;useField[0]=true;}
    private float relationFactor(int x,int y,int z,int k){
        float shape=shapes[x][y][z][k];
        if(direct(x,y,z,k))return 1.0f;
        if(neighbor(x,y,z,k))return Math.max(.25f,Math.min(2.0f*FACTOR*shape,1.0f));
        return Math.min(1.0f,2.0f*FACTOR*shape);
    }
    private float inversePrior(int x,int y,int z,int k){float value=0;for(int m=0;m<C;m++)if(m!=k)value+=ATLAS_PRIOR[m]*SIM[k][m]*shapes[x][y][z][m]*shapes[x][y][z][m];return value;}
    private void computeMemberships(){
        stageOrder.add("computeExactMemberships");
        PowerTable[][] powers=new PowerTable[3][2];
        for(int t=0;t<3;t++)for(int c=0;c<2;c++)powers[t][c]=new PowerTable(0.0f,4.0f,.000001f,exponent(t,c));
        energy=new float[nx][ny][nz][C]; memberships=new float[nx][ny][nz][C]; denominator=new float[nx][ny][nz];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<C;k++)memberships[x][y][z][k]=initial[x][y][z][k];
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++)for(int k=0;k<C;k++){
            if(!direct(x,y,z,k)&&!neighbor(x,y,z,k)&&shapes[x][y][z][k]==0.0f)continue;
            float value=0.0f;
            if(k==0){for(int c=0;c<2;c++){float weight=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];value+=weight*(image[x][y][z][c]==0.0f?0.0f:1.0f);}}
            else for(int c=0;c<2;c++){
                float corrected=attachedFieldByChannel[c][x][y][z]*image[x][y][z][c]-bias[x][y][z][c];
                float d=corrected-CENTROID[c][k],distance=d*d,weight=MODALITY[c]*VARIANCE[c][k]*IMAGE_VARIANCE[c];
                value=(float)(value+(double)weight*powers[category(k)][c].lookup((double)distance,(double)exponent(category(k),c)));
            }
            energy[x][y][z][k]=value;
        }
        int[][] faces={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            float den=0.0f;float[] objectiveAt=new float[C];
            for(int k=0;k<C;k++){
                boolean d=direct(x,y,z,k),n=neighbor(x,y,z,k);float sh=shapes[x][y][z][k];
                if(!d&&!n&&sh==0.0f){memberships[x][y][z][k]=0.0f;continue;}
                float num=energy[x][y][z][k],e=num;
                float neighbors=0.0f;for(int[] face:faces)for(int m=0;m<C;m++)if(m!=k){float q=memberships[x+face[0]][y+face[1]][z+face[2]][m];neighbors+=q*q;}
                num+=SMOOTH_PRIOR[k]*SMOOTH*neighbors;e+=.5f*SMOOTH_PRIOR[k]*SMOOTH*neighbors;
                float inverse=inversePrior(x,y,z,k);num+=COEFFICIENT*inverse;e+=.5f*COEFFICIENT*inverse;
                float factor=relationFactor(x,y,z,k),reciprocal=num>1.0e-10f?1.0f/num:1.0e10f;
                memberships[x][y][z][k]=factor*reciprocal;den+=memberships[x][y][z][k];objectiveAt[k]=e;
            }
            denominator[x][y][z]=den;
            if(den>0.0f)for(int k=0;k<C;k++){
                memberships[x][y][z][k]/=den;objective+=objectiveAt[k]*memberships[x][y][z][k]*memberships[x][y][z][k];
            }
        }
    }
    private int bestMembership(int x,int y,int z){int best=0;float max=0;for(int k=0;k<C;k++)if(memberships[x][y][z][k]>max){best=k;max=memberships[x][y][z][k];}return best;}
    private boolean topology(int x,int y,int z,int target){short[][][] copy=new short[nx][ny][nz];for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int l=0;l<nz;l++)copy[i][j][l]=segmentation[i][j][l];return ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(copy,(short)x,(short)y,(short)z,(short)target);}
    private void initializeThinning(){
        stageOrder.add("initThinningLabels");
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)skeleton[x][y][z]=-1;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            int k=segmentation[x][y][z];skeleton[x][y][z]=(byte)k;volumes[k]+=memberships[x][y][z][k];
            for(int c=0;c<C;c++)available[x][y][z][c]=c==k;
        }
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++){
            int k=segmentation[x][y][z];boolean isBoundary=false;int neighborCount=0,ni=0,nj=0,nl=0;
            for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++){
                int d2=i*i+j*j+l*l;
                if(d2==1&&segmentation[x+i][y+j][z+l]!=k)isBoundary=true;
                else if(d2>0){neighborCount++;ni=x+i;nj=y+j;nl=z+l;}
            }
            if(!isBoundary)continue;available[x][y][z][k]=false;boolean endpoint=false;
            if(neighborCount==1){neighborCount=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(i*i+j*j+l*l>0&&segmentation[ni+i][nj+j][nl+l]==k)neighborCount++;endpoint=neighborCount<=2;}
            float order;
            if(endpoint)order=-0.0f;else{int best=bestMembership(x,y,z);float ratio=k==best?1.0f:(memberships[x][y][z][best]+1.0e-10f)/(memberships[x][y][z][k]+1.0e-10f);order=-1.0f/ratio;}
            queue.addValue(order,x,y,z,k);boundaryCounts[k]++;
        }
        initialVolumes=volumes.clone();
        for(int k=0;k<C;k++)minimums[k]=Math.max(100.0f,.5f*volumes[k]);
    }
    private void processNext(float value){
        int x=queue.getFirstX(),y=queue.getFirstY(),z=queue.getFirstZ(),k=queue.getFirstK();
        processed.add(Float.floatToIntBits(value)+","+x+","+y+","+z+","+k+"\n");queue.removeFirst();
        if(volumes[k]>minimums[k]&&skeleton[x][y][z]==k){
            int id=bestMembership(x,y,z),next=-1;boolean direct=false;float bestPrior=-.1f;
            if(id!=k)for(int[] face:FACES){int xi=x+face[0],yj=y+face[1],zl=z+face[2],neighborLabel=segmentation[xi][yj][zl];
                if(neighborLabel==id)direct=true;
                else if(neighborLabel!=k&&relation[neighborLabel][id]&&shapes[xi][yj][zl][id]>bestPrior){next=neighborLabel;bestPrior=shapes[xi][yj][zl][id];}
            }
            boolean change=id!=k&&(direct||next!=-1);if(change){if(!direct)id=next;change=topology(x,y,z,id);}
            available[x][y][z][k]=false;volumes[k]-=memberships[x][y][z][k];skeleton[x][y][z]=-1;
            if(change)segmentation[x][y][z]=(byte)id;
            for(int[] face:FACES){int xi=x+face[0],yj=y+face[1],zl=z+face[2];if(!available[xi][yj][zl][k])continue;
                int count=0;for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++)for(int l=-1;l<=1;l++)if(face[0]*face[0]+face[1]*face[1]+face[2]*face[2]>0&&segmentation[xi+i][yj+j][zl+l]==k)count++;
                float order;if(count==1)order=value;else{int best=bestMembership(xi,yj,zl);float ratio=k==best?1.0f:(memberships[xi][yj][zl][best]+1.0e-10f)/(memberships[xi][yj][zl][k]+1.0e-10f);order=value-1.0f/ratio;}
                if(volumes[k]>minimums[k]&&order>-3.0f){available[xi][yj][zl][k]=false;queue.addValue(order,xi,yj,zl,k);}
            }
        }
    }
    private void propagateThinning(float limit){stageOrder.add("propagateThinning");float value=0;while(queue.isNotEmpty()&&value>-limit){value=queue.getFirst();processNext(value);}}
    private String processedDigest()throws Exception{return sha(String.join("",processed).getBytes(StandardCharsets.US_ASCII));}
    private String queueDigest()throws Exception{StringBuilder s=new StringBuilder();while(queue.isNotEmpty()){s.append(Float.floatToIntBits(queue.getFirst())).append(',').append(queue.getFirstX()).append(',').append(queue.getFirstY()).append(',').append(queue.getFirstZ()).append(',').append(queue.getFirstK()).append('\n');queue.removeFirst();}return sha(s.toString().getBytes(StandardCharsets.US_ASCII));}
    private static String bits(float[] values){StringBuilder s=new StringBuilder("[");for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(Float.floatToIntBits(values[i]));}return s.append(']').toString();}
    private String ints(int[] values){return java.util.Arrays.toString(values);}
    private String membershipOrder(){StringBuilder order=new StringBuilder("[");for(int i=0;i<stageOrder.size();i++){if(i>0)order.append(',');order.append('"').append(stageOrder.get(i)).append('"');}return order.append(']').toString();}
    private String json()throws Exception{
        int[] classes=new int[C];for(int k=0;k<C;k++)classes[k]=k;
        StringBuilder out=new StringBuilder("{\"source\":\"selected VandyCRUISE startup relation refresh -> exact memberships -> initThinningLabels + propagateThinning\"");
        out.append(",\"toads_source_sha256\":\"").append(TOADS_SHA).append("\",\"powertable_source_sha256\":\"").append(POWER_SHA).append("\",\"lut_sha256\":\"").append(LUT_SHA).append("\",\"heap_sha256\":\"").append(HEAP_SHA).append('\"');
        out.append(",\"shape\":[").append(nx).append(',').append(ny).append(',').append(nz).append(",11],\"relation_limit\":").append(relationLimit).append(",\"volume_ratio_bits\":").append(Float.floatToIntBits(.5f)).append(",\"thinning_limit_bits\":").append(Float.floatToIntBits(3.0f));
        out.append(",\"stage_order\":").append(membershipOrder()).append(",\"field_references_attached\":").append(channelFieldsAttached);
        out.append(",\"relation_map_sha256\":\"").append(short3Digest(relationMap)).append('\"');
        out.append(",\"relation_adjacency_sha256\":\"").append(bool2Digest(relation)).append('\"');
        out.append(",\"boundary_distances_sha256\":\"").append(byte4Digest(boundary)).append('\"');
        out.append(",\"closest_labels_sha256\":\"").append(byte4Digest(closest)).append('\"');
        out.append(",\"starting_membership_sha256\":\"").append(floatDigest(initial)).append('\"');
        out.append(",\"transformed_shape_samples_sha256\":\"").append(floatDigest(shapes)).append('\"');
        out.append(",\"fixed_upstream\":\"class centroids, class variance weights, image variance weights, centroid similarities, modality weights, category exponents, atlas shape priors, and the already-transformed shape sample map\"");
        out.append(",\"registration_transform\":\"not computed or consumed; synthetic transformed-shape samples are directly supplied on the target RAS grid\"");
        out.append(",\"class_energy_sha256\":\"").append(floatDigest(energy)).append('\"');
        out.append(",\"membership_sha256\":\"").append(floatDigest(memberships)).append('\"');
        out.append(",\"denominator_sha256\":\"").append(floatDigest(denominator)).append('\"');
        out.append(",\"objective_bits\":").append(Float.floatToIntBits(objective));
        out.append(",\"boundary_counts\":").append(ints(boundaryCounts));
        out.append(",\"initial_volume_bits\":").append(bits(volumesInitial()));
        out.append(",\"minimum_volume_bits\":").append(bits(minimums));
        out.append(",\"processed_count\":").append(processed.size()).append(",\"processed_sha256\":\"").append(processedDigest()).append('\"');
        out.append(",\"final_volume_bits\":").append(bits(volumes));
        out.append(",\"segmentation_sha256\":\"").append(byteDigest(segmentation)).append('\"');
        out.append(",\"skeleton_sha256\":\"").append(byteDigest(skeleton)).append('\"');
        out.append(",\"available_sha256\":\"").append(bool4Digest(available)).append('\"');
        int queueCount=0;while(queue.isNotEmpty()){queueCount++;queue.removeFirst();}
        out.append(",\"queue_after_count\":").append(queueCount).append(",\"queue_after_sha256\":\"").append(sha(new byte[0])).append("\"}");return out.toString();
    }
    private float[] volumesInitial(){return initialVolumes;}
    public static void main(String[] args)throws Exception{
        if(args.length!=1)throw new IllegalArgumentException("fixture directory required");
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        ToadsInitialThinningStartupReferenceHarness h=new ToadsInitialThinningStartupReferenceHarness(args[0]);
        System.out.println(h.json());
    }
}
