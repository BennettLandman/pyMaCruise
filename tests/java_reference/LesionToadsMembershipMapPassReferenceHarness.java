import java.io.File;
import java.util.Locale;
import java.util.Scanner;

/** Selected Java computeMemberships map sweep, driven by NIfTI-exported floats. */
public final class LesionToadsMembershipMapPassReferenceHarness {
    private static final int NX=7, NY=8, NZ=9, NC=5, CHANNELS=24;
    private static final String[] TYPES={"OBJECT","LESION","BLACKHOLE","OUTMASK","OBJECT"};
    private static final int[] RELATION={0,1,1,1,0}; // 0 DIRECT, 1 PRIOR
    private static final int[][] FACES={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    private static int bits(float v){return Float.floatToIntBits(v);}
    private static float factor(float shape){float twice=2.0f*shape;return 1.0f<twice?1.0f:twice;}
    private static boolean lesionType(int k){return TYPES[k].equals("LESION")||TYPES[k].equals("OPTIMIZED_LESION");}
    public static void main(String[] args)throws Exception{
        Locale.setDefault(Locale.ROOT);
        if(args.length!=1)throw new IllegalArgumentException("pass flat x/y/z/channel float32 text export");
        try(Scanner scan=new Scanner(new File(args[0]))){
            scan.useLocale(Locale.ROOT);
            float[][][][] all=new float[NX][NY][NZ][CHANNELS];
            for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<CHANNELS;c++)all[x][y][z][c]=scan.nextFloat();
            run(all);
        }
    }
    private static void run(float[][][][] all){
        float[][][][] mem=new float[NX][NY][NZ][NC],data=new float[NX][NY][NZ][NC],shape=new float[NX][NY][NZ][NC];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int k=0;k<NC;k++){
            mem[x][y][z][k]=all[x][y][z][k];data[x][y][z][k]=all[x][y][z][8+k];shape[x][y][z][k]=all[x][y][z][16+k];
        }
        float[][][] lesion=new float[NX][NY][NZ],bh=new float[NX][NY][NZ],outlier=new float[NX][NY][NZ];
        float[][][] lesionData=new float[NX][NY][NZ],bhData=new float[NX][NY][NZ],outmaskData=new float[NX][NY][NZ];
        int[][][] vent=new int[NX][NY][NZ],gm=new int[NX][NY][NZ],bstem=new int[NX][NY][NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            lesion[x][y][z]=all[x][y][z][5];bh[x][y][z]=all[x][y][z][6];outlier[x][y][z]=all[x][y][z][7];
            lesionData[x][y][z]=all[x][y][z][13];bhData[x][y][z]=all[x][y][z][14];outmaskData[x][y][z]=all[x][y][z][15];
            vent[x][y][z]=(int)all[x][y][z][21];gm[x][y][z]=(int)all[x][y][z][22];bstem[x][y][z]=(int)all[x][y][z][23];
        }
        final float smoothing=0.2f, coefficient=0.65f, zero=1.0e-10f, infinity=1.0e10f;
        final float[] smoothPrior={0.25f,0.5f,0.75f,0.375f,0.625f};
        final float[] atlasPrior={0.5f,0.75f,1.25f,0.875f,0.625f};
        final float[] lesionSimilarity={0.25f,0.0f,0.5f,0.375f,0.625f};
        final float[] bhSimilarity={0.125f,0.5f,0.0f,0.25f,0.375f};
        final float[][] similarity={{0.0f,0.2f,0.7f,0.3f,0.1f},{0.4f,0.0f,0.3f,0.6f,0.2f},{0.6f,0.1f,0.0f,0.2f,0.5f},{0.2f,0.5f,0.4f,0.0f,0.3f},{0.3f,0.4f,0.1f,0.7f,0.0f}};
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[7,8,9],\"class_order\":[\"OBJECT\",\"LESION\",\"BLACKHOLE\",\"OUTMASK\",\"OBJECT\"],\"voxels\":[");
        boolean first=true;
        for(int x=1;x<NX-1;x++)for(int y=1;y<NY-1;y++)for(int z=1;z<NZ-1;z++){
            float den=0.0f;
            for(int k=0;k<NC;k++){
                if(shape[x][y][z][k]>0.0f){
                    float other=0.0f,lSum=0.0f,bSum=0.0f,current=0.0f;
                    for(int[] f:FACES){int a=x+f[0],b=y+f[1],c=z+f[2];for(int m=0;m<NC;m++)if(m!=k){float v=mem[a][b][c][m];other+=v*v;}float lv=lesion[a][b][c],bv=bh[a][b][c];lSum+=lv*lv;bSum+=bv*bv;}
                    if(lesionType(k)||TYPES[k].equals("BLACKHOLE"))for(int[] f:new int[][]{{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}}){float v=mem[x+f[0]][y+f[1]][z+f[2]][k];current+=v*v;}
                    float n=data[x][y][z][k],nL=(lesionType(k)||TYPES[k].equals("BLACKHOLE"))?lesionData[x][y][z]:0.0f;
                    float nB=TYPES[k].equals("BLACKHOLE")?bhData[x][y][z]:0.0f,nO=TYPES[k].equals("OUTMASK")?outmaskData[x][y][z]:0.0f;
                    float allSmooth=other+lSum;allSmooth+=bSum;n+=smoothPrior[k]*smoothing*allSmooth;
                    if(lesionType(k)||TYPES[k].equals("BLACKHOLE")){float ls=other+current;ls+=bSum;nL+=smoothPrior[k]*smoothing*ls;}
                    if(TYPES[k].equals("BLACKHOLE")){float bs=other+current;bs+=lSum;nB+=smoothPrior[k]*smoothing*bs;}
                    if(TYPES[k].equals("OUTMASK"))nO+=smoothPrior[k]*smoothing*(other+lSum);
                    float inv=0.0f,invL=0.0f,invB=0.0f;
                    for(int m=0;m<NC;m++){
                        if(m!=k){float q=atlasPrior[m]*similarity[k][m];q*=shape[x][y][z][m];q*=shape[x][y][z][m];inv+=q;}
                        if(lesionType(k)||TYPES[k].equals("BLACKHOLE")){float q=atlasPrior[m]*lesionSimilarity[m];q*=shape[x][y][z][m];q*=shape[x][y][z][m];invL+=q;}
                        if(TYPES[k].equals("BLACKHOLE")){float q=atlasPrior[m]*bhSimilarity[m];q*=shape[x][y][z][m];q*=shape[x][y][z][m];invB+=q;}
                    }
                    int lesionClass=2;float selected=shape[x][y][z][lesionClass];float q=atlasPrior[lesionClass]*lesionSimilarity[k];q*=selected;q*=selected;inv+=q;
                    if(TYPES[lesionClass].equals("BLACKHOLE")){q=atlasPrior[lesionClass]*bhSimilarity[k];q*=selected;q*=selected;inv+=q;}
                    n+=coefficient*inv;if(lesionType(k)||TYPES[k].equals("BLACKHOLE"))nL+=coefficient*invL;if(TYPES[k].equals("BLACKHOLE"))nB+=coefficient*invB;if(TYPES[k].equals("OUTMASK"))nO+=coefficient*inv;
                    n=n>zero?1.0f/n:infinity;if(lesionType(k)||TYPES[k].equals("BLACKHOLE"))nL=nL>zero?1.0f/nL:infinity;if(TYPES[k].equals("BLACKHOLE"))nB=nB>zero?1.0f/nB:infinity;if(TYPES[k].equals("OUTMASK"))nO=nO>zero?1.0f/nO:infinity;
                    if(lesionType(k)){
                        mem[x][y][z][k]=RELATION[k]==1?factor(shape[x][y][z][k])*n:n;
                        lesion[x][y][z]=RELATION[k]==1?factor(shape[x][y][z][k])*nL:nL;
                        if(vent[x][y][z]!=-1){float d=vent[x][y][z]==0?2.5f:vent[x][y][z];lesion[x][y][z]=d*lesion[x][y][z]/5.0f;}
                        else if(gm[x][y][z]!=-1){float d=gm[x][y][z]==0?2.5f:gm[x][y][z];lesion[x][y][z]=d*lesion[x][y][z]/8.0f;}
                        if(bstem[x][y][z]!=-1){float d=bstem[x][y][z];lesion[x][y][z]=d*lesion[x][y][z]/4.0f;}
                    }else if(TYPES[k].equals("BLACKHOLE")){
                        mem[x][y][z][k]=RELATION[k]==1?factor(shape[x][y][z][k])*n:n;
                        lesion[x][y][z]=RELATION[k]==1?factor(shape[x][y][z][k])*nL:nL;bh[x][y][z]=RELATION[k]==1?factor(shape[x][y][z][k])*nB:nB;
                        if(vent[x][y][z]!=-1){float d=vent[x][y][z]==0?2.5f:vent[x][y][z];lesion[x][y][z]=d*lesion[x][y][z]/5.0f;}
                        else if(gm[x][y][z]!=-1){float d=gm[x][y][z]==0?2.5f:gm[x][y][z];lesion[x][y][z]=d*lesion[x][y][z]/8.0f;bh[x][y][z]=d*bh[x][y][z]/8.0f;}
                        if(bstem[x][y][z]!=-1){float d=bstem[x][y][z];lesion[x][y][z]=d*lesion[x][y][z]/4.0f;bh[x][y][z]=d*bh[x][y][z]/4.0f;}
                    }else if(TYPES[k].equals("OUTMASK")){
                        float base=(1.0f-outlier[x][y][z])*n+outlier[x][y][z]*nO;float ov=outlier[x][y][z]*nO;
                        if(RELATION[k]==1){float f=factor(shape[x][y][z][k]);base=f*base;ov=(f*outlier[x][y][z])*nO;}
                        mem[x][y][z][k]=base;outlier[x][y][z]=ov;
                    }else mem[x][y][z][k]=n;
                }else{
                    mem[x][y][z][k]=0.0f;if(lesionType(k)||TYPES[k].equals("BLACKHOLE"))lesion[x][y][z]=0.0f;if(TYPES[k].equals("BLACKHOLE"))bh[x][y][z]=0.0f;
                }
                if(lesionType(k))den+=lesion[x][y][z];else if(TYPES[k].equals("BLACKHOLE")){float aux=lesion[x][y][z]+bh[x][y][z];den+=aux;}den+=mem[x][y][z][k];
            }
            if(den>0.0f)for(int k=0;k<NC;k++){mem[x][y][z][k]/=den;if(lesionType(k)||TYPES[k].equals("BLACKHOLE"))lesion[x][y][z]/=den;if(TYPES[k].equals("BLACKHOLE"))bh[x][y][z]/=den;if(TYPES[k].equals("OUTMASK"))outlier[x][y][z]/=den;}
            if(!first)out.append(',');first=false;out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for(int k=0;k<NC;k++){if(k>0)out.append(',');out.append(bits(mem[x][y][z][k]));}
            out.append("],\"lesion_bits\":").append(bits(lesion[x][y][z])).append(",\"black_hole_bits\":").append(bits(bh[x][y][z])).append(",\"outlier_bits\":").append(bits(outlier[x][y][z])).append('}');
        }
        System.out.println(out.append("]}").toString());
    }
}
