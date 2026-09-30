import java.util.Locale;

/** Selected VandyCRUISE OBJECT+LESION+OUTMASK in-place formula harness. */
public final class LesionToadsLesionOutmaskSweepReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }
    private static float factor(float shape) { float doubled=2.0f*shape;return (1.0f<doubled)?1.0f:doubled; }
    private static boolean isLesion(String type) { return type.equals("LESION") || type.equals("OPTIMIZED_LESION"); }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx=5,ny=4,nz=6,nc=4;
        final boolean optimized=args.length>0&&args[0].startsWith("optimized-");
        final boolean outmaskFirst=args.length>0&&(args[0].equals("optimized-outmask-first")||args[0].equals("outmask-first"));
        final String lesionType=optimized?"OPTIMIZED_LESION":"LESION";
        final float smoothing=0.2f,coefficient=0.65f,zero=1.0e-10f,infinity=1.0e10f;
        final float[] smoothPrior={0.25f,0.5f,0.75f,0.375f};
        final float[] atlasPrior={0.5f,0.75f,1.25f,0.875f};
        final float[] lesionSimilarity={0.25f,0.0f,0.5f,0.375f};
        final float[][] similarity={{0.0f,0.2f,0.7f,0.3f},{0.4f,0.0f,0.3f,0.6f},{0.6f,0.1f,0.0f,0.2f},{0.2f,0.5f,0.4f,0.0f}};
        final String[] type=outmaskFirst?new String[]{"OBJECT","OUTMASK",lesionType,"OBJECT"}:new String[]{"OBJECT",lesionType,"OUTMASK","OBJECT"};
        int lesionClass=-1,outClass=-1;
        for(int k=0;k<nc;k++){if(isLesion(type[k]))lesionClass=k;if(type[k].equals("OUTMASK"))outClass=k;}
        final int[] relation={0,1,1,0};
        float[][][][] mems=new float[nx][ny][nz][nc],data=new float[nx][ny][nz][nc],shape=new float[nx][ny][nz][nc];
        float[][][] lesionData=new float[nx][ny][nz],outmaskData=new float[nx][ny][nz],lesion=new float[nx][ny][nz],outlier=new float[nx][ny][nz];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
            lesionData[x][y][z]=0.625f+0.03125f*x+0.015625f*y+0.0078125f*z;
            outmaskData[x][y][z]=0.375f+0.03125f*x+0.015625f*y+0.0078125f*z;
            lesion[x][y][z]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z;
            outlier[x][y][z]=0.125f*x+0.0625f*y+0.03125f*z;
            for(int k=0;k<nc;k++) {
                mems[x][y][z][k]=0.5f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                data[x][y][z][k]=0.75f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                shape[x][y][z][k]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z+0.0625f*k;
            }
        }
        lesion[1][1][1]=0.0f;outlier[1][1][1]=1.0f;outlier[3][2][4]=0.0f;shape[2][2][2][lesionClass]=0.0f;
        int[][] delta={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"classes\":4,\"class_order\":\"").append(outmaskFirst?"outmask-before-lesion":"lesion-before-outmask").append("\",\"optimized_lesion\":").append(optimized).append(",\"smoothing_bits\":1045220557,\"prior_coefficient_bits\":1059481190,\"voxels\":[");
        boolean first=true;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
            if(!first)out.append(',');first=false;float den=0.0f;
            for(int k=0;k<nc;k++) {
                if(shape[x][y][z][k]>0.0f) {
                    float other=0.0f,lesionSum=0.0f,currentSum=0.0f;
                    for(int[] f:delta) {
                        int a=x+f[0],b=y+f[1],c=z+f[2];
                        for(int m=0;m<nc;m++)if(m!=k){float v=mems[a][b][c][m];other+=v*v;}
                        float lv=lesion[a][b][c];lesionSum+=lv*lv;
                    }
                    if(isLesion(type[k])) {
                        int[][] lesionOrder={{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
                        for(int[] f:lesionOrder){float v=mems[x+f[0]][y+f[1]][z+f[2]][k];currentSum+=v*v;}
                    }
                    float num=data[x][y][z][k],numL=isLesion(type[k])?lesionData[x][y][z]:0.0f;
                    float numO=type[k].equals("OUTMASK")?outmaskData[x][y][z]:0.0f;
                    float all=other+lesionSum;
                    num+=smoothPrior[k]*smoothing*all;
                    if(isLesion(type[k]))numL+=smoothPrior[k]*smoothing*(other+currentSum);
                    if(type[k].equals("OUTMASK"))numO+=smoothPrior[k]*smoothing*(other+lesionSum);
                    float invprior=0.0f,invpriorL=0.0f;
                    for(int m=0;m<nc;m++) {
                        if(m!=k){float product=atlasPrior[m]*similarity[k][m];product*=shape[x][y][z][m];product*=shape[x][y][z][m];invprior+=product;}
                        if(isLesion(type[k])){float product=atlasPrior[m]*lesionSimilarity[m];product*=shape[x][y][z][m];product*=shape[x][y][z][m];invpriorL+=product;}
                    }
                    float shapeLesion=shape[x][y][z][lesionClass];
                    float product=atlasPrior[lesionClass]*lesionSimilarity[k];product*=shapeLesion;product*=shapeLesion;invprior+=product;
                    num+=coefficient*invprior;
                    if(isLesion(type[k]))numL+=coefficient*invpriorL;
                    if(type[k].equals("OUTMASK"))numO+=coefficient*invprior;
                    if(num>zero)num=1.0f/num;else num=infinity;
                    if(isLesion(type[k])){if(numL>zero)numL=1.0f/numL;else numL=infinity;}
                    if(type[k].equals("OUTMASK")){if(numO>zero)numO=1.0f/numO;else numO=infinity;}
                    if(isLesion(type[k])) {
                        mems[x][y][z][k]=(relation[k]==1)?factor(shape[x][y][z][k])*num:num;
                        lesion[x][y][z]=(relation[k]==1)?factor(shape[x][y][z][k])*numL:numL;
                    } else if(type[k].equals("OUTMASK")) {
                        float base=(1.0f-outlier[x][y][z])*num+outlier[x][y][z]*numO;
                        if(relation[k]==1){float fac=factor(shape[x][y][z][k]);base=fac*base;outlier[x][y][z]=fac*outlier[x][y][z]*numO;}
                        else outlier[x][y][z]=outlier[x][y][z]*numO;
                        mems[x][y][z][k]=base;
                    } else mems[x][y][z][k]=num;
                } else {
                    mems[x][y][z][k]=0.0f;
                    if(k==lesionClass)lesion[x][y][z]=0.0f;
                }
                if(k==lesionClass)den+=lesion[x][y][z];
                den+=mems[x][y][z][k];
            }
            if(den>0.0f)for(int k=0;k<nc;k++) {
                mems[x][y][z][k]/=den;
                if(k==lesionClass)lesion[x][y][z]/=den;
                else if(k==outClass)outlier[x][y][z]/=den;
            }
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for(int k=0;k<nc;k++){if(k>0)out.append(',');out.append(bits(mems[x][y][z][k]));}
            out.append("],\"lesion_bits\":").append(bits(lesion[x][y][z])).append(",\"outlier_bits\":").append(bits(outlier[x][y][z])).append('}');
        }
        System.out.print(out.append("]}\n"));
    }
}
