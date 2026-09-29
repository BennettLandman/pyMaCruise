import java.util.Locale;

/** Selected VandyCRUISE OBJECT+OUTMASK smoothing and atlas-prior harness. */
public final class LesionToadsOutmaskPriorSweepReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx=5,ny=4,nz=6,nc=3,outClass=1;
        final float smoothing=0.2f,coefficient=0.65f,zero=1.0e-10f,infinity=1.0e10f;
        final float[] smoothPrior={0.25f,0.5f,0.75f},atlasPrior={0.5f,0.75f,1.25f};
        final float[][] similarity={{0.0f,0.2f,0.7f},{0.4f,0.0f,0.3f},{0.6f,0.1f,0.0f}};
        final String[] type={"OBJECT","OUTMASK","OBJECT"};
        final int[] relation={0,1,0};
        float[][][][] mems=new float[nx][ny][nz][nc],data=new float[nx][ny][nz][nc],shape=new float[nx][ny][nz][nc];
        float[][][] outmaskData=new float[nx][ny][nz],outlier=new float[nx][ny][nz];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
            outmaskData[x][y][z]=0.375f+0.03125f*x+0.015625f*y+0.0078125f*z;
            outlier[x][y][z]=0.125f*x+0.0625f*y+0.03125f*z;
            for(int k=0;k<nc;k++) {
                mems[x][y][z][k]=0.5f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                data[x][y][z][k]=0.75f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                shape[x][y][z][k]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z+0.0625f*k;
            }
        }
        outlier[1][1][1]=1.0f;outlier[3][2][4]=0.0f;shape[2][2][2][outClass]=0.0f;
        int[][] delta={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"classes\":3,\"smoothing_bits\":1045220557,\"prior_coefficient_bits\":1059481190,\"voxels\":[");
        boolean first=true;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
            if(!first)out.append(',');first=false;float den=0;
            for(int k=0;k<nc;k++) {
                if(shape[x][y][z][k]>0) {
                    float ngb=0;
                    for(int[] f:delta) {
                        int a=x+f[0],b=y+f[1],c=z+f[2];
                        for(int m=0;m<nc;m++)if(m!=k){float v=mems[a][b][c][m];ngb+=v*v;}
                    }
                    float num=data[x][y][z][k],numO=(k==outClass)?outmaskData[x][y][z]:0.0f;
                    num+=smoothPrior[k]*smoothing*ngb;
                    if(type[k].equals("OUTMASK"))numO+=smoothPrior[k]*smoothing*ngb;
                    float invprior=0.0f;
                    for(int m=0;m<nc;m++)if(m!=k) {
                        float product=atlasPrior[m]*similarity[k][m];product*=shape[x][y][z][m];product*=shape[x][y][z][m];invprior+=product;
                    }
                    num+=coefficient*invprior;
                    if(type[k].equals("OUTMASK"))numO+=coefficient*invprior;
                    if(num>zero)num=1.0f/num;else num=infinity;
                    if(type[k].equals("OUTMASK")) {
                        if(numO>zero)numO=1.0f/numO;else numO=infinity;
                        float base=(1.0f-outlier[x][y][z])*num+outlier[x][y][z]*numO;
                        float outValue=outlier[x][y][z]*numO;
                        if(relation[k]==1) {
                            float twice=2.0f*shape[x][y][z][k];float factor=(1.0f<twice)?1.0f:twice;
                            base=factor*base;outValue=(factor*outlier[x][y][z])*numO;
                        }
                        mems[x][y][z][k]=base;outlier[x][y][z]=outValue;
                    } else {
                        if(relation[k]==1){float twice=2.0f*shape[x][y][z][k];float factor=(1.0f<twice)?1.0f:twice;num=factor*num;}
                        mems[x][y][z][k]=num;
                    }
                } else mems[x][y][z][k]=0.0f;
                den+=mems[x][y][z][k];
            }
            if(den>0)for(int k=0;k<nc;k++){mems[x][y][z][k]/=den;if(k==outClass)outlier[x][y][z]/=den;}
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for(int k=0;k<nc;k++){if(k>0)out.append(',');out.append(bits(mems[x][y][z][k]));}
            out.append("],\"outlier_bits\":").append(bits(outlier[x][y][z])).append('}');
        }
        System.out.print(out.append("]}\n"));
    }
}
