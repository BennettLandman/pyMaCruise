import java.util.Locale;

/** Selected Java BLACKHOLE volume formula with asymmetric distance rules. */
public final class LesionToadsBlackHoleDistanceSweepReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    private static String run(boolean gmFlag) {
        final int nx=5,ny=4,nz=6,nc=3,bhClass=1;
        final float smoothing=0.2f,zero=1.0e-10f,infinity=1.0e10f,spread=1.25f;
        final int maxVent=3,maxGM=4,maxBrain=5;
        final float[] prior={0.25f,0.5f,0.75f};
        final String[] type={"OBJECT","BLACKHOLE","OBJECT"};
        final int[] relation={0,1,0};
        float[][][][] mems=new float[nx][ny][nz][nc],data=new float[nx][ny][nz][nc],shape=new float[nx][ny][nz][nc];
        float[][][] lesion=new float[nx][ny][nz],black=new float[nx][ny][nz],lesionData=new float[nx][ny][nz],blackData=new float[nx][ny][nz];
        int[][][] vent=new int[nx][ny][nz],gm=new int[nx][ny][nz],brain=new int[nx][ny][nz];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
            int sum=x+y+z;
            lesion[x][y][z]=0.25f+0.125f*x+0.0625f*y+0.03125f*z;
            black[x][y][z]=0.5f+0.0625f*x+0.03125f*y+0.015625f*z;
            lesionData[x][y][z]=0.75f+0.0625f*x+0.03125f*y+0.015625f*z;
            blackData[x][y][z]=0.5f+0.03125f*x+0.015625f*y+0.0078125f*z;
            vent[x][y][z]=(sum%3==0)?((x+y)%4):-1;
            gm[x][y][z]=((sum%3!=0)&&((x+2*y+z)%2==0))?(sum%5):-1;
            brain[x][y][z]=((x+2*y+z)%2==0)?(sum%6):-1;
            for(int k=0;k<nc;k++) {
                mems[x][y][z][k]=0.5f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                data[x][y][z][k]=0.75f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                shape[x][y][z][k]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z+0.0625f*k;
            }
        }
        shape[2][2][2][bhClass]=0.0f;
        int[][] delta={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out=new StringBuilder("{\"name\":\"").append(gmFlag?"gm_enabled":"gm_disabled").append("\",\"voxels\":[");
        boolean first=true;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
            if(!first)out.append(',');first=false;float den=0;
            for(int k=0;k<nc;k++) {
                if(shape[x][y][z][k]>0) {
                    float ngb=0,ngbL=0,ngbB=0,ngbS=0;
                    for(int[] f:delta) {
                        int a=x+f[0],b=y+f[1],c=z+f[2];
                        for(int m=0;m<nc;m++)if(m!=k){float v=mems[a][b][c][m];ngb+=v*v;}
                        float lv=lesion[a][b][c];ngbL+=lv*lv;float bv=black[a][b][c];ngbB+=bv*bv;
                    }
                    if(type[k].equals("BLACKHOLE"))for(int f:new int[]{1,0,3,2,5,4}) {
                        int a=x+delta[f][0],b=y+delta[f][1],c=z+delta[f][2];float v=mems[a][b][c][k];ngbS+=v*v;
                    }
                    float num=data[x][y][z][k],numL=(k==bhClass)?lesionData[x][y][z]:0,numB=(k==bhClass)?blackData[x][y][z]:0;
                    num+=prior[k]*smoothing*(ngb+ngbL+ngbB);
                    if(type[k].equals("BLACKHOLE")){numL+=prior[k]*smoothing*(ngb+ngbS+ngbB);numB+=prior[k]*smoothing*(ngb+ngbS+ngbL);}
                    if(num>zero)num=1.0f/num;else num=infinity;
                    if(relation[k]==1){float f=(1.0f<2.0f*shape[x][y][z][k])?1.0f:2.0f*shape[x][y][z][k];num=f*num;}
                    mems[x][y][z][k]=num;
                    if(type[k].equals("BLACKHOLE")) {
                        if(numL>zero)numL=1.0f/numL;else numL=infinity;
                        if(relation[k]==1){float f=(1.0f<2.0f*shape[x][y][z][k])?1.0f:2.0f*shape[x][y][z][k];numL=f*numL;}
                        if(numB>zero)numB=1.0f/numB;else numB=infinity;
                        if(relation[k]==1){float f=(1.0f<2.0f*shape[x][y][z][k])?1.0f:2.0f*shape[x][y][z][k];numB=f*numB;}
                        lesion[x][y][z]=numL;black[x][y][z]=numB;
                        float relationValue=(relation[k]==0)?1.0f:((1.0f<2.0f*shape[x][y][z][k])?1.0f:2.0f*shape[x][y][z][k]);
                        if(maxVent!=0&&vent[x][y][z]!=-1) {
                            float n=(vent[x][y][z]==0)?spread:vent[x][y][z];
                            lesion[x][y][z]=n*lesion[x][y][z]/maxVent;relationValue*=n/maxVent;
                        } else if(gmFlag&&maxGM!=0&&gm[x][y][z]!=-1) {
                            float n=(gm[x][y][z]==0)?spread:gm[x][y][z];
                            lesion[x][y][z]=n*lesion[x][y][z]/maxGM;black[x][y][z]=n*black[x][y][z]/maxGM;relationValue*=n/maxGM;
                        }
                        if(maxBrain!=0&&brain[x][y][z]!=-1) {
                            float n=brain[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/maxBrain;black[x][y][z]=n*black[x][y][z]/maxBrain;relationValue*=n/maxBrain;
                        }
                    }
                } else {mems[x][y][z][k]=0;if(type[k].equals("BLACKHOLE")){lesion[x][y][z]=0;black[x][y][z]=0;}}
                if(k==bhClass)den+=lesion[x][y][z]+black[x][y][z];den+=mems[x][y][z][k];
            }
            if(den>0)for(int k=0;k<nc;k++){mems[x][y][z][k]/=den;if(k==bhClass){lesion[x][y][z]/=den;black[x][y][z]/=den;}}
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for(int k=0;k<nc;k++){if(k>0)out.append(',');out.append(bits(mems[x][y][z][k]));}
            out.append("],\"lesion_bits\":").append(bits(lesion[x][y][z])).append(",\"black_hole_bits\":").append(bits(black[x][y][z])).append('}');
        }
        return out.append("]}").toString();
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"classes\":3,\"smoothing_bits\":1045220557,\"cases\":["+run(true)+","+run(false)+"]}\n");
    }
}
