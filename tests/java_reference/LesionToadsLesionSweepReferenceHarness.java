import java.util.Locale;

/** Selected VandyCRUISE OBJECT+LESION in-place smoothing sweep harness. */
public final class LesionToadsLesionSweepReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx=5, ny=4, nz=6, nc=3, lesionClass=1;
        final float smoothing=0.2f, zero=1.0e-10f, infinity=1.0e10f;
        final float[] prior={0.25f,0.5f,0.75f};
        final String[] type={"OBJECT","LESION","OBJECT"};
        final int[] relation={0,1,0};
        float[][][][] mems=new float[nx][ny][nz][nc];
        float[][][][] data=new float[nx][ny][nz][nc];
        float[][][][] shape=new float[nx][ny][nz][nc];
        float[][][] lesion=new float[nx][ny][nz];
        float[][][] lesionData=new float[nx][ny][nz];
        for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++) {
            lesion[x][y][z]=0.25f+0.125f*x+0.0625f*y+0.03125f*z;
            lesionData[x][y][z]=0.75f+0.0625f*x+0.03125f*y+0.015625f*z;
            for(int k=0;k<nc;k++) {
                mems[x][y][z][k]=0.5f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                data[x][y][z][k]=0.75f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                shape[x][y][z][k]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z+0.0625f*k;
            }
        }
        shape[2][2][2][lesionClass]=0.0f;
        float[][][] blackHole=new float[nx][ny][nz];
        int[][] delta={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out=new StringBuilder();
        out.append("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"classes\":3,\"smoothing_bits\":").append(bits(smoothing)).append(",\"voxels\":[");
        boolean first=true;
        for(int x=1;x<nx-1;x++) for(int y=1;y<ny-1;y++) for(int z=1;z<nz-1;z++) {
            if(!first) out.append(','); first=false;
            float den=0.0f;
            for(int k=0;k<nc;k++) {
                if(shape[x][y][z][k]>0.0f) {
                    float num=data[x][y][z][k];
                    float numL=(k==lesionClass)?lesionData[x][y][z]:0.0f;
                    float ngb=0.0f, ngbL=0.0f, ngbB=0.0f, ngbS=0.0f;
                    for(int[] face:delta) {
                        int a=x+face[0], b=y+face[1], c=z+face[2];
                        for(int m=0;m<nc;m++) if(m!=k) {
                            float value=mems[a][b][c][m]; ngb += value*value;
                        }
                        float l=lesion[a][b][c]; ngbL += l*l;
                    }
                    if(type[k].equals("LESION")) {
                        for(int f:new int[]{1,0,3,2,5,4}) {
                            int a=x+delta[f][0], b=y+delta[f][1], c=z+delta[f][2];
                            float value=mems[a][b][c][k]; ngbS += value*value;
                        }
                    }
                    num += prior[k]*smoothing*(ngb+ngbL+ngbB);
                    if(type[k].equals("LESION")) numL += prior[k]*smoothing*(ngb+ngbS+ngbB);
                    if(num>zero) num=1.0f/num; else num=infinity;
                    if(relation[k]==1) {
                        float twice=2.0f*shape[x][y][z][k];
                        float factor=(1.0f<twice)?1.0f:twice;
                        num=factor*num;
                    }
                    mems[x][y][z][k]=num;
                    if(type[k].equals("LESION")) {
                        if(numL>zero) numL=1.0f/numL; else numL=infinity;
                        if(relation[k]==1) {
                            float twice=2.0f*shape[x][y][z][k];
                            float factor=(1.0f<twice)?1.0f:twice;
                            numL=factor*numL;
                        }
                        lesion[x][y][z]=numL;
                    }
                } else {
                    mems[x][y][z][k]=0.0f;
                    if(type[k].equals("LESION")) lesion[x][y][z]=0.0f;
                }
                if(k==lesionClass) den += lesion[x][y][z]+blackHole[x][y][z];
                den += mems[x][y][z][k];
            }
            if(den>0.0f) for(int k=0;k<nc;k++) {
                mems[x][y][z][k]=mems[x][y][z][k]/den;
                if(k==lesionClass) lesion[x][y][z]=lesion[x][y][z]/den;
            }
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for(int k=0;k<nc;k++) { if(k>0) out.append(','); out.append(bits(mems[x][y][z][k])); }
            out.append("],\"lesion_bits\":").append(bits(lesion[x][y][z])).append('}');
        }
        out.append("]}\n");
        System.out.print(out);
    }
}
