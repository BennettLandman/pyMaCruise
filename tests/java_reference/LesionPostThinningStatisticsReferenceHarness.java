import edu.jhmi.rad.medic.utilities.PowerTable;
import java.io.DataInputStream;
import java.io.FileInputStream;

/**
 * Fixed-state selected-source excerpt for statistics after P39 initial thinning.
 * Input floats/ints are big-endian and use x/y/z/channel/class array order.
 */
public final class LesionPostThinningStatisticsReferenceHarness {
    static final int N=17, NC=2, K=13, LESION_CLASS=10;
    static final String SOURCE_SHA="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    static final String WRAPPER_SHA="556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    static final String POWER_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    static final float EPS=1.0e-30f;
    static final float[][][][] image=new float[N][N][N][NC];
    static final float[][][][] mem=new float[N][N][N][K];
    static final int[][][] seg=new int[N][N][N];
    static final float[][][] lesion=new float[N][N][N], black=new float[N][N][N];
    static final float[][] initial=new float[NC][K];
    static final float[] initialLesion=new float[NC], initialBlack=new float[NC];
    static final float[] lesionPower={0.4f,1.2f};

    static int bits(float value) { return Float.floatToRawIntBits(value); }
    static String vector(float[] values) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(bits(values[i]));}
        return s.append(']').toString();
    }
    static String matrix(float[][] values) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<values.length;i++){if(i>0)s.append(',');s.append(vector(values[i]));}
        return s.append(']').toString();
    }

    static void readInput(String path) throws Exception {
        try(DataInputStream input=new DataInputStream(new FileInputStream(path))) {
            for(int c=0;c<NC;c++) for(int k=0;k<K;k++) initial[c][k]=input.readFloat();
            for(int c=0;c<NC;c++) initialLesion[c]=input.readFloat();
            for(int c=0;c<NC;c++) initialBlack[c]=input.readFloat();
            for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
                for(int c=0;c<NC;c++) image[x][y][z][c]=input.readFloat();
                for(int k=0;k<K;k++) mem[x][y][z][k]=input.readFloat();
                seg[x][y][z]=input.readInt();
                lesion[x][y][z]=input.readFloat();
                black[x][y][z]=input.readFloat();
            }
            if(input.read()!=-1) throw new IllegalArgumentException("trailing bytes in harness input");
        }
    }

    static Object[] refresh() {
        float[][] centroid=new float[NC][K];
        float[] lesionCentroid=initialLesion.clone(), blackCentroid=initialBlack.clone();
        PowerTable lesionTable=new PowerTable(0.0f,4.0f,0.000001f,lesionPower[1]);
        for(int c=0;c<NC;c++) {
            float[] num=new float[K], den=new float[K];
            float numL=0.0f,denL=0.0f,numB=0.0f,denB=0.0f;
            for(int k=0;k<K;k++) for(int x=0;x<N;x++) for(int y=0;y<N;y++) for(int z=0;z<N;z++) if(seg[x][y][z]==k) {
                float m=mem[x][y][z][k], sample=image[x][y][z][c];
                if(k!=LESION_CLASS) { num[k]+=m*m*sample; den[k]+=m*m; }
                else {
                    float lm=lesion[x][y][z], bm=black[x][y][z], winner=Math.max(m,Math.max(lm,bm));
                    float factor=1.0f;
                    // computeCentroids routes by max(mem, lesion, black-hole),
                    // with ties resolved in this order: class, lesion, black-hole.
                    if(winner==m) { num[k]+=m*m*sample; den[k]+=m*m; }
                    else if(winner==lm && c==1) {
                        float d=sample-initialLesion[c]; d*=d;
                        if(d>1e-10f) factor=(float)(lesionTable.lookup((double)d,(double)lesionPower[c])/d);
                    }
                    if(winner==lm) { numL+=factor*lm*lm*sample; denL+=factor*lm*lm; }
                    else if(winner==bm) { numB+=bm*bm*sample; denB+=bm*bm; }
                }
            }
            for(int k=0;k<K;k++) centroid[c][k]=den[k]>0.0f?num[k]/den[k]:0.0f;
            lesionCentroid[c]=denL>0.0f?numL/denL:0.0f;
            blackCentroid[c]=denB>0.0f?numB/denB:0.0f;
        }

        float[][] std=new float[NC][K]; float[] lesionStd=new float[NC], blackStd=new float[NC];
        for(int c=0;c<NC;c++) {
            float lnum=0.0f,lden=0.0f,bnum=0.0f,bden=0.0f;
            for(int k=0;k<K;k++) {
                float num=0.0f,den=0.0f;
                for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) if(seg[x][y][z]==k) {
                    float m=mem[x][y][z][k],lm=lesion[x][y][z],bm=black[x][y][z],sample=image[x][y][z][c];
                    if(k==LESION_CLASS && m<lm) {
                        if(lm>=bm) { float d=sample-lesionCentroid[c],w=lm*lm+EPS; lnum+=w*d*d; lden+=w; }
                        else { float d=sample-blackCentroid[c],w=lm*lm+EPS; bnum+=w*d*d; bden+=w; }
                    } else {
                        float d=sample-centroid[c][k],w=m*m+EPS; num+=w*d*d; den+=w;
                    }
                }
                std[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
            }
            lesionStd[c]=lden>0.0f?(float)Math.sqrt(lnum/lden):0.0f;
            blackStd[c]=bden>0.0f?(float)Math.sqrt(bnum/bden):0.0f;
        }

        float[] channelWeight=new float[NC]; float sum=0.0f;
        for(int c=0;c<NC;c++) {
            float num=0.0f,den=0.0f;
            for(int k=0;k<K;k++) for(int x=0;x<N;x++)for(int y=0;y<N;y++)for(int z=0;z<N;z++) {
                float m=mem[x][y][z][k],sample=image[x][y][z][c];
                if(k!=LESION_CLASS) {
                    float d=sample-centroid[c][k]; num+=(m*m+EPS)*d*d; den+=(m*m+EPS);
                } else {
                    float lm=lesion[x][y][z],bm=black[x][y][z];
                    float d=sample-centroid[c][k],dl=sample-lesionCentroid[c],db=sample-blackCentroid[c];
                    num+=(m*m+EPS)*d+(lm*lm+EPS)*dl+(bm*bm+EPS)*db;
                    den+=(m*m+EPS)+(lm*lm+EPS)+(bm*bm+EPS);
                }
            }
            channelWeight[c]=num>0.0f?den/num:1.0e10f; sum+=channelWeight[c];
        }
        for(int c=0;c<NC;c++) channelWeight[c]/=sum;
        return new Object[]{centroid, new float[][]{lesionCentroid,blackCentroid}, std, lesionStd, blackStd, channelWeight};
    }

    public static void main(String[] args) throws Exception {
        if(args.length!=1) throw new IllegalArgumentException("expected binary input path");
        readInput(args[0]); Object[] result=refresh();
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"wrapper_sha256\":\""+WRAPPER_SHA+"\",\"powertable_sha256\":\""+POWER_SHA+"\"");
        System.out.print(",\"class\":"+matrix((float[][])result[0]));
        System.out.print(",\"aux_centroids\":"+matrix((float[][])result[1]));
        System.out.print(",\"class_std\":"+vector(flatten((float[][])result[2])));
        System.out.print(",\"lesion_std\":"+vector((float[])result[3]));
        System.out.print(",\"blackhole_std\":"+vector((float[])result[4]));
        System.out.print(",\"weights\":"+vector((float[])result[5])+"}");
    }
    static float[] flatten(float[][] a) { float[] v=new float[a.length*a[0].length];int n=0;for(float[] row:a)for(float x:row)v[n++]=x;return v; }
}
