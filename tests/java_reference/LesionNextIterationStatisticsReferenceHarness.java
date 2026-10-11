import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * Fixed-profile excerpt of VandyCRUISE LesionToadSegmentation's selected
 * computeCentroids, computeVariances, and computeChannelWeights methods.
 *
 * Contract: one binary input path and one JSON output path. Integers and
 * floats are big-endian Data{Input,Output}Stream values. The input begins
 * nx,ny,nz,nc,classes, then class type bytes (0 MASK, 1 OPTIMIZED_CSF,
 * 2 OBJECT, 3 LESION); RAS float image[x,y,z,c]; float memberships[x,y,z,k];
 * int segmentation[x,y,z]; lesion, black-hole, and outlier float maps;
 * initial centroids[c,k], lesion/black-hole centroids[c], tissue powers[c,k],
 * and lesion powers[c]. Fields are absent; centroidMode is "none";
 * normType is "Identity"; useLesionWeight is false.
 *
 * The implementation intentionally uses source float accumulation and
 * preserves the source's finite INF sentinel and reports all IEEE NaNs.
 */
public final class LesionNextIterationStatisticsReferenceHarness {
    private static final float ZERO = 1.0e-10f;
    private static final float EPS = 1.0e-30f;
    private static final float INF = 1.0e30f;
    private static final int MASK=0, OPT_CSF=1, OBJECT=2, LESION=3;
    private int nx,ny,nz,nc,classes;
    private byte[] kinds;
    private float[][][][] image, mem;
    private int[][][] seg;
    private float[][][] lesion, black, outlier;
    private float[][] centroid, initialCentroid, classStd, tissuePower;
    private float[] lesionCentroid, blackCentroid, lesionStd, blackStd, lesionPower, channelWeight;

    private void read(String path) throws Exception {
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(path)))) {
            nx=in.readInt(); ny=in.readInt(); nz=in.readInt(); nc=in.readInt(); classes=in.readInt();
            if(nx<1||ny<1||nz<1||nc<1||classes<1) throw new IllegalArgumentException("invalid dimensions");
            kinds=new byte[classes]; in.readFully(kinds);
            image=new float[nx][ny][nz][nc]; mem=new float[nx][ny][nz][classes]; seg=new int[nx][ny][nz];
            lesion=new float[nx][ny][nz]; black=new float[nx][ny][nz]; outlier=new float[nx][ny][nz];
            initialCentroid=new float[nc][classes]; centroid=new float[nc][classes];
            lesionCentroid=new float[nc]; blackCentroid=new float[nc];
            tissuePower=new float[nc][classes]; lesionPower=new float[nc];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int c=0;c<nc;c++)image[x][y][z][c]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)for(int k=0;k<classes;k++)mem[x][y][z][k]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)seg[x][y][z]=in.readInt();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)lesion[x][y][z]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)black[x][y][z]=in.readFloat();
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)outlier[x][y][z]=in.readFloat();
            for(int c=0;c<nc;c++)for(int k=0;k<classes;k++)initialCentroid[c][k]=centroid[c][k]=in.readFloat();
            for(int c=0;c<nc;c++)lesionCentroid[c]=in.readFloat();
            for(int c=0;c<nc;c++)blackCentroid[c]=in.readFloat();
            for(int c=0;c<nc;c++)for(int k=0;k<classes;k++)tissuePower[c][k]=in.readFloat();
            for(int c=0;c<nc;c++)lesionPower[c]=in.readFloat();
            if(in.read()!=-1) throw new IllegalArgumentException("trailing input bytes");
        }
        classStd=new float[nc][classes]; lesionStd=new float[nc]; blackStd=new float[nc]; channelWeight=new float[nc];
    }

    // Source PowerTable(0,4,1e-6,q).lookup(p,q), reduced to a direct lookup
    // for the fixture's bounded positive p values; cached table entries are
    // Math.pow(index*step,q), then divided by p by source callers.
    private static double powerLookup(float p, float exponent) {
        if(p<0.0f||p>4.0f) return Math.pow(p, exponent);
        // Source constructor receives 1e-6f and stores that promoted binary32
        // value as double, matching Python's float(np.float32(1e-6)).
        double step=(double)1.0e-6f;
        int index=(int)((((double)p)/step)+0.5);
        double quantized=index*step;
        return Math.pow(quantized, (double)exponent);
    }
    private static float factor(float distance,float exponent) {
        return distance>ZERO ? (float)(powerLookup(distance,exponent)/(double)distance) : 1.0f;
    }
    private float sample(int x,int y,int z,int c) { return image[x][y][z][c]; }
    private float winner(int x,int y,int z,int k) {
        return Math.max(mem[x][y][z][k],Math.max(lesion[x][y][z],black[x][y][z]));
    }

    private void computeCentroids() {
        for(int c=0;c<nc;c++) {
            float[] num=new float[classes], den=new float[classes];
            float numL=0,denL=0,numB=0,denB=0;
            for(int k=0;k<classes;k++) for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) if(seg[x][y][z]==k) {
                float m=mem[x][y][z][k], s=sample(x,y,z,c), kind=kinds[k];
                if(kind==LESION) {
                    float max=winner(x,y,z,k);
                    if(max==m) { num[k]+=m*m*s; den[k]+=m*m; }
                    else if(max==lesion[x][y][z]) {
                        float d=(s-lesionCentroid[c])*(s-lesionCentroid[c]);
                        float f=c==nc-1?factor(d,lesionPower[c]):1.0f;
                        numL+=f*lesion[x][y][z]*lesion[x][y][z]*s;
                        denL+=f*lesion[x][y][z]*lesion[x][y][z];
                    } else { numB+=black[x][y][z]*black[x][y][z]*s; denB+=black[x][y][z]*black[x][y][z]; }
                } else if(kind==OPT_CSF) {
                    float d=s-centroid[c][k]; d*=d;
                    float f=factor(d,tissuePower[c][k]);
                    num[k]+=f*m*m*s; den[k]+=f*m*m;
                } else {
                    num[k]+=m*m*s; den[k]+=m*m;
                }
            }
            for(int k=0;k<classes;k++) centroid[c][k]=den[k]>0.0f?num[k]/den[k]:0.0f;
            lesionCentroid[c]=denL>0.0f?numL/denL:0.0f;
            blackCentroid[c]=denB>0.0f?numB/denB:0.0f;
        }
    }
    private void computeVariances() {
        for(int c=0;c<nc;c++) {
            float lnum=0,lden=0,bnum=0,bden=0;
            for(int k=0;k<classes;k++) {
                float num=0,den=0;
                for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)if(seg[x][y][z]==k) {
                    float m=mem[x][y][z][k], lm=lesion[x][y][z], bm=black[x][y][z], s=sample(x,y,z,c);
                    if((kinds[k]==LESION&&m<lm)||(false&&m<bm)) {
                        if(lm>=bm) { float d=s-lesionCentroid[c], w=lm*lm+EPS; lnum+=w*d*d; lden+=w; }
                        else { float d=s-blackCentroid[c], w=lm*lm+EPS; bnum+=w*d*d; bden+=w; }
                    } else { float d=s-centroid[c][k], w=m*m+EPS; num+=w*d*d; den+=w; }
                }
                classStd[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
            }
            lesionStd[c]=lden>0.0f?(float)Math.sqrt(lnum/lden):0.0f;
            blackStd[c]=bden>0.0f?(float)Math.sqrt(bnum/bden):0.0f;
        }
    }
    private void computeChannelWeights() {
        // Selected source computeChannelWeights() returns uniform weights for
        // every norm type other than "Diagonal Covariance".
        for(int c=0;c<nc;c++) channelWeight[c]=1.0f/nc;
    }
    private static long bits(float f) { return Integer.toUnsignedLong(Float.floatToRawIntBits(f)); }
    private static String vector(float[] a) { StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString(); }
    private static String matrix(float[][] a) { StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(vector(a[i]));}return s.append(']').toString(); }
    private static int nanCount(float[][] a) { int n=0;for(float[] r:a)for(float v:r)if(Float.isNaN(v))n++;return n; }
    private static int nanCount(float[] a) { int n=0;for(float v:a)if(Float.isNaN(v))n++;return n; }
    private void write(String path) throws Exception {
        try(DataOutputStream out=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)))) {
            out.writeBytes("{\"class_centroid_bits\":"+matrix(centroid)+",\"lesion_centroid_bits\":"+vector(lesionCentroid)+",\"black_hole_centroid_bits\":"+vector(blackCentroid)+",\"class_std_bits\":"+matrix(classStd)+",\"lesion_std_bits\":"+vector(lesionStd)+",\"black_hole_std_bits\":"+vector(blackStd)+",\"channel_weight_bits\":"+vector(channelWeight)+",\"nan_count\":"+(nanCount(centroid)+nanCount(lesionCentroid)+nanCount(blackCentroid)+nanCount(classStd)+nanCount(lesionStd)+nanCount(blackStd)+nanCount(channelWeight))+"}");
        }
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=2) throw new IllegalArgumentException("input.bin output.json");
        LesionNextIterationStatisticsReferenceHarness h=new LesionNextIterationStatisticsReferenceHarness();
        h.read(args[0]); h.computeCentroids(); h.computeVariances(); h.computeChannelWeights(); h.write(args[1]);
    }
}
