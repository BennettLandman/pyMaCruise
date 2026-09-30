import edu.jhmi.rad.medic.utilities.PowerTable;

/** Fixed-input reference for the selected post-thinning statistics call order. */
public final class LesionStatisticsRefreshReferenceHarness {
    static final int NX=6, NY=5, NZ=4, NC=2, K=5;
    static final String SOURCE_SHA="545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    static final String WRAPPER_SHA="556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb";
    static final String POWER_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    static final String[] TYPE={"object","optimized_gm","lesion","blackhole","object"};
    static final float[][][][] image=new float[NX][NY][NZ][NC];
    static final float[][][][] field=new float[NX][NY][NZ][NC];
    static final float[][][][] mem=new float[NX][NY][NZ][K];
    static final int[][][] seg=new int[NX][NY][NZ];
    static final float[][][] outlier=new float[NX][NY][NZ], lesion=new float[NX][NY][NZ], black=new float[NX][NY][NZ];
    static final float[][] initial=new float[NC][K];
    static final float[] initialLesion={0.85f,0.95f}, initialBlack={0.45f,0.52f};
    static final float[] gmPower={0.9f,1.05f}, lesionPower={1.2f,0.8f};
    static final float[] lesionNumerator=new float[NC], lesionDenominator=new float[NC];

    static int bits(float x) { return Float.floatToRawIntBits(x); }

    static void makeData() {
        for (int c=0;c<NC;c++) for (int k=0;k<K;k++) initial[c][k]=0.2f+0.15f*k+0.11f*c;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            seg[x][y][z]=(x+2*y+z+((x*y+2*z)%3))%4; // class 4 intentionally absent
            for (int k=0;k<K;k++) mem[x][y][z][k]=0.08f+((x*3+y*5+z*7+k*2)%9)*0.055f;
            for (int c=0;c<NC;c++) {
                image[x][y][z][c]=0.25f+((x*11+y*7+z*3+c*5)%23)*0.12f;
                field[x][y][z][c]=0.75f+((x*2+y*3+z*4+c)%5)*0.11f;
            }
            outlier[x][y][z]=((x+3*y+z)%5)*0.13f;
            lesion[x][y][z]=0.04f+((x+2*y+z)%9)*0.07f;
            black[x][y][z]=0.03f+((2*x+y+3*z)%9)*0.075f;
        }
        int[][][] points=new int[4][2][3]; int[] count=new int[4];
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            int k=seg[x][y][z]; if(count[k]<2) points[k][count[k]++]=new int[]{x,y,z};
        }
        setAux(points[2][0],2,0.30f,0.80f,0.20f);
        setAux(points[2][1],2,0.25f,0.70f,0.70f);
        setAux(points[3][0],3,0.25f,0.20f,0.90f);
    }

    static void setAux(int[] p,int k,float m,float l,float b) {
        int x=p[0],y=p[1],z=p[2]; mem[x][y][z][k]=m; lesion[x][y][z]=l; black[x][y][z]=b;
    }

    static float[][][] refresh(boolean useField) {
        float[][] centroid=new float[NC][K]; for(int c=0;c<NC;c++) centroid[c]=initial[c].clone();
        float[] lesionCentroid=initialLesion.clone(), blackCentroid=initialBlack.clone();
        PowerTable[] gmTables=new PowerTable[NC], lesionTables=new PowerTable[NC];
        for(int c=0;c<NC;c++) {
            gmTables[c]=new PowerTable(0.0f,4.0f,0.000001f,gmPower[c]);
            lesionTables[c]=new PowerTable(0.0f,4.0f,0.000001f,lesionPower[c]);
        }

        // Selected computeCentroids body excerpt (wrapper invokes this first).
        for (int c=0;c<NC;c++) {
            float[] num=new float[K], den=new float[K]; float numL=0,denL=0,numB=0,denB=0;
            for (int k=0;k<K;k++) for(int x=0;x<NX;x++) for(int y=0;y<NY;y++) for(int z=0;z<NZ;z++) if(seg[x][y][z]==k) {
                float m=mem[x][y][z][k], sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                if(TYPE[k].equals("object")) { num[k]+=m*m*sample; den[k]+=m*m; }
                else if(TYPE[k].equals("optimized_gm")) {
                    float d=sample-initial[c][k]; d*=d; float factor=1.0f;
                    if(d>1e-10f) factor=(float)(gmTables[c].lookup((double)d,(double)gmPower[c])/d);
                    num[k]+=factor*m*m*sample; den[k]+=factor*m*m;
                } else if(TYPE[k].equals("lesion")||TYPE[k].equals("blackhole")) {
                    float lm=lesion[x][y][z], bm=black[x][y][z], winner=Math.max(m,Math.max(lm,bm));
                    float factor=1.0f;
                    // useLesionWeight defaults false: selected source weights
                    // unfielded FLAIR lesion means, but not fielded means.
                    if(winner!=m && winner==lm && c==1 && !useField) {
                        float d=sample-initialLesion[c]; d*=d;
                        if(d>1e-10f) factor=(float)(lesionTables[c].lookup((double)d,(double)lesionPower[c])/d);
                    }
                    if(winner==m) { num[k]+=factor*m*m*sample; den[k]+=factor*m*m; }
                    else if(winner==lm) { numL+=factor*lm*lm*sample; denL+=factor*lm*lm; }
                    else { numB+=factor*bm*bm*sample; denB+=factor*bm*bm; }
                }
            }
            for(int k=0;k<K;k++) centroid[c][k]=den[k]>0.0f?num[k]/den[k]:0.0f;
            lesionNumerator[c]=numL; lesionDenominator[c]=denL;
            lesionCentroid[c]=denL>0.0f?numL/denL:0.0f;
            blackCentroid[c]=denB>0.0f?numB/denB:0.0f;
        }

        // Selected computeVariances body excerpt (consumes updated centroids).
        float[][] std=new float[NC][K]; float[] lesionStd=new float[NC], blackStd=new float[NC];
        for(int c=0;c<NC;c++) {
            float lnum=0,lden=0,bnum=0,bden=0;
            for(int k=0;k<K;k++) {
                float num=0,den=0;
                for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) if(seg[x][y][z]==k) {
                    float m=mem[x][y][z][k], lm=lesion[x][y][z], bm=black[x][y][z];
                    float sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                    if((TYPE[k].equals("lesion")&&m<lm)||(TYPE[k].equals("blackhole")&&m<bm)) {
                        if(lm>=bm) { float d=sample-lesionCentroid[c]; float w=lm*lm+1e-30f; lnum+=w*d*d; lden+=w; }
                        else { float d=sample-blackCentroid[c]; float w=lm*lm+1e-30f; bnum+=w*d*d; bden+=w; }
                    } else {
                        float d=sample-centroid[c][k]; float w=m*m+1e-30f; num+=w*d*d; den+=w;
                    }
                }
                std[c][k]=den>0.0f?(float)Math.sqrt(num/den):0.0f;
            }
            lesionStd[c]=lden>0?(float)Math.sqrt(lnum/lden):0.0f;
            blackStd[c]=bden>0?(float)Math.sqrt(bnum/bden):0.0f;
        }

        // Selected computeChannelWeights body excerpt (consumes updated centroids).
        float[] channelWeight=new float[NC]; float sum=0;
        for(int c=0;c<NC;c++) {
            float num=0,den=0;
            for(int k=0;k<K;k++) for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++) {
                float m=mem[x][y][z][k], sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                if(TYPE[k].equals("object")) {
                    float d=sample-centroid[c][k]; num+=(m*m+1e-30f)*d*d; den+=(m*m+1e-30f);
                } else if(TYPE[k].equals("lesion")||TYPE[k].equals("blackhole")) {
                    float lm=lesion[x][y][z], bm=black[x][y][z];
                    float d=sample-centroid[c][k], dl=sample-lesionCentroid[c], db=sample-blackCentroid[c];
                    num+=(m*m+1e-30f)*d+(lm*lm+1e-30f)*dl+(bm*bm+1e-30f)*db;
                    den+=(m*m+1e-30f)+(lm*lm+1e-30f)+(bm*bm+1e-30f);
                }
            }
            channelWeight[c]=num>0.0f?den/num:1e10f; sum+=channelWeight[c];
        }
        for(int c=0;c<NC;c++) channelWeight[c]/=sum;
        return new float[][][]{centroid,{lesionCentroid,blackCentroid},{flatten(std),lesionStd,blackStd},{channelWeight}};
    }

    static float[] flatten(float[][] a) { float[] v=new float[a.length*a[0].length];int n=0;for(float[] row:a)for(float x:row)v[n++]=x;return v; }
    static String vector(float[] a) { StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString(); }
    static String matrix(float[][] a) { StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(vector(a[i]));}return s.append(']').toString(); }

    public static void main(String[] args) {
        makeData();
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"wrapper_sha256\":\""+WRAPPER_SHA+"\",\"powertable_sha256\":\""+POWER_SHA+"\"");
        for(int f=0;f<2;f++) {
            float[][][] o=refresh(f==1); String key=f==1?"with_field":"no_field";
            System.out.print(",\""+key+"\":{\"class\":"+matrix(o[0])+",\"aux_centroids\":"+matrix(o[1])+",\"lesion_num_den\":"+matrix(new float[][]{lesionNumerator,lesionDenominator})+",\"class_std\":"+vector(o[2][0])+",\"lesion_std\":"+vector(o[2][1])+",\"blackhole_std\":"+vector(o[2][2])+",\"weights\":"+vector(o[3][0])+"}");
        }
        float[] fallback=new float[NC];for(int c=0;c<NC;c++)fallback[c]=(float)1/NC;
        System.out.print(",\"uniform_weight_fallback\":"+vector(fallback));
        System.out.print("}");
    }
}
