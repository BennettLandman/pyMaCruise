import java.util.Random;

/** Small-grid formula excerpt for one selected propagateThinning queue item. */
public final class ToadsThinningTransitionReferenceHarness {
    private static final int C=4, EMPTY=-1;
    private static final int[][] FACE={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};

    private static final class Case {
        int[][][] seg=new int[3][3][3], skel=new int[3][3][3];
        float[] mem={0.1f,0.8f,0.05f,0.05f};
        boolean[][] relation=new boolean[C][C];
        float[][][][] shape=new float[3][3][3][C];
        float[] volume={200.0f,210.0f,220.0f,230.0f}, minimum={100.0f,100.0f,100.0f,100.0f};
        int old=0, target=1, candidate=-1, proposed=-1;
        int[][][] inputSeg,inputSkel;
        float[] inputVolume;
        boolean topology=false, topologyChecked=false, eligible=true;
        Case(Random r) {
            for(int x=0;x<3;x++)for(int y=0;y<3;y++)for(int z=0;z<3;z++){
                seg[x][y][z]=r.nextInt(C);skel[x][y][z]=seg[x][y][z];
                for(int k=0;k<C;k++)shape[x][y][z][k]=(1+r.nextInt(9))/10.0f;
            }
            seg[1][1][1]=old;skel[1][1][1]=old;
        }
        void calculate() {
            if(skel[1][1][1]!=old){eligible=false;return;}
            int id=0;float max=0.0f;for(int k=0;k<C;k++)if(mem[k]>max){id=k;max=mem[k];}
            int next=-1;float bestPrior=-0.1f;boolean direct=false;
            if(id!=old) for(int[] d:FACE){int nx=1+d[0],ny=1+d[1],nz=1+d[2],label=seg[nx][ny][nz];
                if(label==id)direct=true;
                else if(label!=old&&relation[label][id]) {float prior=shape[nx][ny][nz][id];if(prior>bestPrior){next=label;bestPrior=prior;}}
            }
            boolean change=id!=old&&(direct||next!=-1);
            candidate=id;if(!direct&&next!=-1)candidate=next;
            proposed=change?candidate:-1;
            topology=false;topologyChecked=change;
            if(change) topology=ToadsHomeomorphicChangeReferenceHarness.isHomeomorphic(toShort(seg),(short)1,(short)1,(short)1,(short)candidate);
            eligible=skel[1][1][1]==old&&volume[old]>minimum[old];
            if(!eligible) {candidate=-1;return;}
            volume[old]-=mem[old];skel[1][1][1]=EMPTY;
            if(change&&topology)seg[1][1][1]=candidate;else candidate=-1;
        }
    }
    private static int[][][] copy(int[][][] a){int[][][] b=new int[3][3][3];for(int x=0;x<3;x++)for(int y=0;y<3;y++)for(int z=0;z<3;z++)b[x][y][z]=a[x][y][z];return b;}
    private static short[][][] toShort(int[][][] a){short[][][] b=new short[3][3][3];for(int x=0;x<3;x++)for(int y=0;y<3;y++)for(int z=0;z<3;z++)b[x][y][z]=(short)a[x][y][z];return b;}
    private static String labels(int[][][] a){StringBuilder b=new StringBuilder();for(int[][] x:a)for(int[] y:x)for(int v:y)b.append((char)('0'+v));return b.toString();}
    private static String floats(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(Float.toString(a[i]));}return b.append(']').toString();}
    private static String bools(boolean[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(a[i]);}return b.append(']').toString();}
    private static String matrix(boolean[][] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bools(a[i]));}return b.append(']').toString();}
    private static String shapes(float[][][][] a){StringBuilder b=new StringBuilder("[");boolean first=true;for(int x=0;x<3;x++)for(int y=0;y<3;y++)for(int z=0;z<3;z++)for(int k=0;k<C;k++){if(!first)b.append(',');first=false;b.append(Float.toString(a[x][y][z][k]));}return b.append(']').toString();}
    private static void fill(Case c,int label){for(int x=0;x<3;x++)for(int y=0;y<3;y++)for(int z=0;z<3;z++){c.seg[x][y][z]=label;c.skel[x][y][z]=label;}}
    private static void decode(Case c,String text){for(int i=0;i<27;i++){int x=i/9,y=(i/3)%3,z=i%3;c.seg[x][y][z]=text.charAt(i)-'0';c.skel[x][y][z]=c.seg[x][y][z];}}
    private static Case make(String kind){
        Case c=new Case(new Random(0));fill(c,0);
        if(kind.equals("same")){c.mem=new float[]{0.8f,0.1f,0.05f,0.05f};}
        else if(kind.equals("low_volume")){c.seg[0][1][1]=1;c.volume[0]=100.0f;c.minimum[0]=100.0f;}
        else if(kind.equals("not_skeleton")){c.skel[1][1][1]=1;c.seg[0][1][1]=1;}
        else if(kind.equals("no_connection")){c.relation=new boolean[C][C];}
        else if(kind.equals("direct_allowed")){c.seg[0][1][1]=1;}
        else if(kind.equals("direct_rejected")){decode(c,"020033121011213000023100323");c.old=1;c.target=3;c.mem=new float[]{0.1f,0.1f,0.1f,0.8f};}
        else if(kind.equals("related_allowed")){c.seg[0][1][1]=2;c.relation[2][1]=true;c.shape[0][1][1][1]=0.7f;}
        else if(kind.equals("related_rejected")){decode(c,"020323122333000001220101021");c.mem=new float[]{0.1f,0.8f,0.05f,0.05f};c.relation[2][1]=true;c.relation[3][1]=true;c.shape[0][1][1][1]=0.7f;c.shape[1][0][1][1]=0.3f;}
        else if(kind.equals("related_tie")){c.seg[0][1][1]=2;c.seg[1][0][1]=3;c.relation[2][1]=c.relation[3][1]=true;c.shape[0][1][1][1]=0.5f;c.shape[1][0][1][1]=0.5f;}
        c.inputSeg=copy(c.seg);c.inputSkel=copy(c.skel);c.inputVolume=c.volume.clone();
        c.calculate();return c;
    }
    private static void emit(String kind,Case c,boolean first){
        if(!first)System.out.print(',');
        boolean[] available={false,true,true,true};
        boolean[] outAvailable=available.clone();
        if(c.eligible)outAvailable[c.old]=false;
        System.out.print("{\"kind\":\""+kind+"\",\"segmentation\":\""+labels(c.inputSeg)+"\",\"skeleton\":\""+labels(c.inputSkel)+"\",\"memberships\":"+floats(c.mem)+",\"relations\":"+matrix(c.relation)+",\"shape\":"+shapes(c.shape)+",\"volumes\":"+floats(c.inputVolume)+",\"minimums\":"+floats(c.minimum)+",\"available\":"+bools(available)+",\"expected\":{\"segmentation_label\":"+c.seg[1][1][1]+",\"skeleton_label\":"+c.skel[1][1][1]+",\"volumes\":"+floats(c.volume)+",\"available\":"+bools(outAvailable)+",\"candidate_label\":"+(c.proposed<0?"null":c.proposed)+",\"topology_accepted\":"+(c.topologyChecked?(c.topology?"true":"false"):"null")+",\"eligible\":"+c.eligible+"}}");
    }
    public static void main(String[] args){
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        String[] kinds={"same","low_volume","not_skeleton","no_connection","direct_allowed","direct_rejected","related_allowed","related_rejected","related_tie"};
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.propagateThinning local queue item\",\"source_sha256\":\"63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"lut_source_sha256\":\"50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f\",\"seed\":0,\"cases\":[");
        for(int i=0;i<kinds.length;i++){Case c=make(kinds[i]);emit(kinds[i],c,i==0);}
        System.out.println("]}");
    }
}
