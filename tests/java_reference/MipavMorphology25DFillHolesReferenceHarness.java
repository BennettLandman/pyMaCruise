/* Test-only excerpt from MIPAV AlgorithmMorphology25D.FILL_HOLES and its
 * per-slice fillHolesRegion flood. Input/output are coordinate ordered x/y/z. */
import java.util.ArrayDeque;
import java.util.Queue;
public final class MipavMorphology25DFillHolesReferenceHarness {
    static final String SOURCE_SHA="d3e271151bce7e635777764648ef21c0253485699e577897c4ff4821f75c70f6";
    static final int X=8,Y=7,Z=3;
    static final int[][][] a=new int[X][Y][Z];
    static void init(){
        for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++)
            a[x][y][z]=((x*7+y*11+z*13)%9<5)?1:0;
        for(int x=1;x<X-1;x++)for(int y=1;y<Y-1;y++)a[x][y][0]=1;
        for(int x=2;x<6;x++)for(int y=2;y<5;y++)a[x][y][1]=1;
        a[3][3][1]=0; // enclosed hole
        a[4][3][1]=0; // 4-connected component
        a[1][1][2]=1; // isolated foreground
    }
    static void fill(){
        for(int z=0;z<Z;z++){
            for(int x=0;x<X;x++){a[x][0][z]=0;a[x][Y-1][z]=0;}
            for(int y=0;y<Y;y++){a[0][y][z]=0;a[X-1][y][z]=0;}
            boolean[][] seen=new boolean[X][Y];Queue<int[]> q=new ArrayDeque<>();
            q.add(new int[]{0,0});seen[0][0]=true;
            int[][] d={{1,0},{-1,0},{0,1},{0,-1}};
            while(!q.isEmpty()){
                int[] p=q.remove();int x=p[0],y=p[1];
                for(int[] v:d){int xx=x+v[0],yy=y+v[1];
                    if(xx>=0&&xx<X&&yy>=0&&yy<Y&&!seen[xx][yy]&&a[xx][yy][z]==0){seen[xx][yy]=true;q.add(new int[]{xx,yy});}
                }
            }
            for(int x=0;x<X;x++)for(int y=0;y<Y;y++)a[x][y][z]=seen[x][y]?0:1;
        }
    }
    public static void main(String[] args){
        init();fill();StringBuilder s=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA).append("\",\"data\":[");
        boolean first=true;for(int x=0;x<X;x++)for(int y=0;y<Y;y++)for(int z=0;z<Z;z++){if(!first)s.append(',');first=false;s.append(a[x][y][z]);}
        System.out.print(s.append("]}").toString());
    }
}
