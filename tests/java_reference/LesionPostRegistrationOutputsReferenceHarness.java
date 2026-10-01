import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Selected-method reference for LesionTOADS terminal label/subtype updates
 * and the continuing-branch relation refresh, driven by test-generated state.
 * Source: VandyCRUISE LesionToadSegmentation.java SHA-256
 * 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 */
public final class LesionPostRegistrationOutputsReferenceHarness {
    private static final int MAX_RELATIONS=500;
    private static final int[][] BUILD_FACES={{-1,0,0},{0,-1,0},{0,0,-1},{0,0,1},{0,1,0},{1,0,0}};

    private static final class Input {
        int nx,ny,nz,classes,lesionClass,limit;
        byte[] templateLabels;
        byte[][][] segmentation,previousClassification;
        float[][][][] memberships;
        float[][][] lesion,blackHole;
    }

    private static Input read(String path) throws Exception {
        try(DataInputStream in=new DataInputStream(new FileInputStream(path))) {
            Input a=new Input();
            a.nx=in.readInt();a.ny=in.readInt();a.nz=in.readInt();
            a.classes=in.readInt();a.lesionClass=in.readInt();a.limit=in.readInt();
            a.templateLabels=new byte[a.classes];
            for(int k=0;k<a.classes;k++)a.templateLabels[k]=in.readByte();
            a.segmentation=new byte[a.nx][a.ny][a.nz];
            a.previousClassification=new byte[a.nx][a.ny][a.nz];
            for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++)a.segmentation[x][y][z]=in.readByte();
            for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++)a.previousClassification[x][y][z]=in.readByte();
            a.memberships=new float[a.nx][a.ny][a.nz][a.classes];
            for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++)for(int k=0;k<a.classes;k++)a.memberships[x][y][z][k]=in.readFloat();
            a.lesion=new float[a.nx][a.ny][a.nz];
            a.blackHole=new float[a.nx][a.ny][a.nz];
            for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++)a.lesion[x][y][z]=in.readFloat();
            for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++)a.blackHole[x][y][z]=in.readFloat();
            if(in.read()!=-1)throw new IllegalArgumentException("unexpected input trailing bytes");
            return a;
        }
    }

    private static void updateClassification(Input a, byte[][][] labels, short[][][] subtype) {
        for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++) {
            byte classLabel=a.templateLabels[a.segmentation[x][y][z]];
            labels[x][y][z]=classLabel;
            if(classLabel==a.templateLabels[a.lesionClass]
                    &&(a.lesion[x][y][z]>a.memberships[x][y][z][a.lesionClass]
                    ||a.blackHole[x][y][z]>a.memberships[x][y][z][a.lesionClass])) {
                subtype[x][y][z]=a.lesion[x][y][z]<a.blackHole[x][y][z]?(short)2:(short)1;
            }
        }
    }

    private static final class Relations {
        short[][][] map;
        boolean[][] rows;
    }

    private static Relations computeRelations(Input a) {
        short[][][] map=new short[a.nx][a.ny][a.nz];
        List<Set<Integer>> relationSets=new ArrayList<>();
        for(int k=0;k<a.classes;k++)relationSets.add(new LinkedHashSet<>(List.of(k)));
        for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++)map[x][y][z]=a.segmentation[x][y][z];
        for(int x=1;x<a.nx-1;x++)for(int y=1;y<a.ny-1;y++)for(int z=1;z<a.nz-1;z++) {
            Set<Integer> members=new LinkedHashSet<>();
            int center=a.segmentation[x][y][z];
            for(int[] d:BUILD_FACES) {
                int neighbor=a.segmentation[x+d[0]][y+d[1]][z+d[2]];
                if(neighbor!=center)members.add(neighbor);
            }
            if(!members.isEmpty()) {
                members.add(center);
                int id=-1;
                for(int r=0;r<relationSets.size();r++)if(relationSets.get(r).equals(members)){id=r;break;}
                if(id<0) {
                    if(relationSets.size()>=MAX_RELATIONS)throw new IllegalStateException("relation capacity");
                    id=relationSets.size();relationSets.add(new LinkedHashSet<>(members));
                }
                map[x][y][z]=(short)id;
            }
        }
        boolean[][] rows=new boolean[relationSets.size()][a.classes];
        for(int r=0;r<relationSets.size();r++)for(int k:relationSets.get(r))rows[r][k]=true;
        for(int k=0;k<a.classes;k++)for(int r=a.classes;r<relationSets.size();r++)if(rows[r][k])
            for(int m=0;m<a.classes;m++)if(rows[r][m])rows[k][m]=true;
        for(int t=1;t<5*a.limit;t++) {
            short[][][] tmp=new short[a.nx][a.ny][a.nz];
            for(int x=1;x<a.nx-1;x++)for(int y=1;y<a.ny-1;y++)for(int z=1;z<a.nz-1;z++) {
                if(map[x][y][z]<a.classes) {
                    if      (map[x-1][y][z]>=a.classes)tmp[x][y][z]=map[x-1][y][z];
                    else if (map[x+1][y][z]>=a.classes)tmp[x][y][z]=map[x+1][y][z];
                    else if (map[x][y-1][z]>=a.classes)tmp[x][y][z]=map[x][y-1][z];
                    else if (map[x][y+1][z]>=a.classes)tmp[x][y][z]=map[x][y+1][z];
                    else if (map[x][y][z-1]>=a.classes)tmp[x][y][z]=map[x][y][z-1];
                    else if (map[x][y][z+1]>=a.classes)tmp[x][y][z]=map[x][y][z+1];
                    else tmp[x][y][z]=map[x][y][z];
                } else tmp[x][y][z]=map[x][y][z];
            }
            for(int x=1;x<a.nx-1;x++)for(int y=1;y<a.ny-1;y++)for(int z=1;z<a.nz-1;z++)map[x][y][z]=tmp[x][y][z];
        }
        Relations out=new Relations();out.map=map;out.rows=rows;return out;
    }

    private static void appendInts(StringBuilder out,int[] values) {
        out.append('[');for(int i=0;i<values.length;i++){if(i>0)out.append(',');out.append(values[i]);}out.append(']');
    }
    private static void appendBool(StringBuilder out,boolean[][] values) {
        out.append('[');boolean first=true;for(boolean[] row:values)for(boolean value:row){if(!first)out.append(',');first=false;out.append(value);}out.append(']');
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=1)throw new IllegalArgumentException("expected state binary path");
        Input a=read(args[0]);
        byte[][][] labels=new byte[a.nx][a.ny][a.nz];
        short[][][] subtype=new short[a.nx][a.ny][a.nz];
        updateClassification(a,labels,subtype);
        Relations relations=computeRelations(a);
        int voxels=a.nx*a.ny*a.nz;
        int[] segmentation=new int[voxels],classification=new int[voxels],subtypes=new int[voxels],relationMap=new int[voxels];
        int i=0;for(int x=0;x<a.nx;x++)for(int y=0;y<a.ny;y++)for(int z=0;z<a.nz;z++) {
            segmentation[i]=a.segmentation[x][y][z];classification[i]=labels[x][y][z];
            subtypes[i]=subtype[x][y][z];relationMap[i]=relations.map[x][y][z];i++;
        }
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"segmentation\":");
        appendInts(out,segmentation);out.append(",\"classification\":");appendInts(out,classification);
        out.append(",\"lesion_classification\":");appendInts(out,subtypes);
        out.append(",\"relation_map\":");appendInts(out,relationMap);
        out.append(",\"relation_count\":").append(relations.rows.length).append(",\"relation_membership\":");appendBool(out,relations.rows);
        System.out.println(out.append('}'));
    }
}
