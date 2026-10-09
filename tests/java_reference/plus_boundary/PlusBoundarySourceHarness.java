package edu.vanderbilt.masi.plugins.CRUISE.utilities;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/** Executes the preserved MinFilter3D source on a held-fixed label volume. */
public final class PlusBoundarySourceHarness {
    private PlusBoundarySourceHarness() {}

    public static void main(String[] args) throws Exception {
        try (DataInputStream input = new DataInputStream(new FileInputStream(args[0]))) {
            int nx=input.readInt(), ny=input.readInt(), nz=input.readInt();
            int[][][] labels=new int[nx][ny][nz];
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)labels[x][y][z]=input.readInt();
            if(input.read()!=-1)throw new IllegalArgumentException("trailing input bytes");
            int[] overlooked={120,121,202,203,156,157};
            for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)
                for(int label:overlooked)if(labels[x][y][z]==label){labels[x][y][z]=0;break;}
            Set<Integer> unique=new HashSet<>();
            for(int[][] plane:labels)for(int[] row:plane)for(int value:row)unique.add(value);
            int[][][] boundary=new int[nx][ny][nz];
            MinFilter3D source=new MinFilter3D();
            Integer[] ordered=unique.toArray(new Integer[0]);Arrays.sort(ordered);
            for(int label:ordered)if(label>=100)boundary=source.GetBoundrayForOneLabel(boundary,labels,label,0);
            try(DataOutputStream output=new DataOutputStream(new FileOutputStream(args[1]))) {
                output.writeInt(nx);output.writeInt(ny);output.writeInt(nz);
                for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++)output.writeInt(boundary[x][y][z]);
            }
        }
    }
}
