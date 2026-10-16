import edu.jhmi.rad.medic.structures.RotationMatrix;
import java.security.MessageDigest;

/** Selected AtlasEMPreprocess crop/uncrop source excerpt on a small synthetic case. */
public final class AtlasEMCropReferenceHarness {
    private static final int X=13,Y=11,Z=9,C=2,BORDER=2;
    private static final float[] SPACING={0.8f,-1.1f,1.4f};
    private static final float[] INITIAL_TRANSFORM={0.12f,-0.08f,0.05f,0.4f,-0.3f,0.2f};
    private static final float BG_RATIO=0.15f, ILOW=0.1f, IHIGH=1.1f;

    private static float[][][][] images() {
        float[][][][] a=new float[C][X][Y][Z];
        for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
            a[0][x][y][z]=0.02f+((x*7+y*5+z*3)%13)*0.001f;
            a[1][x][y][z]=0.01f+((x*3+y*11+z*2)%9)*0.001f;
            if(x>=3&&x<=9&&y>=2&&y<=8&&z>=1&&z<=7)
                a[0][x][y][z]=0.4f+(x+y+z)*0.01f;
            if(x>=1&&x<=4&&y>=5&&y<=9&&z>=3&&z<=7)
                a[1][x][y][z]=0.35f+(x+y+z)*0.012f;
        }
        return a;
    }

    private static String digest(float[][][][] a) throws Exception {
        MessageDigest d=MessageDigest.getInstance("SHA-256");
        for(int x=0;x<a[0].length;x++) for(int y=0;y<a[0][x].length;y++) for(int z=0;z<a[0][x][y].length;z++) for(int c=0;c<a.length;c++) {
            int b=Float.floatToIntBits(a[c][x][y][z]);
            d.update((byte)(b>>>24));d.update((byte)(b>>>16));d.update((byte)(b>>>8));d.update((byte)b);
        }
        StringBuilder h=new StringBuilder();for(byte b:d.digest())h.append(String.format("%02x",b));return h.toString();
    }

    private static String bits(float[] a) {
        StringBuilder s=new StringBuilder("[");
        for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Integer.toUnsignedString(Float.floatToIntBits(a[i])));}
        return s.append(']').toString();
    }

    public static void main(String[] args) throws Exception {
        float[][][][] original=images();
        int[] lo={X-1,Y-1,Z-1},hi={0,0,0};
        float threshold=BG_RATIO*(IHIGH-ILOW);
        for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) for(int c=0;c<C;c++)
            if(original[c][x][y][z]>threshold){lo[0]=Math.min(lo[0],x);lo[1]=Math.min(lo[1],y);lo[2]=Math.min(lo[2],z);hi[0]=Math.max(hi[0],x);hi[1]=Math.max(hi[1],y);hi[2]=Math.max(hi[2],z);}
        int mx=hi[0]-lo[0]+1+2*BORDER,my=hi[1]-lo[1]+1+2*BORDER,mz=hi[2]-lo[2]+1+2*BORDER;
        float[][][][] cropped=new float[C][mx][my][mz];
        for(int c=0;c<C;c++) for(int x=lo[0];x<=hi[0];x++) for(int y=lo[1];y<=hi[1];y++) for(int z=lo[2];z<=hi[2];z++)
            cropped[c][x-lo[0]+BORDER][y-lo[1]+BORDER][z-lo[2]+BORDER]=original[c][x][y][z];

        float[] trans=INITIAL_TRANSFORM.clone();
        RotationMatrix rot=new RotationMatrix();rot.setParameters(trans[0],trans[1],trans[2]);float[][] r=rot.getMatrix();
        float dx=(mx/2.0f+lo[0]-BORDER-X/2.0f)*SPACING[0];
        float dy=(my/2.0f+lo[1]-BORDER-Y/2.0f)*SPACING[1];
        float dz=(mz/2.0f+lo[2]-BORDER-Z/2.0f)*SPACING[2];
        for(int a=0;a<3;a++) trans[a+3]+=r[a][0]*dx+r[a][1]*dy+r[a][2]*dz;
        float[] croppedTransform=trans.clone();

        float[][][][] uncropped=new float[C][X][Y][Z];
        for(int c=0;c<C;c++) for(int x=Math.max(lo[0]-BORDER,0);x<=Math.min(hi[0]+BORDER,X-1);x++)
            for(int y=Math.max(lo[1]-BORDER,0);y<=Math.min(hi[1]+BORDER,Y-1);y++)
                for(int z=Math.max(lo[2]-BORDER,0);z<=Math.min(hi[2]+BORDER,Z-1);z++)
                    uncropped[c][x][y][z]=cropped[c][x-lo[0]+BORDER][y-lo[1]+BORDER][z-lo[2]+BORDER];
        trans[3]-=r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        trans[4]-=r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        trans[5]-=r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;

        float[] registered=croppedTransform.clone();
        registered[3]+=0.17f;registered[4]-=0.09f;registered[5]+=0.04f;
        float[][][][] postRegistrationNative=new float[C][X][Y][Z];
        for(int c=0;c<C;c++) for(int x=Math.max(lo[0]-BORDER,0);x<=Math.min(hi[0]+BORDER,X-1);x++)
            for(int y=Math.max(lo[1]-BORDER,0);y<=Math.min(hi[1]+BORDER,Y-1);y++)
                for(int z=Math.max(lo[2]-BORDER,0);z<=Math.min(hi[2]+BORDER,Z-1);z++)
                    postRegistrationNative[c][x][y][z]=cropped[c][x-lo[0]+BORDER][y-lo[1]+BORDER][z-lo[2]+BORDER];
        registered[3]-=r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        registered[4]-=r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        registered[5]-=r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;
        int[] lo2={X-1,Y-1,Z-1},hi2={0,0,0};
        for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) for(int c=0;c<C;c++)
            if(postRegistrationNative[c][x][y][z]>threshold){lo2[0]=Math.min(lo2[0],x);lo2[1]=Math.min(lo2[1],y);lo2[2]=Math.min(lo2[2],z);hi2[0]=Math.max(hi2[0],x);hi2[1]=Math.max(hi2[1],y);hi2[2]=Math.max(hi2[2],z);}
        int mx2=hi2[0]-lo2[0]+1+2*BORDER,my2=hi2[1]-lo2[1]+1+2*BORDER,mz2=hi2[2]-lo2[2]+1+2*BORDER;
        float[][][][] recropped=new float[C][mx2][my2][mz2];
        for(int c=0;c<C;c++) for(int x=lo2[0];x<=hi2[0];x++) for(int y=lo2[1];y<=hi2[1];y++) for(int z=lo2[2];z<=hi2[2];z++)
            recropped[c][x-lo2[0]+BORDER][y-lo2[1]+BORDER][z-lo2[2]+BORDER]=postRegistrationNative[c][x][y][z];
        float d2x=(mx2/2.0f+lo2[0]-BORDER-X/2.0f)*SPACING[0];
        float d2y=(my2/2.0f+lo2[1]-BORDER-Y/2.0f)*SPACING[1];
        float d2z=(mz2/2.0f+lo2[2]-BORDER-Z/2.0f)*SPACING[2];
        for(int a=0;a<3;a++) registered[a+3]+=r[a][0]*d2x+r[a][1]*d2y+r[a][2]*d2z;

        System.out.printf("{\"atlas_em_preprocess_sha256\":\"acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b\",\"rotation_matrix_sha256\":\"144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0\",\"lower\":[%d,%d,%d],\"upper\":[%d,%d,%d],\"cropped_dims\":[%d,%d,%d],\"crop_sha256\":\"%s\",\"uncrop_sha256\":\"%s\",\"crop_transform_bits\":%s,\"restored_transform_bits\":%s,\"recrop_lower\":[%d,%d,%d],\"recrop_upper\":[%d,%d,%d],\"recrop_dims\":[%d,%d,%d],\"recrop_sha256\":\"%s\",\"registered_transform_bits\":%s}",
                lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],mx,my,mz,digest(cropped),digest(uncropped),bits(croppedTransform),bits(trans),lo2[0],lo2[1],lo2[2],hi2[0],hi2[1],hi2[2],mx2,my2,mz2,digest(recropped),bits(registered));
    }
}
