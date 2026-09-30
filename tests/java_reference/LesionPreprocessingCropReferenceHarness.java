/*
 * Selected VandyCRUISE LesionToadPreprocess.findCroppingBoundaries() and
 * cropImages() excerpt. Source SHA-256:
 * 72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08.
 * The active wrapper calls normalizeImages first and constructs the preprocessor
 * with bgLabel=1, bgRatio=0, border=10. This RAS fixture uses axial identity
 * geometry; image and atlas supports are asymmetric and the topology expands
 * the intensity-only crop. Normalization itself is excluded.
 */
public final class LesionPreprocessingCropReferenceHarness {
    static final int NX=11,NY=10,NZ=9,NC=2,BS=10,BG=1;
    static final float[] RES={0.8f,1.2f,1.5f};
    static float image(int c,int x,int y,int z) {
        if(!(x>=3&&x<=7&&y>=2&&y<=7&&z>=2&&z<=6)||((x+2*y+3*z+c)%6)==0)return 0.0f;
        return ((x*11+y*7+z*3+c*5)%17+1)/20.0f;
    }
    static byte topology(int x,int y,int z) {
        if(x>=7&&x<10&&y>=1&&y<6&&z>=3&&z<8)return 3;
        if(x>=1&&x<5&&y>=3&&y<8&&z>=1&&z<6)return 2;
        return 1;
    }
    static int bits(float v){return Float.floatToIntBits(v);}
    static int[] bits(float[] values){int[] result=new int[values.length];for(int i=0;i<values.length;i++)result[i]=bits(values[i]);return result;}
    static void append(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    public static void main(String[] args) {
        int x0=NX-1,xN=0,y0=NY-1,yN=0,z0=NZ-1,zN=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<NC;c++) {
            float value=image(c,x,y,z);
            if(value>0.0f||value>0.0f){x0=Math.min(x0,x);xN=Math.max(xN,x);y0=Math.min(y0,y);yN=Math.max(yN,y);z0=Math.min(z0,z);zN=Math.max(zN,z);}
        }
        // Axial RAS identity: computeImageCoordinates maps topology indices to themselves.
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(topology(x,y,z)>BG) {
            x0=Math.min(x0,x);xN=Math.max(xN,x);y0=Math.min(y0,y);yN=Math.max(yN,y);z0=Math.min(z0,z);zN=Math.max(zN,z);
        }
        int mx=xN-x0+1+2*BS,my=yN-y0+1+2*BS,mz=zN-z0+1+2*BS;
        float[] transform=new float[6];
        transform[3]=(mx/2.0f+x0-BS-NX/2.0f)*RES[0];
        transform[4]=(my/2.0f+y0-BS-NY/2.0f)*RES[1];
        transform[5]=(mz/2.0f+z0-BS-NZ/2.0f)*RES[2];
        int[] cropped=new int[NC*mx*my*mz];
        for(int c=0;c<NC;c++)for(int x=x0;x<=xN;x++)for(int y=y0;y<=yN;y++)for(int z=z0;z<=zN;z++) {
            int a=x-x0+BS,b=y-y0+BS,d=z-z0+BS;
            float value=(x<0||x>=NX||y<0||y>=NY||z<0||z>=NZ)?0.0f:image(c,x,y,z);
            cropped[((c*mx+a)*my+b)*mz+d]=bits(value);
        }
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08\",\"dims_xyz\":["+NX+","+NY+","+NZ+"],\"channels\":"+NC+",\"bounds_inclusive\":["+x0+","+xN+","+y0+","+yN+","+z0+","+zN+"],\"padding\":"+BS+",\"background_ratio_bits\":"+bits(0.0f)+",\"cropped_dims_xyz\":["+mx+","+my+","+mz+"],\"cropped_transform_bits\":");
        append(out,bits(transform));out.append(",\"cropped_image_bits\":");append(out,cropped);
        out.append(",\"call_order\":[\"normalizeImages\",\"findCroppingBoundaries\",\"cropImages\"]}");
        System.out.println(out.toString());
    }
}
