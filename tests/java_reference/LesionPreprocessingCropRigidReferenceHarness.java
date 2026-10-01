/*
 * Fixed-input excerpt of VandyCRUISE LesionToadPreprocess crop methods.
 * LesionToadPreprocess.java SHA-256:
 * 72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08
 * Nonzero Rodrigues rotation/translation exercises computeImageCoordinates
 * before inclusive bounds and the rotation-aware crop translation update.
 */
public final class LesionPreprocessingCropRigidReferenceHarness {
    static final int NX=11, NY=10, NZ=9, NC=2, BS=10, BG=1;
    static final float[] RES={0.8f,1.2f,1.5f}, TRES={0.8f,1.2f,1.5f};
    static final float[] INITIAL={0.04f,-0.03f,0.08f,-1.1f,-2.6f,-1.5f};
    static float image(int c,int x,int y,int z) {
        if(!(x>=3&&x<=7&&y>=2&&y<=7&&z>=2&&z<=6)||((x+2*y+3*z+c)%6)==0)return 0.0f;
        return ((x*11+y*7+z*3+c*5)%17+1)/20.0f;
    }
    static byte topology(int x,int y,int z) {
        if(x>=7&&x<10&&y>=1&&y<6&&z>=3&&z<8)return 3;
        if(x>=1&&x<5&&y>=3&&y<8&&z>=1&&z<6)return 2;
        return 1;
    }
    static float[][] rotation(float a,float b,float c) {
        double n=Math.sqrt(a*a+b*b+c*c),an=a,bn=b,cn=c;
        if(n>1.0){double m=n%2.0;if(m>1.0){an=m*a/n-2*a/n;bn=m*b/n-2*b/n;cn=m*c/n-2*c/n;}else{an=a*m/n;bn=b*m/n;cn=c*m/n;}}
        double a2=an*an,b2=bn*bn,c2=cn*cn,w=Math.sqrt(1-a2-b2-c2);
        return new float[][]{{(float)(1-2*b2-2*c2),(float)(2*an*bn-2*w*cn),(float)(2*an*cn+2*w*bn)},
          {(float)(2*an*bn+2*w*cn),(float)(1-2*a2-2*c2),(float)(2*bn*cn-2*w*an)},
          {(float)(2*an*cn-2*w*bn),(float)(2*bn*cn+2*w*an),(float)(1-2*a2-2*b2)}};
    }
    static float[] map(int x,int y,int z,float[] t,float[][] r) {
        float[] out=new float[3];
        float[] q={(x-NX/2.0f)*TRES[0]-t[3],(y-NY/2.0f)*TRES[1]-t[4],(z-NZ/2.0f)*TRES[2]-t[5]};
        for(int axis=0;axis<3;axis++)out[axis]=(r[0][axis]*q[0]+r[1][axis]*q[1]+r[2][axis]*q[2])/RES[axis]+(axis==0?NX:axis==1?NY:NZ)/2.0f;
        return out;
    }
    static int floor(float f){return f==(int)f?(int)f:(f>0?(int)f:(int)f-1);}
    static int ceil(float f){return f==(int)f?(int)f:(f>0?(int)f+1:(int)f);}
    static int bits(float f){return Float.floatToIntBits(f);}
    static int[] bits(float[] a){int[] out=new int[a.length];for(int i=0;i<a.length;i++)out[i]=bits(a[i]);return out;}
    static void append(StringBuilder s,int[] a){s.append('[');for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(a[i]);}s.append(']');}
    public static void main(String[] args) {
        int x0=NX-1,xN=0,y0=NY-1,yN=0,z0=NZ-1,zN=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<NC;c++)if(image(c,x,y,z)>0){x0=Math.min(x0,x);xN=Math.max(xN,x);y0=Math.min(y0,y);yN=Math.max(yN,y);z0=Math.min(z0,z);zN=Math.max(zN,z);}
        float[][] r=rotation(INITIAL[0],INITIAL[1],INITIAL[2]);
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(topology(x,y,z)>BG){float[] v=map(x,y,z,INITIAL,r);x0=Math.min(x0,floor(v[0]));xN=Math.max(xN,ceil(v[0]));y0=Math.min(y0,floor(v[1]));yN=Math.max(yN,ceil(v[1]));z0=Math.min(z0,floor(v[2]));zN=Math.max(zN,ceil(v[2]));}
        int mx=xN-x0+1+2*BS,my=yN-y0+1+2*BS,mz=zN-z0+1+2*BS;
        float[] updated=INITIAL.clone();
        float dx=(mx/2.0f+x0-BS-NX/2.0f)*RES[0],dy=(my/2.0f+y0-BS-NY/2.0f)*RES[1],dz=(mz/2.0f+z0-BS-NZ/2.0f)*RES[2];
        updated[3]+=r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        updated[4]+=r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        updated[5]+=r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;
        int[] crop=new int[NC*mx*my*mz];
        for(int c=0;c<NC;c++)for(int x=x0;x<=xN;x++)for(int y=y0;y<=yN;y++)for(int z=z0;z<=zN;z++){
            int a=x-x0+BS,b=y-y0+BS,d=z-z0+BS;
            float value=(x<0||x>=NX||y<0||y>=NY||z<0||z>=NZ)?0:image(c,x,y,z);
            crop[((c*mx+a)*my+b)*mz+d]=bits(value);
        }
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08\",\"bounds_inclusive\":["+x0+","+xN+","+y0+","+yN+","+z0+","+zN+"],\"padding\":"+BS+",\"cropped_dims_xyz\":["+mx+","+my+","+mz+"],\"cropped_transform_bits\":");
        append(out,bits(updated));out.append(",\"cropped_image_bits\":");append(out,crop);out.append(",\"call_order\":[\"normalizeImages\",\"findCroppingBoundaries\",\"cropImages\"]}");
        System.out.println(out);
    }
}
