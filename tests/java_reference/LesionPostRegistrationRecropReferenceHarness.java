/*
 * Test-only excerpt harness for LesionToadPreprocess crop lifecycle methods.
 * It runs the selected source sequence twice: initial crop, uncrop after atlas
 * registration, then topology-aware recrop. MIPAV/JIST dependencies are not
 * required. The copied method scope and source hashes are recorded beside the
 * generated reference JSON.
 */
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.MessageDigest;
import java.util.Locale;

public final class LesionPostRegistrationRecropReferenceHarness {
    static final int NX=11, NY=10, NZ=9, NC=2, BORDER=2, BG=1;
    static final float[] IRES={0.8f,1.2f,1.5f}, TRES={0.8f,1.2f,1.5f};
    static final float[] ORIENT={1.0f,1.0f,1.0f};
    static final float[] ROT={0.04f,-0.03f,0.08f};
    static final float[] INITIAL={-1.1f,-2.6f,-1.5f};
    static final float[][][][] images=new float[NC][NX][NY][NZ];
    static final byte[][][] topology=new byte[NX][NY][NZ];
    static float[] transform=new float[6];
    static int x0,xN,y0,yN,z0,zN,mx,my,mz;

    static float[][] rotation(float a,float b,float c) {
        double n=Math.sqrt(a*a+b*b+c*c), an=a,bn=b,cn=c;
        if(n>1.0){double m=n%2.0;if(m>1.0){an=m*a/n-2*a/n;bn=m*b/n-2*b/n;cn=m*c/n-2*c/n;}else{an=a*m/n;bn=b*m/n;cn=c*m/n;}}
        double a2=an*an,b2=bn*bn,c2=cn*cn,w=Math.sqrt(1-a2-b2-c2);
        return new float[][]{
            {(float)(1-2*b2-2*c2),(float)(2*an*bn-2*w*cn),(float)(2*an*cn+2*w*bn)},
            {(float)(2*an*bn+2*w*cn),(float)(1-2*a2-2*c2),(float)(2*bn*cn-2*w*an)},
            {(float)(2*an*cn-2*w*bn),(float)(2*bn*cn+2*w*an),(float)(1-2*a2-2*b2)}
        };
    }

    static float[] imageCoordinates(int x,int y,int z,float[][] r){
        float[] p=new float[3];float[] d={(x-NX/2.0f)*TRES[0]-transform[3],(y-NY/2.0f)*TRES[1]-transform[4],(z-NZ/2.0f)*TRES[2]-transform[5]};
        for(int a=0;a<3;a++)p[a]=(r[0][a]*d[0]+r[1][a]*d[1]+r[2][a]*d[2])*ORIENT[a]/IRES[a]+(a==0?NX:a==1?NY:NZ)/2.0f;
        return p;
    }
    static int floor(float v){return (int)Math.floor(v);}
    static int ceil(float v){return (int)Math.ceil(v);}

    static void findCroppingBoundaries(){
        x0=NX-1;xN=0;y0=NY-1;yN=0;z0=NZ-1;zN=0;
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)for(int c=0;c<NC;c++)if(images[c][x][y][z]>0.0f){
            x0=Math.min(x0,x);xN=Math.max(xN,x);y0=Math.min(y0,y);yN=Math.max(yN,y);z0=Math.min(z0,z);zN=Math.max(zN,z);
        }
        float[][] r=rotation(transform[0],transform[1],transform[2]);
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(topology[x][y][z]>BG){
            float[] p=imageCoordinates(x,y,z,r);
            if(p[0]<x0)x0=floor(p[0]);if(p[0]>xN)xN=ceil(p[0]);
            if(p[1]<y0)y0=floor(p[1]);if(p[1]>yN)yN=ceil(p[1]);
            if(p[2]<z0)z0=floor(p[2]);if(p[2]>zN)zN=ceil(p[2]);
        }
    }

    static float[][][][] crop(){
        mx=xN-x0+1+2*BORDER;my=yN-y0+1+2*BORDER;mz=zN-z0+1+2*BORDER;
        float[][] r=rotation(transform[0],transform[1],transform[2]);
        float dx=(mx/2.0f+x0-BORDER-NX/2.0f)*IRES[0]/ORIENT[0];
        float dy=(my/2.0f+y0-BORDER-NY/2.0f)*IRES[1]/ORIENT[1];
        float dz=(mz/2.0f+z0-BORDER-NZ/2.0f)*IRES[2]/ORIENT[2];
        transform[3]+=r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        transform[4]+=r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        transform[5]+=r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;
        float[][][][] out=new float[NC][mx][my][mz];
        for(int c=0;c<NC;c++)for(int x=x0;x<=xN;x++)for(int y=y0;y<=yN;y++)for(int z=z0;z<=zN;z++){
            int xx=x-x0+BORDER,yy=y-y0+BORDER,zz=z-z0+BORDER;
            if(x>=0&&x<NX&&y>=0&&y<NY&&z>=0&&z<NZ)out[c][xx][yy][zz]=images[c][x][y][z];
        }
        return out;
    }

    static float[][][][] uncrop(float[][][][] source){
        float[][] r=rotation(transform[0],transform[1],transform[2]);
        float dx=(mx/2.0f+x0-BORDER-NX/2.0f)*IRES[0]/ORIENT[0];
        float dy=(my/2.0f+y0-BORDER-NY/2.0f)*IRES[1]/ORIENT[1];
        float dz=(mz/2.0f+z0-BORDER-NZ/2.0f)*IRES[2]/ORIENT[2];
        transform[3]-=r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        transform[4]-=r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        transform[5]-=r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;
        float[][][][] out=new float[NC][NX][NY][NZ];
        for(int c=0;c<NC;c++)for(int x=Math.max(x0-BORDER,0);x<=Math.min(xN+BORDER,NX-1);x++)for(int y=Math.max(y0-BORDER,0);y<=Math.min(yN+BORDER,NY-1);y++)for(int z=Math.max(z0-BORDER,0);z<=Math.min(zN+BORDER,NZ-1);z++)out[c][x][y][z]=source[c][x-x0+BORDER][y-y0+BORDER][z-z0+BORDER];
        return out;
    }

    static String bits(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(Float.floatToIntBits(a[i]));}return s.append(']').toString();}
    static String digest(float[][][][] a)throws Exception{
        MessageDigest md=MessageDigest.getInstance("SHA-256");ByteBuffer b=ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        for(int c=0;c<NC;c++)for(int x=0;x<a[c].length;x++)for(int y=0;y<a[c][x].length;y++)for(int z=0;z<a[c][x][y].length;z++){b.clear();b.putFloat(a[c][x][y][z]);md.update(b.array());}
        StringBuilder s=new StringBuilder();for(byte v:md.digest())s.append(String.format(Locale.ROOT,"%02x",v&255));return s.toString();
    }
    static String dims(){return "["+mx+","+my+","+mz+"]";}
    static String bounds(){return "["+x0+","+xN+","+y0+","+yN+","+z0+","+zN+"]";}

    public static void main(String[] args)throws Exception{
        for(int i=0;i<3;i++){transform[i]=ROT[i];transform[i+3]=INITIAL[i];}
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            for(int c=0;c<NC;c++)if(x>=3&&x<=7&&y>=2&&y<=7&&z>=2&&z<=6&&(x+2*y+3*z+c)%6!=0)images[c][x][y][z]=(float)(((x*11+y*7+z*3+c*5)%17+1)/20.0f);
            topology[x][y][z]=1;
        }
        for(int x=1;x<5;x++)for(int y=3;y<8;y++)for(int z=1;z<6;z++)topology[x][y][z]=5;
        for(int x=7;x<10;x++)for(int y=1;y<6;y++)for(int z=3;z<8;z++)topology[x][y][z]=6;
        findCroppingBoundaries();String firstBounds=bounds();float[][][][] firstCrop=crop();String firstDims=dims();float[] firstTransform=transform.clone();String firstDigest=digest(firstCrop);
        float[][][][] restored=uncrop(firstCrop);
        for(int c=0;c<NC;c++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)images[c][x][y][z]=restored[c][x][y][z];
        float[] restoredTransform=transform.clone();
        findCroppingBoundaries();String secondBounds=bounds();float[][][][] secondCrop=crop();
        System.out.println("{\"first_bounds\":"+firstBounds+",\"first_dimensions\":"+firstDims+",\"first_transform_bits\":"+bits(firstTransform)+",\"first_crop_sha256\":\""+firstDigest+"\",\"restored_transform_bits\":"+bits(restoredTransform)+",\"second_bounds\":"+secondBounds+",\"second_dimensions\":"+dims()+",\"second_transform_bits\":"+bits(transform)+",\"second_crop_sha256\":\""+digest(secondCrop)+"\"}");
    }
}
