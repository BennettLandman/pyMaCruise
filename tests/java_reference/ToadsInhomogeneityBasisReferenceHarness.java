import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

/** Selected InhomogeneityCorrection2.compute3DChebyshevCoefficients formula excerpt. */
public final class ToadsInhomogeneityBasisReferenceHarness {
    private static final String SOURCE_SHA256="6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8";
    private static int bits(float value){return Float.floatToIntBits(value);}
    private static String basis(int x,int y,int z,int nx,int ny,int nz,int degree){
        float[] pol=new float[20];float rx=x/(float)(nx-1),ry=y/(float)(ny-1),rz=z/(float)(nz-1);
        pol[0]=1.0f;pol[1]=rx;pol[2]=ry;pol[3]=rz;
        if(degree>=2){pol[4]=2.0f*rx*rx-1.0f;pol[5]=2.0f*ry*ry-1.0f;pol[6]=2.0f*rz*rz-1.0f;pol[7]=pol[1]*pol[2];pol[8]=pol[2]*pol[3];pol[9]=pol[3]*pol[1];}
        if(degree>=3){pol[10]=4.0f*rx*rx*rx-3.0f*rx;pol[11]=4.0f*ry*ry*ry-3.0f*ry;pol[12]=4.0f*rz*rz*rz-3.0f*rz;pol[13]=pol[4]*pol[2];pol[14]=pol[1]*pol[5];pol[15]=pol[5]*pol[3];pol[16]=pol[2]*pol[6];pol[17]=pol[6]*pol[1];pol[18]=pol[3]*pol[4];pol[19]=pol[1]*pol[2]*pol[3];}
        int count=degree==1?4:degree==2?10:20;StringBuilder out=new StringBuilder("[");for(int i=0;i<count;i++){if(i>0)out.append(',');out.append(bits(pol[i]));}return out.append(']').toString();
    }
    public static void main(String[] args)throws IOException{
        byte[] bytes;try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(Path.of(args[0])));ByteArrayOutputStream out=new ByteArrayOutputStream()){in.transferTo(out);bytes=out.toByteArray();}
        ByteBuffer h=ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);if(h.getInt(0)!=348)throw new IOException("expected NIfTI-1");int nx=h.getShort(42),ny=h.getShort(44),nz=h.getShort(46);
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA256).append("\",\"dimensions\":[").append(nx).append(',').append(ny).append(',').append(nz).append("],\"degrees\":[");
        for(int degree=1;degree<=3;degree++){if(degree>1)out.append(',');out.append("{\"degree\":").append(degree).append(",\"values\":[");boolean first=true;for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++){if(!first)out.append(',');first=false;out.append(basis(x,y,z,nx,ny,nz,degree));}out.append("]}");}
        System.out.println(out.append("]}").toString());
    }
}
