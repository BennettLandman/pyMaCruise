/*
 * Test-only source excerpt harness for selected toads-cruise methods.
 * The body of findCroppingBoundaries and computeImageCoordinates below follows
 * the selected source lines directly; dependencies and provenance are recorded
 * in tests/data/synthetic/toads_crop_geometry/README.md.
 */
import java.util.Locale;

public final class ToadsCropGeometryReferenceHarness {
    static final int NX = 7, NY = 8, NZ = 6;
    static final int TX = 7, TY = 8, TZ = 6;
    static final int BORDER = 2, BG = 1;
    static final float[] IMAGE_RES = {1.0f, 1.25f, 1.5f};
    static final float[] TEMPLATE_RES = {1.0f, 1.0f, 1.0f};
    static final float[] ORIENT = {-1.0f, -1.0f, 1.0f};
    static final float[] ROT = {0.12f, -0.08f, 0.05f};
    static final float[] TRANS = {0.4f, -0.25f, 0.3f};

    static final float[][][] image = new float[NX][NY][NZ];
    static final byte[][][] topology = new byte[TX][TY][TZ];

    static int x0, xN, y0, yN, z0, zN;
    static int mx, my, mz;
    static float[] transform = new float[6];
    static float[] low = {0.0f}, high = {1.0f};
    static boolean normalized = true;
    static float bgRatio = 0.2f;

    static float[][] rotation(float a, float b, float c) {
        // Transcribed from RotationMatrix.setParameters(float,float,float).
        double n = Math.sqrt(a*a + b*b + c*c), an=a, bn=b, cn=c;
        if (n > 1.0) {
            double m = n % 2.0;
            if (m > 1.0) {
                an = m*a/n - 2*a/n; bn = m*b/n - 2*b/n; cn = m*c/n - 2*c/n;
            } else { an = a*m/n; bn = b*m/n; cn = c*m/n; }
        }
        double a2=an*an, b2=bn*bn, c2=cn*cn;
        double w=Math.sqrt(1.0-a2-b2-c2);
        float[][] r = new float[3][3];
        r[0][0]=(float)(1-2*b2-2*c2); r[1][0]=(float)(2*an*bn+2*w*cn); r[2][0]=(float)(2*an*cn-2*w*bn);
        r[0][1]=(float)(2*an*bn-2*w*cn); r[1][1]=(float)(1-2*a2-2*c2); r[2][1]=(float)(2*bn*cn+2*w*an);
        r[0][2]=(float)(2*an*cn+2*w*bn); r[1][2]=(float)(2*bn*cn-2*w*an); r[2][2]=(float)(1-2*a2-2*b2);
        return r;
    }

    static float[] computeImageCoordinates(int x,int y,int z,float[] trans,float[][] rot) {
        // Selected-source method body: MultipleToadPreprocess.java:653-663.
        float xCt=TX/2.0f,yCt=TY/2.0f,zCt=TZ/2.0f;
        float xCi=NX/2.0f,yCi=NY/2.0f,zCi=NZ/2.0f;
        float rtx=TEMPLATE_RES[0],rty=TEMPLATE_RES[1],rtz=TEMPLATE_RES[2];
        float rix=IMAGE_RES[0],riy=IMAGE_RES[1],riz=IMAGE_RES[2];
        float orx=ORIENT[0],ory=ORIENT[1],orz=ORIENT[2];
        float[] XI = new float[3];
        XI[0] = (rot[0][0]*((x-xCt)*rtx-trans[3])+rot[1][0]*((y-yCt)*rty-trans[4])+rot[2][0]*((z-zCt)*rtz-trans[5]))*orx/rix + xCi;
        XI[1] = (rot[0][1]*((x-xCt)*rtx-trans[3])+rot[1][1]*((y-yCt)*rty-trans[4])+rot[2][1]*((z-zCt)*rtz-trans[5]))*ory/riy + yCi;
        XI[2] = (rot[0][2]*((x-xCt)*rtx-trans[3])+rot[1][2]*((y-yCt)*rty-trans[4])+rot[2][2]*((z-zCt)*rtz-trans[5]))*orz/riz + zCi;
        return XI;
    }

    static float[] templateCoordinates(int x, int y, int z, float[][] r) {
        // Transcribed from MultipleToadPreprocess.computeTemplateCoordinates.
        float xci=mx/2.0f, yci=my/2.0f, zci=mz/2.0f;
        return new float[] {
            (r[0][0]*(x-xci)*IMAGE_RES[0]/ORIENT[0]+r[0][1]*(y-yci)*IMAGE_RES[1]/ORIENT[1]+r[0][2]*(z-zci)*IMAGE_RES[2]/ORIENT[2]+transform[3])/TEMPLATE_RES[0]+TX/2.0f,
            (r[1][0]*(x-xci)*IMAGE_RES[0]/ORIENT[0]+r[1][1]*(y-yci)*IMAGE_RES[1]/ORIENT[1]+r[1][2]*(z-zci)*IMAGE_RES[2]/ORIENT[2]+transform[4])/TEMPLATE_RES[1]+TY/2.0f,
            (r[2][0]*(x-xci)*IMAGE_RES[0]/ORIENT[0]+r[2][1]*(y-yci)*IMAGE_RES[1]/ORIENT[1]+r[2][2]*(z-zci)*IMAGE_RES[2]/ORIENT[2]+transform[5])/TEMPLATE_RES[2]+TZ/2.0f
        };
    }

    static int floor(float f) { return f==(int)f ? (int)f : (f>0 ? (int)f : (int)f-1); }
    static int ceil(float f) { return f==(int)f ? (int)f : (f>0 ? (int)f+1 : (int)f); }
    static int round(float f) { return f==(int)f ? (int)f : (f>0 ? (int)(f+0.5f) : (int)(f-0.5f)); }
    static byte sampleTopology(float x, float y, float z) {
        if (x<0 || x>TX-1 || y<0 || y>TY-1 || z<0 || z>TZ-1) return (byte)BG;
        return topology[round(x)][round(y)][round(z)];
    }

    static void findCroppingBoundaries() {
        // Selected-source loop body: MultipleToadPreprocess.java:302-344.
        x0=NX-1; xN=0; y0=NY-1; yN=0; z0=NZ-1; zN=0;
        // Source scans channels inside each voxel; this fixture has nc=1.
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            for (int c=0;c<1;c++) {
                if (((normalized) && (image[x][y][z] > bgRatio))
                    || (image[x][y][z] > bgRatio*(high[c]-low[c]))) {
                    if (x<x0) x0=x; if (x>xN) xN=x;
                    if (y<y0) y0=y; if (y>yN) yN=y;
                    if (z<z0) z0=z; if (z>zN) zN=z;
                }
            }
        }
        float[][] r=rotation(ROT[0],ROT[1],ROT[2]);
        for (int x=0;x<TX;x++) for (int y=0;y<TY;y++) for (int z=0;z<TZ;z++) if (topology[x][y][z] > BG) {
            float[] v=computeImageCoordinates(x,y,z,transform,r);
            if (v[0]<x0) x0=floor(v[0]); if (v[0]>xN) xN=ceil(v[0]);
            if (v[1]<y0) y0=floor(v[1]); if (v[1]>yN) yN=ceil(v[1]);
            if (v[2]<z0) z0=floor(v[2]); if (v[2]>zN) zN=ceil(v[2]);
        }
    }

    static void crop() {
        mx=xN-x0+1+2*BORDER; my=yN-y0+1+2*BORDER; mz=zN-z0+1+2*BORDER;
        float[][] r=rotation(ROT[0],ROT[1],ROT[2]);
        float dx=(mx/2.0f+x0-BORDER-NX/2.0f)*IMAGE_RES[0]/ORIENT[0];
        float dy=(my/2.0f+y0-BORDER-NY/2.0f)*IMAGE_RES[1]/ORIENT[1];
        float dz=(mz/2.0f+z0-BORDER-NZ/2.0f)*IMAGE_RES[2]/ORIENT[2];
        // Java compound assignment rounds the float result after RHS double operations.
        transform[3] += r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        transform[4] += r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        transform[5] += r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;
        float[][][] smaller=new float[mx][my][mz];
        for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
            int a=x-x0+BORDER,b=y-y0+BORDER,c=z-z0+BORDER;
            smaller[a][b][c]=(x<0||x>=NX||y<0||y>=NY||z<0||z>=NZ)?0:image[x][y][z];
        }
        croppedImage=smaller;
    }
    static float[][][] croppedImage;

    static byte[][][] transformTopology() {
        byte[][][] out=new byte[mx][my][mz];
        for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) out[x][y][z]=(byte)BG;
        float[][] r=rotation(ROT[0],ROT[1],ROT[2]);
        for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
            float[] p=templateCoordinates(x,y,z,r);
            out[x-x0+BORDER][y-y0+BORDER][z-z0+BORDER]=sampleTopology(p[0],p[1],p[2]);
        }
        return out;
    }

    static float[][][] uncrop(float[][][] cropped) {
        // Inverse transform update and clipped padded copy from uncropImages.
        float[][] r=rotation(ROT[0],ROT[1],ROT[2]);
        float dx=(mx/2.0f+x0-BORDER-NX/2.0f)*IMAGE_RES[0]/ORIENT[0];
        float dy=(my/2.0f+y0-BORDER-NY/2.0f)*IMAGE_RES[1]/ORIENT[1];
        float dz=(mz/2.0f+z0-BORDER-NZ/2.0f)*IMAGE_RES[2]/ORIENT[2];
        transform[3] -= r[0][0]*dx+r[0][1]*dy+r[0][2]*dz;
        transform[4] -= r[1][0]*dx+r[1][1]*dy+r[1][2]*dz;
        transform[5] -= r[2][0]*dx+r[2][1]*dy+r[2][2]*dz;
        float[][][] larger=new float[NX][NY][NZ];
        for(int x=Math.max(x0-BORDER,0);x<=Math.min(xN+BORDER,NX-1);x++)
            for(int y=Math.max(y0-BORDER,0);y<=Math.min(yN+BORDER,NY-1);y++)
                for(int z=Math.max(z0-BORDER,0);z<=Math.min(zN+BORDER,NZ-1);z++)
                    larger[x][y][z]=cropped[x-x0+BORDER][y-y0+BORDER][z-z0+BORDER];
        return larger;
    }

    static String floats(float[] a) { StringBuilder s=new StringBuilder("["); for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(String.format(Locale.ROOT,"%.9g",a[i]));} return s.append(']').toString(); }
    static void array3(StringBuilder s, float[][][] a) { s.append('['); for(int x=0;x<a.length;x++){if(x>0)s.append(',');s.append('[');for(int y=0;y<a[x].length;y++){if(y>0)s.append(',');s.append('[');for(int z=0;z<a[x][y].length;z++){if(z>0)s.append(',');s.append(String.format(Locale.ROOT,"%.9g",a[x][y][z]));}s.append(']');}s.append(']');}s.append(']'); }
    static void array3(StringBuilder s, byte[][][] a) { s.append('['); for(int x=0;x<a.length;x++){if(x>0)s.append(',');s.append('[');for(int y=0;y<a[x].length;y++){if(y>0)s.append(',');s.append('[');for(int z=0;z<a[x][y].length;z++){if(z>0)s.append(',');s.append(a[x][y][z]);}s.append(']');}s.append(']');}s.append(']'); }

    public static void main(String[] args) {
        for (int i=0;i<3;i++) { transform[i]=ROT[i]; transform[i+3]=TRANS[i]; }
        for (int x=0;x<TX;x++) for(int y=0;y<TY;y++) for(int z=0;z<TZ;z++) topology[x][y][z]=(byte)BG;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            if ((x>=1&&x<=5&&y>=2&&y<=6&&z>=1&&z<=4) && ((x+2*y+3*z)%5!=0)) image[x][y][z]=(float)((x*17+y*7+z*3)%13)/10.0f;
        }
        // Distinct labels exercise categorical nearest-neighbor transformation.
        for (int x=1;x<=3;x++) for(int y=2;y<=4;y++) for(int z=1;z<=3;z++) if(x==1||x==3||y==2||y==4||z==1||z==3) topology[x][y][z]=(byte)2;
        for (int x=3;x<=6;x++) for(int y=3;y<=6;y++) for(int z=2;z<=5;z++) topology[x][y][z]=(byte)3;
        findCroppingBoundaries();
        int[] bounds={x0,xN,y0,yN,z0,zN};
        crop();
        float[] croppedTransform=transform.clone();
        byte[][][] transformedTopology=transformTopology();
        float[][][] uncroppedImage=uncrop(croppedImage);
        float[] uncroppedTransform=transform.clone();
        StringBuilder s=new StringBuilder();
        s.append("{\n  \"dims_image\":[").append(NX).append(',').append(NY).append(',').append(NZ).append("],\n");
        s.append("  \"dims_template\":[").append(TX).append(',').append(TY).append(',').append(TZ).append("],\n");
        s.append("  \"image_resolution\":").append(floats(IMAGE_RES)).append(",\n  \"template_resolution\":").append(floats(TEMPLATE_RES)).append(",\n  \"image_orientation_signs\":").append(floats(ORIENT)).append(",\n");
        s.append("  \"initial_transform\":").append(floats(new float[]{ROT[0],ROT[1],ROT[2],TRANS[0],TRANS[1],TRANS[2]})).append(",\n");
        s.append("  \"bounds_inclusive\":[").append(bounds[0]).append(',').append(bounds[1]).append(',').append(bounds[2]).append(',').append(bounds[3]).append(',').append(bounds[4]).append(',').append(bounds[5]).append("],\n");
        s.append("  \"padding\":").append(BORDER).append(",\n  \"cropped_dims\":[").append(mx).append(',').append(my).append(',').append(mz).append("],\n");
        s.append("  \"cropped_transform\":").append(floats(croppedTransform)).append(",\n  \"cropped_image\":"); array3(s,croppedImage);
        s.append(",\n  \"transformed_topology\":"); array3(s,transformedTopology);
        s.append(",\n  \"uncropped_transform\":").append(floats(uncroppedTransform)).append(",\n  \"uncropped_image\":"); array3(s,uncroppedImage);
        s.append("\n}\n");
        System.out.print(s);
    }
}
