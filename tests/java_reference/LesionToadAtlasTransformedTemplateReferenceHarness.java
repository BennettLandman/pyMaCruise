/**
 * Pinned numerical excerpt of LesionToadAtlas.computeTransformedTemplate,
 * ParametricTransform.imageToTemplate (rigid), RotationMatrix.setParameters,
 * ImageFunctions.nearestNeighborInterpolation(byte[][][], byte, ...), and
 * Numerics.round from selected VandyCRUISE source.
 */
public final class LesionToadAtlasTransformedTemplateReferenceHarness {
    private static final int NX=7, NY=8, NZ=9, NTX=9, NTY=10, NTZ=11;
    private static final float[] P={0.02f,-0.03f,0.01f,0.3f,-0.15f,0.4f};
    private static float[][] rotation(float a, float b, float c) {
        double norm=Math.sqrt(a*a+b*b+c*c), an=a, bn=b, cn=c;
        if(norm>1.0){double m=norm%2;if(m>1.0){an=m*a/norm-2*a/norm;bn=m*b/norm-2*b/norm;cn=m*c/norm-2*c/norm;}else{an=a*m/norm;bn=b*m/norm;cn=c*m/norm;}}
        double a2=an*an,b2=bn*bn,c2=cn*cn,w2=1-a2-b2-c2;
        if(w2<0&&w2>-1e-10)w2=0; double w=Math.sqrt(w2);
        return new float[][]{{(float)(1-2*b2-2*c2),(float)(2*an*bn-2*w*cn),(float)(2*an*cn+2*w*bn)},
            {(float)(2*an*bn+2*w*cn),(float)(1-2*a2-2*c2),(float)(2*bn*cn-2*w*an)},
            {(float)(2*an*cn-2*w*bn),(float)(2*bn*cn+2*w*an),(float)(1-2*a2-2*b2)}};
    }
    private static int round(float n) { if(n==(int)n)return (int)n; return n>0?(int)(n+0.5f):(int)(n-0.5f); }
    private static byte sample(byte[][][] atlas, byte zero, float x, float y, float z) {
        if(x<0||x>NTX-1||y<0||y>NTY-1||z<0||z>NTZ-1)return zero;
        return atlas[round(x)][round(y)][round(z)];
    }
    private static byte[][][] map(byte[][][] atlas) {
        float[] imageRes={1.1f,1.3f,0.9f}, templateRes={0.9f,1.0f,1.1f};
        float x0i=NX/2.0f,y0i=NY/2.0f,z0i=NZ/2.0f;
        float x0t=NTX/2.0f,y0t=NTY/2.0f,z0t=NTZ/2.0f;
        float[][] rot=rotation(P[0],P[1],P[2]);
        byte[][][] out=new byte[NX][NY][NZ];
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            float Xi=(x*1.0f-x0i)*imageRes[0], Yi=(y*1.0f-y0i)*imageRes[1], Zi=(z*1.0f-z0i)*imageRes[2];
            float xt=((rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+P[3])/templateRes[0]+x0t)/1.0f;
            float yt=((rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+P[4])/templateRes[1]+y0t)/1.0f;
            float zt=((rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+P[5])/templateRes[2]+z0t)/1.0f;
            out[x][y][z]=sample(atlas,(byte)1,xt,yt,zt);
        }
        return out;
    }
    public static void main(String[] args) {
        byte[][][] atlas=new byte[NTX][NTY][NTZ];
        for(int x=0;x<NTX;x++)for(int y=0;y<NTY;y++)for(int z=0;z<NTZ;z++) atlas[x][y][z]=(byte)(1+((x*7+y*3+z*5+x*y+2*z*y)%17));
        // Sentinel distinguishes class metadata label[0] from template[0][0][0].
        atlas[0][0][0]=(byte)42;
        byte[][][] out=map(atlas);
        System.out.print("{\"source_method\":\"LesionToadAtlas.computeTransformedTemplate\",\"transform_mode\":\"rigid\",\"template_shape\":[9,10,11],\"image_shape\":[7,8,9],\"background_label\":1,\"output_labels\":[");
        boolean first=true; for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){if(!first)System.out.print(",");first=false;System.out.print(out[x][y][z]);}
        System.out.println("]}");
    }
}
