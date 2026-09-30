import edu.jhmi.rad.medic.utilities.PowerTable;

/** Selected mixed-category computeCentroids() branches in one class vector. */
public final class LesionMixedCentroidsReferenceHarness {
    private static final int NX=6,NY=5,NZ=4,NC=2,C=13;
    private static final float ZERO=1e-10f,SMOOTH=0.35f;
    private static final String SOURCE_SHA="77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0";
    private static final String POWER_SHA="2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae";
    private static final String[] TYPE={"object","mask","outmask","outlier","optimized_csf","optimized_gm","optimized_wm","lesion","blackhole","optimized_lesion","optimized_blackhole","rayleigh","object"};
    private static final String[] MODALITY={"T1_MPRAGE","FLAIR"};
    private static final float[][] INITIAL=new float[NC][C],PRIOR=new float[NC][C],TISSUE_POWER=new float[NC][C];
    private static final float[] LESION_C={0.85f,0.95f},BLACK_C={0.45f,0.52f};
    private static final float[] LESION_POWER={1.2f,0.8f},OPT_AUX_POWER={1.25f,0.75f};
    private static final int[][] LINKS=new int[NC][C];
    private static final int[][][] seg=new int[NX][NY][NZ];
    private static final float[][][][] mem=new float[NX][NY][NZ][C],image=new float[NX][NY][NZ][NC],field=new float[NX][NY][NZ][NC];
    private static final float[][][] outlier=new float[NX][NY][NZ],lesion=new float[NX][NY][NZ],black=new float[NX][NY][NZ];

    private static int bits(float x){return Float.floatToRawIntBits(x);}
    private static float[][][] update(String mode,boolean useField){
        float[][] cls=new float[NC][C];for(int c=0;c<NC;c++)cls[c]=INITIAL[c].clone();
        float[] lesionC=LESION_C.clone(),blackC=BLACK_C.clone();
        PowerTable[][] tables=new PowerTable[NC][C];
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++)tables[c][k]=new PowerTable(0.0f,4.0f,0.000001f,TISSUE_POWER[c][k]);
        PowerTable[] lesionTables=new PowerTable[NC],auxTables=new PowerTable[NC];
        for(int c=0;c<NC;c++){lesionTables[c]=new PowerTable(0.0f,4.0f,0.000001f,LESION_POWER[c]);auxTables[c]=new PowerTable(0.0f,4.0f,0.000001f,OPT_AUX_POWER[c]);}
        for(int c=0;c<NC;c++){
            float[] num=new float[C],den=new float[C];float numL=0,denL=0,numB=0,denB=0;
            for(int k=0;k<C;k++)for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++)if(seg[x][y][z]==k){
                float m=mem[x][y][z][k],sample=useField?field[x][y][z][c]*image[x][y][z][c]:image[x][y][z][c];
                float factor=1.0f;
                if(TYPE[k].equals("outlier"))m*=1.0f-outlier[x][y][z];
                if(TYPE[k].startsWith("optimized_")){
                    boolean tissue=TYPE[k].equals("optimized_csf")||TYPE[k].equals("optimized_gm")||TYPE[k].equals("optimized_wm");
                    if(tissue && (!TYPE[k].equals("optimized_csf")||MODALITY[c].equals("T1_MPRAGE"))){
                        float d=sample-INITIAL[c][k];d*=d;if(d>ZERO)factor=(float)(tables[c][k].lookup((double)d,(double)TISSUE_POWER[c][k])/d);
                    }else if(!tissue){
                        float winner=Math.max(m,Math.max(lesion[x][y][z],black[x][y][z]));float d=winner*winner*sample;
                        if(d>ZERO)factor=(float)(auxTables[c].lookup((double)d,(double)OPT_AUX_POWER[c])/d);
                        if(winner==m){num[k]+=factor*m*m*sample;den[k]+=factor*m*m;}
                        else if(winner==lesion[x][y][z]){float a=lesion[x][y][z];numL+=factor*a*a*sample;denL+=factor*a*a;}
                        else{float a=black[x][y][z];numB+=factor*a*a*sample;denB+=factor*a*a;}
                        continue;
                    }
                }else if(TYPE[k].equals("lesion")||TYPE[k].equals("blackhole")){
                    float lm=lesion[x][y][z],bm=black[x][y][z],winner=Math.max(m,Math.max(lm,bm));
                    if(winner!=m&&winner==lm){float d=sample-LESION_C[c];d*=d;if(d>ZERO)factor=(float)(lesionTables[c].lookup((double)d,(double)LESION_POWER[c])/d);}
                    if(winner==m){num[k]+=factor*m*m*sample;den[k]+=factor*m*m;}
                    else if(winner==lm){numL+=factor*lm*lm*sample;denL+=factor*lm*lm;}
                    else{numB+=factor*bm*bm*sample;denB+=factor*bm*bm;}
                    continue;
                }
                float square=m*m,weighted=TYPE[k].startsWith("optimized_")?factor*m*m:square;
                float term=weighted*sample;
                if(TYPE[k].equals("rayleigh"))term*=sample;
                num[k]+=term;den[k]+=weighted;
            }
            for(int k=0;k<C;k++){
                boolean ray=TYPE[k].equals("rayleigh")&&!mode.equals("prior");
                if(mode.equals("none")){cls[c][k]=den[k]>0?(ray?(float)Math.sqrt(num[k]/(2.0f*den[k])):num[k]/den[k]):0;}
                else if(mode.equals("linked")){float ns=num[k],ds=den[k];for(int l=0;l<C;l++)if(l!=k&&LINKS[c][k]==LINKS[c][l]){ns+=SMOOTH*num[l];ds+=SMOOTH*den[l];}cls[c][k]=ds>0?ns/ds:0;}
                else if(mode.equals("damped")){if(den[k]>0)cls[c][k]=(1-SMOOTH)*num[k]/den[k]+SMOOTH*cls[c][k];}
                else if(mode.equals("prior")){float pw=SMOOTH*1e4f,ns=num[k]+pw*PRIOR[c][k],ds=den[k]+pw;cls[c][k]=ds>0?ns/ds:0;}
                else{cls[c][k]=den[k]>0?num[k]/den[k]:0;}
            }
            if(denL>0)lesionC[c]=0.2f*(numL/denL)+0.8f*lesionC[c];else lesionC[c]=0;
            blackC[c]=denB>0?numB/denB:0;
        }
        return new float[][][]{cls,{lesionC},{blackC}};
    }
    private static String matrix(float[][] a){StringBuilder b=new StringBuilder("[");for(int c=0;c<a.length;c++){if(c>0)b.append(',');b.append('[');for(int k=0;k<a[c].length;k++){if(k>0)b.append(',');b.append(bits(a[c][k]));}b.append(']');}return b.append(']').toString();}
    private static String vector(float[] a){StringBuilder b=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)b.append(',');b.append(bits(a[i]));}return b.append(']').toString();}
    public static void main(String[] args){
        for(int c=0;c<NC;c++)for(int k=0;k<C;k++){INITIAL[c][k]=0.15f+0.17f*k+0.15f*c;PRIOR[c][k]=0.2f+0.13f*k+0.1f*c;TISSUE_POWER[c][k]=0.75f+0.11f*((k+2*c)%7);LINKS[c][k]=(k+2*c)%5;}
        for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){
            seg[x][y][z]=(x+2*y+z)%12;
            for(int k=0;k<C;k++)mem[x][y][z][k]=0.08f+((x*3+y*5+z*7+k*2)%9)*0.055f;
            for(int c=0;c<NC;c++){image[x][y][z][c]=0.25f+((x*11+y*7+z*3+c*5)%23)*0.12f;field[x][y][z][c]=0.75f+((x*2+y*3+z*4+c)%5)*0.11f;}
            outlier[x][y][z]=((x+3*y+z)%5)*0.13f;lesion[x][y][z]=0.04f+((x+2*y+z)%9)*0.07f;black[x][y][z]=0.03f+((2*x+y+3*z)%9)*0.075f;
            if((x+y+z)%6==0)black[x][y][z]=lesion[x][y][z];
        }
        int[][] p=new int[12][3];boolean[] set=new boolean[12];for(int x=0;x<NX;x++)for(int y=0;y<NY;y++)for(int z=0;z<NZ;z++){int k=seg[x][y][z];if(!set[k]){p[k]=new int[]{x,y,z};set[k]=true;}}
        int x=p[7][0],y=p[7][1],z=p[7][2];mem[x][y][z][7]=0.8f;lesion[x][y][z]=0.8f;black[x][y][z]=0.2f;
        x=p[8][0];y=p[8][1];z=p[8][2];mem[x][y][z][8]=0.3f;lesion[x][y][z]=0.8f;black[x][y][z]=0.8f;
        x=p[9][0];y=p[9][1];z=p[9][2];mem[x][y][z][9]=0.6f;lesion[x][y][z]=0.6f;black[x][y][z]=0.6f;
        x=p[10][0];y=p[10][1];z=p[10][2];mem[x][y][z][10]=0.3f;lesion[x][y][z]=0.2f;black[x][y][z]=0.9f;
        int seen=0;for(int a=0;a<NX;a++)for(int b=0;b<NY;b++)for(int d=0;d<NZ;d++)if(seg[a][b][d]==9&&seen++==1){mem[a][b][d][9]=1e-6f;lesion[a][b][d]=black[a][b][d]=0;}
        System.out.print("{\"source_sha256\":\""+SOURCE_SHA+"\",\"powertable_sha256\":\""+POWER_SHA+"\"");
        for(String mode:new String[]{"none","linked","damped","prior","fallback"})for(int f=0;f<2;f++){float[][][] o=update(mode,f==1);String suffix=f==1?"with_field":"no_field";System.out.print(",\""+mode+"_"+suffix+"\":{\"class\":"+matrix(o[0])+",\"lesion\":"+vector(o[1][0])+",\"blackhole\":"+vector(o[2][0])+"}");}
        System.out.print("}");
    }
}
