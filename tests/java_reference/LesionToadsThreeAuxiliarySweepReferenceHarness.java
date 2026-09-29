import java.util.Locale;

/** Selected VandyCRUISE OBJECT+LESION+BLACKHOLE+OUTMASK formula harness. */
public final class LesionToadsThreeAuxiliarySweepReferenceHarness {
    private static int bits(float value) { return Float.floatToIntBits(value); }
    private static float factor(float shape) { float doubled=2.0f*shape;return (1.0f<doubled)?1.0f:doubled; }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        final int nx=5,ny=4,nz=6,nc=5;
        final float smoothing=0.2f,coefficient=0.65f,zero=1.0e-10f,infinity=1.0e10f;
        final float[] smoothPrior={0.25f,0.5f,0.75f,0.375f,0.625f};
        final float[] atlasPrior={0.5f,0.75f,1.25f,0.875f,0.625f};
        final float[] lesionSimilarity={0.25f,0.0f,0.5f,0.375f,0.625f};
        final float[] blackHoleSimilarity={0.125f,0.5f,0.0f,0.25f,0.375f};
        final float[][] similarity={{0.0f,0.2f,0.7f,0.3f,0.1f},{0.4f,0.0f,0.3f,0.6f,0.2f},{0.6f,0.1f,0.0f,0.2f,0.5f},{0.2f,0.5f,0.4f,0.0f,0.3f},{0.3f,0.4f,0.1f,0.7f,0.0f}};
        final String[] type=args.length>0&&args[0].equals("lesion-last")
            ?new String[]{"OBJECT","BLACKHOLE","OUTMASK","LESION","OBJECT"}
            :new String[]{"OBJECT","LESION","BLACKHOLE","OUTMASK","OBJECT"};
        int lesionClass=-1,lesionTypeClass=-1,blackHoleClass=-1,outmaskClass=-1;
        for(int k=0;k<nc;k++) {
            if(type[k].equals("LESION")||type[k].equals("BLACKHOLE"))lesionClass=k;
            if(type[k].equals("LESION"))lesionTypeClass=k;
            if(type[k].equals("BLACKHOLE"))blackHoleClass=k;
            if(type[k].equals("OUTMASK"))outmaskClass=k;
        }
        final int[] relation={0,1,1,1,0};
        float[][][][] mems=new float[nx][ny][nz][nc],data=new float[nx][ny][nz][nc],shape=new float[nx][ny][nz][nc];
        float[][][] lesionData=new float[nx][ny][nz],blackHoleData=new float[nx][ny][nz],outmaskData=new float[nx][ny][nz];
        float[][][] lesion=new float[nx][ny][nz],blackHole=new float[nx][ny][nz],outlier=new float[nx][ny][nz];
        byte[][][] vent=new byte[nx][ny][nz],gm=new byte[nx][ny][nz],brainstem=new byte[nx][ny][nz];
        for(int x=0;x<nx;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++) {
            lesionData[x][y][z]=0.625f+0.03125f*x+0.015625f*y+0.0078125f*z;
            blackHoleData[x][y][z]=0.5f+0.03125f*x+0.015625f*y+0.0078125f*z;
            outmaskData[x][y][z]=0.375f+0.03125f*x+0.015625f*y+0.0078125f*z;
            lesion[x][y][z]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z;
            blackHole[x][y][z]=0.125f+0.03125f*x+0.015625f*y+0.0078125f*z;
            outlier[x][y][z]=0.125f*x+0.0625f*y+0.03125f*z;
            vent[x][y][z]=(byte)(x==1?0:x==2?3:-1);gm[x][y][z]=(byte)(y==1?0:4);brainstem[x][y][z]=(byte)(z==2?0:2);
            for(int k=0;k<nc;k++) {
                mems[x][y][z][k]=0.5f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                data[x][y][z][k]=0.75f+0.125f*x+0.0625f*y+0.03125f*z+0.0625f*k;
                shape[x][y][z][k]=0.25f+0.0625f*x+0.03125f*y+0.015625f*z+0.0625f*k;
            }
        }
        lesion[1][1][1]=0.0f;outlier[1][1][1]=1.0f;outlier[3][2][4]=0.0f;
        shape[2][2][2][blackHoleClass]=0.0f;shape[1][2][3][lesionTypeClass]=0.0f;
        int[][] delta={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"shape\":[5,4,6],\"classes\":5,\"class_order\":\"").append(type[lesionClass].equals("BLACKHOLE")?"blackhole-last":"lesion-last").append("\",\"voxels\":[");
        boolean first=true;
        for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
            if(!first)out.append(',');first=false;float den=0.0f;
            for(int k=0;k<nc;k++) {
                if(shape[x][y][z][k]>0.0f) {
                    float other=0.0f,lesionSum=0.0f,blackHoleSum=0.0f,currentSum=0.0f;
                    for(int[] f:delta) {
                        int a=x+f[0],b=y+f[1],c=z+f[2];
                        for(int m=0;m<nc;m++)if(m!=k){float v=mems[a][b][c][m];other+=v*v;}
                        float lv=lesion[a][b][c],bv=blackHole[a][b][c];lesionSum+=lv*lv;blackHoleSum+=bv*bv;
                    }
                    if(type[k].equals("LESION")||type[k].equals("BLACKHOLE")) {
                        int[][] order={{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
                        for(int[] f:order){float v=mems[x+f[0]][y+f[1]][z+f[2]][k];currentSum+=v*v;}
                    }
                    float num=data[x][y][z][k],numL=(type[k].equals("LESION")||type[k].equals("BLACKHOLE"))?lesionData[x][y][z]:0.0f;
                    float numB=type[k].equals("BLACKHOLE")?blackHoleData[x][y][z]:0.0f,numO=type[k].equals("OUTMASK")?outmaskData[x][y][z]:0.0f;
                    float all=other+lesionSum;all+=blackHoleSum;num+=smoothPrior[k]*smoothing*all;
                    if(type[k].equals("LESION")||type[k].equals("BLACKHOLE")){float ls=other+currentSum;ls+=blackHoleSum;numL+=smoothPrior[k]*smoothing*ls;}
                    if(type[k].equals("BLACKHOLE")){float bs=other+currentSum;bs+=lesionSum;numB+=smoothPrior[k]*smoothing*bs;}
                    if(type[k].equals("OUTMASK"))numO+=smoothPrior[k]*smoothing*(other+lesionSum);
                    float invprior=0.0f,invpriorL=0.0f,invpriorB=0.0f;
                    for(int m=0;m<nc;m++) {
                        if(m!=k){float q=atlasPrior[m]*similarity[k][m];q*=shape[x][y][z][m];q*=shape[x][y][z][m];invprior+=q;}
                        if(type[k].equals("LESION")||type[k].equals("BLACKHOLE")){float q=atlasPrior[m]*lesionSimilarity[m];q*=shape[x][y][z][m];q*=shape[x][y][z][m];invpriorL+=q;}
                        if(type[k].equals("BLACKHOLE")){float q=atlasPrior[m]*blackHoleSimilarity[m];q*=shape[x][y][z][m];q*=shape[x][y][z][m];invpriorB+=q;}
                    }
                    float selectedShape=shape[x][y][z][lesionClass];
                    float q=atlasPrior[lesionClass]*lesionSimilarity[k];q*=selectedShape;q*=selectedShape;invprior+=q;
                    if(type[lesionClass].equals("BLACKHOLE")){q=atlasPrior[lesionClass]*blackHoleSimilarity[k];q*=selectedShape;q*=selectedShape;invprior+=q;}
                    num+=coefficient*invprior;
                    if(type[k].equals("LESION")||type[k].equals("BLACKHOLE"))numL+=coefficient*invpriorL;
                    if(type[k].equals("BLACKHOLE"))numB+=coefficient*invpriorB;
                    if(type[k].equals("OUTMASK"))numO+=coefficient*invprior;
                    if(num>zero)num=1.0f/num;else num=infinity;
                    if(type[k].equals("LESION")||type[k].equals("BLACKHOLE")){if(numL>zero)numL=1.0f/numL;else numL=infinity;}
                    if(type[k].equals("BLACKHOLE")){if(numB>zero)numB=1.0f/numB;else numB=infinity;}
                    if(type[k].equals("OUTMASK")){if(numO>zero)numO=1.0f/numO;else numO=infinity;}
                    if(type[k].equals("LESION")) {
                        mems[x][y][z][k]=relation[k]==1?factor(shape[x][y][z][k])*num:num;
                        lesion[x][y][z]=relation[k]==1?factor(shape[x][y][z][k])*numL:numL;
                        float fac=relation[k]==1?factor(shape[x][y][z][k]):1.0f;
                        if(vent[x][y][z]!=-1){float n=vent[x][y][z]==0?2.5f:vent[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/5.0f;fac*=n/5.0f;}
                        else if(gm[x][y][z]!=-1){float n=gm[x][y][z]==0?2.5f:gm[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/8.0f;fac*=n/8.0f;}
                        if(brainstem[x][y][z]!=-1){float n=brainstem[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/4.0f;fac*=n/4.0f;}
                    } else if(type[k].equals("BLACKHOLE")) {
                        mems[x][y][z][k]=relation[k]==1?factor(shape[x][y][z][k])*num:num;
                        lesion[x][y][z]=relation[k]==1?factor(shape[x][y][z][k])*numL:numL;
                        blackHole[x][y][z]=relation[k]==1?factor(shape[x][y][z][k])*numB:numB;
                        float fac=relation[k]==1?factor(shape[x][y][z][k]):1.0f;
                        if(vent[x][y][z]!=-1){float n=vent[x][y][z]==0?2.5f:vent[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/5.0f;fac*=n/5.0f;}
                        else if(gm[x][y][z]!=-1){float n=gm[x][y][z]==0?2.5f:gm[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/8.0f;blackHole[x][y][z]=n*blackHole[x][y][z]/8.0f;fac*=n/8.0f;}
                        if(brainstem[x][y][z]!=-1){float n=brainstem[x][y][z];lesion[x][y][z]=n*lesion[x][y][z]/4.0f;blackHole[x][y][z]=n*blackHole[x][y][z]/4.0f;fac*=n/4.0f;}
                    } else if(type[k].equals("OUTMASK")) {
                        float base=(1.0f-outlier[x][y][z])*num+outlier[x][y][z]*numO;
                        float outValue=outlier[x][y][z]*numO;
                        if(relation[k]==1){float f=factor(shape[x][y][z][k]);base=f*base;outValue=(f*outlier[x][y][z])*numO;}
                        mems[x][y][z][k]=base;outlier[x][y][z]=outValue;
                    } else mems[x][y][z][k]=num;
                } else {
                    mems[x][y][z][k]=0.0f;
                    if(type[k].equals("LESION")||type[k].equals("BLACKHOLE"))lesion[x][y][z]=0.0f;
                    if(type[k].equals("BLACKHOLE"))blackHole[x][y][z]=0.0f;
                }
                if(type[k].equals("LESION"))den+=lesion[x][y][z];
                else if(type[k].equals("BLACKHOLE")){float aux=lesion[x][y][z]+blackHole[x][y][z];den+=aux;}
                den+=mems[x][y][z][k];
            }
            if(den>0.0f)for(int k=0;k<nc;k++) {
                mems[x][y][z][k]/=den;
                if(type[k].equals("LESION")||type[k].equals("BLACKHOLE"))lesion[x][y][z]/=den;
                if(type[k].equals("BLACKHOLE"))blackHole[x][y][z]/=den;
                if(type[k].equals("OUTMASK"))outlier[x][y][z]/=den;
            }
            out.append("{\"coordinate\":[").append(x).append(',').append(y).append(',').append(z).append("],\"class_bits\":[");
            for(int k=0;k<nc;k++){if(k>0)out.append(',');out.append(bits(mems[x][y][z][k]));}
            out.append("],\"lesion_bits\":").append(bits(lesion[x][y][z])).append(",\"black_hole_bits\":").append(bits(blackHole[x][y][z])).append(",\"outlier_bits\":").append(bits(outlier[x][y][z])).append('}');
        }
        System.out.print(out.append("]}\n"));
    }
}
