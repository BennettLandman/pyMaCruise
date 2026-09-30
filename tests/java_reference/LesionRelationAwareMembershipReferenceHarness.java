/** Selected LesionTOADS relation dispatch used by computeMemberships().
 * The selected VandyCRUISE method resolves each class from
 * relationList[relationMap[x][y][z]][k], then applies the existing direct or
 * prior inverse-membership rule. Source SHA-256:
 * 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 */
public final class LesionRelationAwareMembershipReferenceHarness {
    private static boolean direct(boolean[][] relationList, short[][][] relationMap,
                                  int x, int y, int z, int k) {
        return relationList[relationMap[x][y][z]][k];
    }
    private static float membership(float energy, float shape, boolean direct) {
        float value = energy > 1.0e-10f ? 1.0f / energy : 1.0e10f;
        if (!direct) {
            float doubled = 2.0f * shape;
            float factor = 1.0f < doubled ? 1.0f : doubled;
            value = factor * value;
        }
        return value;
    }
    private static int bits(float value) { return Float.floatToRawIntBits(value); }
    private static String singleVoxel() {
        boolean[] direct={true,true,false,false};
        float[] shape={0.30f,0.25f,0.20f,0.45f};
        float[] raw={2.0f,3.0f,4.0f,5.0f};
        float rawL=6.0f, rawB=7.0f, rawO=8.0f, outlier=0.25f;
        float[] mem=new float[4]; float lesion=0.0f,bh=0.0f,out=0.0f,den=0.0f;
        float[] energy=raw.clone();
        energy[3]=(1.0f-outlier)*raw[3]+outlier*rawO;
        for(int k=0;k<4;k++){
            mem[k]=membership(raw[k],shape[k],direct[k]);
            if(k==1||k==2) lesion=membership(rawL,shape[k],direct[k]);
            if(k==2) bh=membership(rawB,shape[k],direct[k]);
            if(k==3){
                float classOut=membership(raw[3],1.0f,true), outMem=membership(rawO,1.0f,true);
                float factor=direct[k]?1.0f:(2.0f*shape[k]<1.0f?2.0f*shape[k]:1.0f);
                mem[k]=factor*((1.0f-outlier)*classOut+outlier*outMem);
                out=factor*outlier*outMem;
            }
            if(k==1) den+=lesion;
            else if(k==2) den+=lesion+bh;
            den+=mem[k];
        }
        for(int k=0;k<4;k++){
            mem[k]/=den;
            if(k==1||k==2) lesion/=den;
            if(k==2) bh/=den;
            else if(k==3) out/=den;
        }
        float objective=0.0f;
        for(int k=0;k<4;k++) objective+=energy[k]*mem[k]*mem[k];
        objective+=rawL*lesion*lesion+rawB*bh*bh;
        StringBuilder s=new StringBuilder("{\"class_bits\":[");
        for(int k=0;k<4;k++){if(k>0)s.append(',');s.append(bits(mem[k]));}
        return s.append("],\"lesion_bits\":").append(bits(lesion))
            .append(",\"black_hole_bits\":").append(bits(bh))
            .append(",\"outlier_bits\":").append(bits(out))
            .append(",\"denominator_bits\":").append(bits(den))
            .append(",\"objective_bits\":").append(bits(objective)).append('}').toString();
    }
    public static void main(String[] args) {
        boolean[][] relationList={{true,false,true,false},{true,true,false,false}};
        short[][][] relationMap={{{0,1},{1,0}},{{0,1},{1,0}}};
        boolean[] class1=new boolean[8]; int index=0;
        for(int x=0;x<2;x++)for(int y=0;y<2;y++)for(int z=0;z<2;z++)
            class1[index++]=direct(relationList,relationMap,x,y,z,1);
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"class1_dispatch\":[");
        for(int i=0;i<class1.length;i++){if(i>0)System.out.print(",");System.out.print(class1[i]);}
        System.out.println("],\"direct_bits\":"+Float.floatToRawIntBits(membership(2.0f,0.25f,true))+
            ",\"prior_bits\":"+Float.floatToRawIntBits(membership(2.0f,0.25f,false))+",\"single_voxel\":"+singleVoxel()+"}");
    }
}
