/** Selected one-voxel image residuals through smoothing/prior membership normalization. */
public final class LesionToadsDataEnergyMembershipReferenceHarness {
    private static final float[][] NEIGHBORS={{0.10f,0.25f},{0.15f,0.30f},{0.20f,0.35f},{0.22f,0.37f},{0.12f,0.27f},{0.18f,0.33f}};
    private static final float[] NLESION={0.10f,0.20f,0.30f,0.40f,0.50f,0.60f};
    private static final float[] NBH={0.05f,0.15f,0.25f,0.35f,0.45f,0.55f};
    private static final float[][] SHAPE_PRIORS={{0.4f,0.6f},{0.25f,0.75f}};
    private static final float[][][] SIMILARITIES={{{0.9f,0.4f},{0.3f,0.8f}},{{0.7f,0.2f},{0.4f,0.9f}}};
    private static final float[][] LESION_SIM={{0.2f,0.7f},{0.3f,0.8f}};
    private static final float[][] BH_SIM={{0.6f,0.1f},{0.9f,0.2f}};
    private static final float[][] SHAPES={{0.45f,0.3f},{0.4f,0.2f}};
    private static final float[][] IMAGES={{0.35f,0.75f},{0.65f,0.25f}};
    private static final float[][][] CENTROIDS={{{0.1f,0.3f},{0.6f,0.2f}},{{0.2f,0.5f},{0.7f,0.15f}}};
    private static final float[][] LESION_CENTROIDS={{0.5f,0.45f},{0.3f,0.4f}};
    private static final float[][] BH_CENTROIDS={{0.8f,0.4f},{0.9f,0.6f}};
    private static final float[][] WEIGHTS={{1.2f,0.8f},{0.75f,1.25f}};
    private static final String[] TARGETS={"LESION","BLACKHOLE"};
    private static final float[] SMOOTHING={0.6f,0.25f},PRIOR={0.4f,0.8f},PRIOR_COEFFICIENT={0.35f,0.45f};
    private static final int[] VENT={2,-1},GM={0,0},BSTEM={3,1},MAX_VENT={4,9},MAX_GM={5,5},MAX_BSTEM={5,4};
    private static final float[] SPREAD={0.5f,0.375f};
    private static final float ZERO=1.0e-10f,INF=1.0e10f;
    private static long bits(float x){return Float.floatToIntBits(x)&0xffffffffL;}
    private static String vec(float[] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(bits(a[i]));}return s.append(']').toString();}
    private static String matrix(float[][] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(vec(a[i]));}return s.append(']').toString();}
    private static String cube(float[][][] a){StringBuilder s=new StringBuilder("[");for(int i=0;i<a.length;i++){if(i>0)s.append(',');s.append(matrix(a[i]));}return s.append(']').toString();}
    private static float sf(float x){float twice=2.0f*x;return 1.0f<twice?1.0f:twice;}
    private static float inv(float e,String rel,float shape){float n=e>ZERO?1.0f/e:INF;return rel.equals("prior")?sf(shape)*n:n;}
    private static float[] data(int scenario,int k){
        float ce=0,le=0,be=0;String type=k==0?"OBJECT":TARGETS[scenario];
        for(int c=0;c<2;c++){
            float d=IMAGES[scenario][c]-CENTROIDS[scenario][k][c];float wd=WEIGHTS[scenario][c]*d;ce+=wd*d;
            if(type.equals("LESION")||type.equals("BLACKHOLE")){
                float dl=IMAGES[scenario][c]-LESION_CENTROIDS[scenario][c];float wl=WEIGHTS[scenario][c]*dl;le+=wl*dl;
            }
            if(type.equals("BLACKHOLE")){
                float db=IMAGES[scenario][c]-BH_CENTROIDS[scenario][c];float wb=WEIGHTS[scenario][c]*db;be+=wb*db;
            }
        }
        return new float[]{ce,le,be};
    }
    private static float[] update(int scenario,int k,float[] raw){
        String type=k==0?"OBJECT":TARGETS[scenario],rel=k==0?"direct":(scenario==0?"direct":"prior");
        float p=PRIOR[scenario],smooth=SMOOTHING[scenario],pc=PRIOR_COEFFICIENT[scenario];
        float ngb=0,ngbL=0,ngbB=0,ngbS=0;
        for(int f=0;f<6;f++){for(int m=0;m<2;m++)if(m!=k)ngb+=NEIGHBORS[f][m]*NEIGHBORS[f][m];ngbL+=NLESION[f]*NLESION[f];ngbB+=NBH[f]*NBH[f];}
        if(type.equals("LESION")||type.equals("BLACKHOLE"))for(int f:new int[]{1,0,3,2,5,4})ngbS+=NEIGHBORS[f][k]*NEIGHBORS[f][k];
        float common=ngb+ngbL+ngbB;float classTerm=raw[0]+p*smooth*common;float classEnergy=raw[0]+0.5f*p*smooth*common;
        float lesionTerm=raw[1],bhTerm=raw[2],le=raw[1],be=raw[2];
        if(type.equals("LESION")||type.equals("BLACKHOLE")){float v=ngb+ngbS+ngbB;lesionTerm+=p*smooth*v;le+=0.5f*p*smooth*v;}
        if(type.equals("BLACKHOLE")){float v=ngb+ngbS+ngbL;bhTerm+=p*smooth*v;be+=0.5f*p*smooth*v;}
        float ip=0,il=0,ib=0;
        for(int m=0;m<2;m++){
            float shape=SHAPES[scenario][m];if(m!=k)ip+=SHAPE_PRIORS[scenario][m]*SIMILARITIES[scenario][k][m]*shape*shape;
            if(type.equals("LESION")||type.equals("BLACKHOLE"))il+=SHAPE_PRIORS[scenario][m]*LESION_SIM[scenario][m]*shape*shape;
            if(type.equals("BLACKHOLE"))ib+=SHAPE_PRIORS[scenario][m]*BH_SIM[scenario][m]*shape*shape;
        }
        ip+=SHAPE_PRIORS[scenario][1]*LESION_SIM[scenario][k]*SHAPES[scenario][1]*SHAPES[scenario][1];
        if(TARGETS[scenario].equals("BLACKHOLE"))ip+=SHAPE_PRIORS[scenario][1]*BH_SIM[scenario][k]*SHAPES[scenario][1]*SHAPES[scenario][1];
        classTerm+=pc*ip;classEnergy+=0.5f*pc*ip;
        if(type.equals("LESION")||type.equals("BLACKHOLE")){lesionTerm+=pc*il;le+=0.5f*pc*il;}
        if(type.equals("BLACKHOLE")){bhTerm+=pc*ib;be+=0.5f*pc*ib;}
        return new float[]{inv(classTerm,rel,SHAPES[scenario][k]),classEnergy,lesionTerm,le,bhTerm,be};
    }
    private static void run(int s,boolean first){
        String target=TARGETS[s],rel=s==0?"direct":"prior";String[] types={"OBJECT",target};float[] mem=new float[2],energy=new float[2];
        float lesion=0,bh=0,relation=0,le=0,be=0;float[][] rawData={data(s,0),data(s,1)};
        for(int k=0;k<2;k++){
            float[] v=update(s,k,rawData[k]);mem[k]=v[0];energy[k]=v[1];
            if(k==1){lesion=inv(v[2],rel,SHAPES[s][1]);le=v[3];be=v[5];if(target.equals("BLACKHOLE"))bh=inv(v[4],rel,SHAPES[s][1]);relation=rel.equals("direct")?1.0f:sf(SHAPES[s][1]);
                if(target.equals("LESION")){
                    if(MAX_VENT[s]!=0&&VENT[s]!=-1){float n=VENT[s]==0?SPREAD[s]:VENT[s];lesion=n*lesion/MAX_VENT[s];relation*=n/MAX_VENT[s];}
                    else if(MAX_GM[s]!=0&&GM[s]!=-1){float n=GM[s]==0?SPREAD[s]:GM[s];lesion=n*lesion/MAX_GM[s];relation*=n/MAX_GM[s];}
                    if(MAX_BSTEM[s]!=0&&BSTEM[s]!=-1){lesion=BSTEM[s]*lesion/MAX_BSTEM[s];relation*=(float)BSTEM[s]/MAX_BSTEM[s];}
                }else{
                    if(MAX_VENT[s]!=0&&VENT[s]!=-1){float n=VENT[s]==0?SPREAD[s]:VENT[s];lesion=n*lesion/MAX_VENT[s];relation*=n/MAX_VENT[s];}
                    else if(MAX_GM[s]!=0&&GM[s]!=-1){float n=GM[s]==0?SPREAD[s]:GM[s];lesion=n*lesion/MAX_GM[s];bh=n*bh/MAX_GM[s];relation*=n/MAX_GM[s];}
                    if(MAX_BSTEM[s]!=0&&BSTEM[s]!=-1){lesion=BSTEM[s]*lesion/MAX_BSTEM[s];bh=BSTEM[s]*bh/MAX_BSTEM[s];relation*=(float)BSTEM[s]/MAX_BSTEM[s];}
                }
            }
        }
        float den=0;for(int k=0;k<2;k++){if(types[k].equals("LESION"))den+=lesion;else if(types[k].equals("BLACKHOLE"))den+=lesion+bh;den+=mem[k];}
        float dist=0;if(den>0)for(int k=0;k<2;k++){mem[k]/=den;if(types[k].equals("LESION"))lesion/=den;else if(types[k].equals("BLACKHOLE")){lesion/=den;bh/=den;}dist+=energy[k]*mem[k]*mem[k];}
        dist+=le*lesion*lesion+be*bh*bh;
        if(!first)System.out.print(',');
        System.out.print("{\"target_type\":\""+target+"\",\"relation_type\":\""+rel+"\",\"image_bits\":"+vec(IMAGES[s])
            +",\"centroid_bits\":"+matrix(CENTROIDS[s])+",\"lesion_centroid_bits\":"+vec(LESION_CENTROIDS[s])
            +",\"black_hole_centroid_bits\":"+vec(BH_CENTROIDS[s])+",\"channel_weight_bits\":"+vec(WEIGHTS[s])
            +",\"shape_prior_bits\":"+vec(SHAPE_PRIORS[s])+",\"similarity_bits\":"+matrix(SIMILARITIES[s])
            +",\"lesion_similarity_bits\":"+vec(LESION_SIM[s])+",\"black_hole_similarity_bits\":"+vec(BH_SIM[s])
            +",\"shape_values_bits\":"+vec(SHAPES[s])+",\"smoothing_bits\":"+bits(SMOOTHING[s])+",\"prior_bits\":"+bits(PRIOR[s])
            +",\"prior_coefficient_bits\":"+bits(PRIOR_COEFFICIENT[s])+",\"ventricle_distance\":"+VENT[s]+",\"gray_matter_distance\":"+GM[s]
            +",\"brainstem_distance\":"+BSTEM[s]+",\"max_ventricle_distance\":"+MAX_VENT[s]+",\"max_gray_matter_distance\":"+MAX_GM[s]
            +",\"max_brainstem_distance\":"+MAX_BSTEM[s]+",\"spread_bits\":"+bits(SPREAD[s])+",\"neighbor_class_bits\":"+matrix(NEIGHBORS)
            +",\"neighbor_lesion_bits\":"+vec(NLESION)+",\"neighbor_black_hole_bits\":"+vec(NBH)
            +",\"class_data_bits\":"+matrix(rawData)+",\"class_memberships_bits\":"+vec(mem)+",\"class_energies_bits\":"+vec(energy)
            +",\"lesion_energy_bits\":"+bits(le)+",\"black_hole_energy_bits\":"+bits(be)+",\"lesion_membership_bits\":"+bits(lesion)
            +",\"black_hole_membership_bits\":"+bits(bh)+",\"lesion_relation_bits\":"+bits(relation)+",\"denominator_bits\":"+bits(den)
            +",\"local_energy_bits\":"+bits(dist)+"}");
    }
    public static void main(String[] args){
        System.out.print("{\"source\":\"VandyCRUISE LesionToadSegmentation.computeMemberships selected data-energy composition\",\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"cases\":[");
        run(0,true);run(1,false);System.out.println("]}");
    }
}
