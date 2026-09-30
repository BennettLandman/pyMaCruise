import edu.jhmi.rad.medic.utilities.PowerTable;

/*
 * Selected VandyCRUISE LesionToadSegmentation.computeMemberships() energy
 * excerpt, source SHA-256 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99.
 * PowerTable.java SHA-256 2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae.
 * The excerpt retains x/y/z, class, and channel loops, field multiplication,
 * the Sulcal-CSF block/comment behavior, lesion FLAIR weighting, and shared
 * auxiliary energy assignment. It intentionally omits spatial/update terms.
 */
public final class LesionDataEnergyMapsSpecialReferenceHarness {
    static final int NX = 7, NY = 8, NZ = 9, NC = 3, K = 4;
    static final int OBJECT = 1, LESION = 4, BLACKHOLE = 6;
    static final float[] WEIGHTS = {0.5f, 1.25f, 0.75f};
    static final int[] TYPES = {OBJECT, OBJECT, BLACKHOLE, LESION};
    static final float CSF_EXPONENT = 0.61f, LESION_EXPONENT = 0.73f;

    static final class Case {
        final String name, sulcalName;
        final boolean useField, hasMprage, useLesionWeight;
        final int t1, flair;
        Case(String name, String sulcalName, boolean useField, boolean hasMprage,
             int t1, int flair, boolean useLesionWeight) {
            this.name=name; this.sulcalName=sulcalName; this.useField=useField;
            this.hasMprage=hasMprage; this.t1=t1; this.flair=flair;
            this.useLesionWeight=useLesionWeight;
        }
    }

    static float image(int x, int y, int z, int c) {
        int fold = (3*x + 2*y + z) % 5 - 2;
        return (13*x - 9*y + 7*z + 5*fold + 11*c - 20) / 32.0f;
    }
    static float field(int x, int y, int z, int c) {
        return 0.625f + ((x + 2*y + 3*z + c) % 7) / 16.0f;
    }
    static float centroid(int k, int c) { return (2 + 5*k - 3*c) / 16.0f; }
    static float lesionCentroid(int c) { return (3 + 4*c) / 16.0f; }
    static float blackHoleCentroid(int c) { return (-5 + 2*c) / 16.0f; }
    static int bits(float value) { return Float.floatToIntBits(value); }
    static String name(int type) {
        if (type == OBJECT) return "OBJECT";
        if (type == LESION) return "LESION";
        if (type == BLACKHOLE) return "BLACKHOLE";
        throw new IllegalArgumentException("unsupported class type");
    }
    static void ints(StringBuilder out, int[] values) {
        out.append('[');
        for (int i=0; i<values.length; i++) { if (i>0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    static int[] centroidBits() {
        int[] result = new int[K*NC];
        for (int k=0;k<K;k++) for(int c=0;c<NC;c++) result[k*NC+c]=bits(centroid(k,c));
        return result;
    }
    static String run(Case test) {
        int voxels=NX*NY*NZ;
        int[] classes=new int[voxels*K], lesionByClass=new int[voxels*K], bhByClass=new int[voxels*K];
        int[] lesionOut=new int[voxels], bhOut=new int[voxels];
        PowerTable csfTable=new PowerTable(0.0f,4.0f,0.000001f,CSF_EXPONENT);
        PowerTable lesionTable=new PowerTable(0.0f,4.0f,0.000001f,LESION_EXPONENT);
        for(int x=1;x<NX-1;x++) for(int y=1;y<NY-1;y++) for(int z=1;z<NZ-1;z++) {
            int voxel=(x*NY+y)*NZ+z;
            float lesionEnergy=0.0f,blackHoleEnergy=0.0f;
            for(int k=0;k<K;k++) {
                float num=0.0f,numL=0.0f,numB=0.0f;
                if(k==1 && test.t1!=-1 && test.sulcalName.equals("Sulcal-CSF")) {
                    float value=test.useField ? field(x,y,z,test.t1)*image(x,y,z,test.t1) : image(x,y,z,test.t1);
                    float difference=value-centroid(k,test.t1), squared=difference*difference;
                    num += test.hasMprage ? csfTable.lookup(squared,CSF_EXPONENT) : squared;
                    // The source's generic-loop else is commented out in the
                    // unfielded branch, but active in the fielded branch.
                    if(!test.useField) {
                        for(int c=0;c<NC;c++) {
                            float d=image(x,y,z,c)-centroid(k,c);
                            num += WEIGHTS[c]*d*d;
                        }
                    }
                } else if(k==0 || k==1 || k==2 || k==3) {
                    for(int c=0;c<NC;c++) {
                        float value=test.useField ? field(x,y,z,c)*image(x,y,z,c) : image(x,y,z,c);
                        float d=value-centroid(k,c);
                        num += WEIGHTS[c]*d*d;
                    }
                }
                if(TYPES[k]==BLACKHOLE) {
                    for(int c=0;c<NC;c++) {
                        float value=test.useField ? field(x,y,z,c)*image(x,y,z,c) : image(x,y,z,c);
                        float dl=value-lesionCentroid(c), db=value-blackHoleCentroid(c);
                        numL += WEIGHTS[c]*dl*dl;
                        numB += WEIGHTS[c]*db*db;
                    }
                } else if(TYPES[k]==LESION) {
                    if(test.useLesionWeight && test.flair!=-1) {
                        float value=test.useField ? field(x,y,z,test.flair)*image(x,y,z,test.flair) : image(x,y,z,test.flair);
                        float d=value-lesionCentroid(test.flair);
                        numL += lesionTable.lookup(d*d,LESION_EXPONENT);
                    } else {
                        for(int c=0;c<NC;c++) {
                            float value=test.useField ? field(x,y,z,c)*image(x,y,z,c) : image(x,y,z,c);
                            float d=value-lesionCentroid(c);
                            numL += WEIGHTS[c]*d*d;
                        }
                    }
                }
                classes[voxel*K+k]=bits(num);
                if(TYPES[k]==LESION || TYPES[k]==BLACKHOLE) {
                    lesionByClass[voxel*K+k]=bits(numL);
                    lesionEnergy=numL;
                }
                if(TYPES[k]==BLACKHOLE) {
                    bhByClass[voxel*K+k]=bits(numB);
                    blackHoleEnergy=numB;
                }
            }
            lesionOut[voxel]=bits(lesionEnergy); bhOut[voxel]=bits(blackHoleEnergy);
        }
        StringBuilder out=new StringBuilder("{\"name\":\"").append(test.name).append("\",\"use_field\":").append(test.useField)
            .append(",\"has_mprage\":").append(test.hasMprage).append(",\"t1_channel\":").append(test.t1)
            .append(",\"flair_channel\":").append(test.flair).append(",\"use_lesion_weight\":").append(test.useLesionWeight)
            .append(",\"sulcal_name\":\"").append(test.sulcalName).append("\",\"class_types\":[");
        for(int k=0;k<K;k++){if(k>0)out.append(',');out.append('\"').append(name(TYPES[k])).append('\"');}
        out.append("],\"structure_names\":[\"Cerebrum-GM\",\"").append(test.sulcalName).append("\",\"Blackhole\",\"Lesion\"]");
        out.append(",\"class_centroid_bits\":");ints(out,centroidBits());
        out.append(",\"lesion_centroid_bits\":[");for(int c=0;c<NC;c++){if(c>0)out.append(',');out.append(bits(lesionCentroid(c)));}out.append(']');
        out.append(",\"black_hole_centroid_bits\":[");for(int c=0;c<NC;c++){if(c>0)out.append(',');out.append(bits(blackHoleCentroid(c)));}out.append(']');
        out.append(",\"channel_weight_bits\":[");for(int c=0;c<NC;c++){if(c>0)out.append(',');out.append(bits(WEIGHTS[c]));}out.append(']');
        out.append(",\"csf_power_exponent_bits\":").append(bits(CSF_EXPONENT)).append(",\"lesion_power_exponent_bits\":").append(bits(LESION_EXPONENT));
        out.append(",\"class_energy_bits\":");ints(out,classes);out.append(",\"lesion_energy_by_class_bits\":");ints(out,lesionByClass);
        out.append(",\"black_hole_energy_by_class_bits\":");ints(out,bhByClass);out.append(",\"lesion_energy_bits\":");ints(out,lesionOut);
        out.append(",\"black_hole_energy_bits\":");ints(out,bhOut);return out.append('}').toString();
    }

    public static void main(String[] args) {
        Case[] cases={
            new Case("unfielded_sulcal_direct","Sulcal-CSF",false,false,0,-1,false),
            new Case("unfielded_sulcal_mprage","Sulcal-CSF",false,true,0,-1,false),
            new Case("fielded_sulcal_direct","Sulcal-CSF",true,false,0,-1,false),
            new Case("fielded_sulcal_mprage","Sulcal-CSF",true,true,0,-1,false),
            new Case("sulcal_no_t1","Sulcal-CSF",false,true,-1,-1,false),
            new Case("sulcal_case_mismatch","sulcal-csf",false,true,0,-1,false),
            new Case("lesion_flair_unweighted","Cerebrum-GM",false,false,0,-1,true),
            new Case("lesion_flair_power_unfielded","Cerebrum-GM",false,false,0,2,true),
            new Case("lesion_flair_power_fielded","Cerebrum-GM",true,false,0,2,true)
        };
        System.out.print("{\"source_sha256\":\"545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99\",\"power_table_sha256\":\"2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae\",\"shape\":["+NX+","+NY+","+NZ+"],\"channels\":"+NC+",\"cases\":[");
        for(int i=0;i<cases.length;i++){if(i>0)System.out.print(',');System.out.print(run(cases[i]));}
        System.out.println("]}");
    }
}
