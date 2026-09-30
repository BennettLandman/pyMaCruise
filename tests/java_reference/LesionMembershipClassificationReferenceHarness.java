/**
 * Selected VandyCRUISE LesionToadSegmentation membership-map excerpt.
 *
 * Pins the behavior of updateLesionClassifcationformMemberships: scan names
 * with equalsIgnoreCase (last Cerebellum-GM match wins; default index 0), then
 * apply the strict Java-float membership comparison for each voxel.
 */
public final class LesionMembershipClassificationReferenceHarness {
    private static final String SOURCE_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final int NX=3, NY=2, NZ=2, CLASSES=4;
    private static final byte[] CLASSIFICATION = {42,42,42,42,7,42,42,42,42,42,42,42};
    private static final byte[] TEMPLATE_LABELS = {0,42,43,44};
    // Voxel-major, class-minor values, equivalent to mems[x][y][z][k].
    private static final float[] MEMBERSHIPS = {
        .05f,.30f,.07f,.20f,  .10f,.29f,.30f,.20f,
        .10f,.31f,.10f,.20f,  .10f,.20f,.10f,.20f,
        .01f,.10f,.02f,.02f,  .10f,.10f,.10f,.10f,
        .10f,.10f,.20f,.40f,  .10f,.15f,.10f,.20f,
        .05f,.25f,.05f,.15f,  .10f,.30f,.10f,.10f,
        .10f,.24f,.10f,.20f,  .10f,.20f,.10f,.25f
    };
    private static final float[] LESION = {
        .10f,.10f,.10f,.10f, .10f,.10f,.10f,.10f,
        .10f,.10f,.10f,.10f
    };

    private static final class Result {
        final short[] map;
        final int volume;
        final int gm;
        Result(short[] map, int volume, int gm) { this.map=map; this.volume=volume; this.gm=gm; }
    }

    private static Result execute(String[] names, boolean hasLesion) {
        int gm=0;
        for (int k=0; k<CLASSES; k++)
            if (names[k].equalsIgnoreCase("Cerebellum-GM")) gm=k;
        short[] map = new short[NX*NY*NZ];
        int volume=0, i=0;
        for (int x=0; x<NX; x++) for (int y=0; y<NY; y++) for (int z=0; z<NZ; z++) {
            if (hasLesion && CLASSIFICATION[i] == TEMPLATE_LABELS[1]
                    && MEMBERSHIPS[i*CLASSES+1] < (LESION[i] + MEMBERSHIPS[i*CLASSES+gm])) {
                map[i]=(short)1;
                volume++;
            } else {
                map[i]=(short)0;
            }
            i++;
        }
        return new Result(map, volume, gm);
    }

    private static void emit(StringBuilder out, byte[] v) {
        out.append('['); for(int i=0;i<v.length;i++){if(i>0)out.append(',');out.append(v[i]);} out.append(']');
    }
    private static void emit(StringBuilder out, short[] v) {
        out.append('['); for(int i=0;i<v.length;i++){if(i>0)out.append(',');out.append(v[i]);} out.append(']');
    }
    private static void emit(StringBuilder out, float[] v) {
        out.append('['); for(int i=0;i<v.length;i++){if(i>0)out.append(',');out.append(Float.toString(v[i]));} out.append(']');
    }
    private static void emitCase(StringBuilder out, String name, String[] names, boolean hasLesion) {
        Result r=execute(names,hasLesion);
        out.append("{\"name\":\"").append(name).append("\",\"atlas_names\":[");
        for(int i=0;i<names.length;i++){if(i>0)out.append(',');out.append('"').append(names[i]).append('"');}
        out.append("],\"has_lesion\":").append(hasLesion).append(",\"gm_class_index\":")
            .append(r.gm).append(",\"map\":"); emit(out,r.map);
        out.append(",\"lesion_volume\":").append(r.volume).append('}');
    }
    public static void main(String[] args) {
        StringBuilder out=new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"source_method\":\"LesionToadSegmentation.updateLesionClassifcationformMemberships\"")
            .append(",\"shape\":[").append(NX).append(',').append(NY).append(',').append(NZ)
            .append("],\"lesion_class\":1,\"classification\":"); emit(out,CLASSIFICATION);
        out.append(",\"template_labels\":"); emit(out,TEMPLATE_LABELS);
        out.append(",\"memberships\":"); emit(out,MEMBERSHIPS);
        out.append(",\"lesion_membership\":"); emit(out,LESION);
        out.append(",\"cases\":[");
        emitCase(out,"duplicate_last_case_insensitive_match",new String[]{"Cerebellum-GM","Target","cErEbElLuM-gm","CEREBELLUM-GM"},true);
        out.append(',');
        emitCase(out,"absent_name_uses_default_zero",new String[]{"ReferenceZero","Target","Other","Elsewhere"},true);
        out.append(',');
        emitCase(out,"no_lesion_forces_zero",new String[]{"Cerebellum-GM","Target","cErEbElLuM-gm","CEREBELLUM-GM"},false);
        System.out.println(out.append("]}").toString());
    }
}
