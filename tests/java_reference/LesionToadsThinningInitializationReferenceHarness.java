import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Selected-method excerpt of LesionToadSegmentation.initThinningLabels. */
public final class LesionToadsThinningInitializationReferenceHarness {
    private static final String SOURCE_SHA =
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99";
    private static final String HEAP_SOURCE_SHA =
        "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a";
    private static final int NX = 4, NY = 4, NZ = 4, CLASSES = 3;
    private static final float ZERO = 1.0e-10f;

    private static final class Inputs {
        final byte[][][] segmentation = new byte[NX][NY][NZ];
        final float[][][][] memberships = new float[NX][NY][NZ][CLASSES];
        final float[][][] lesion = new float[NX][NY][NZ];
        final float[][][] blackHole = new float[NX][NY][NZ];
        final byte[][][] skeleton = new byte[NX][NY][NZ];
        final boolean[][][][] available = new boolean[NX][NY][NZ][CLASSES];

        Inputs() {
            for (int x=0; x<NX; x++) for (int y=0; y<NY; y++) for (int z=0; z<NZ; z++) {
                int label = (x + y + z) % CLASSES;
                if ((x==0 && y==1 && z==1) || (x==2 && y==1 && z==1)
                        || (x==1 && y==0 && z==1) || (x==1 && y==2 && z==1)
                        || (x==1 && y==1 && z==0) || (x==1 && y==1 && z==2)) label = 0;
                segmentation[x][y][z] = (byte)label;
                skeleton[x][y][z] = (byte)(x==0 || x==NX-1 || y==0 || y==NY-1 || z==0 || z==NZ-1 ? 7 : -1);
                lesion[x][y][z] = (1 + ((x*3 + y*5 + z*7) % 9)) / 32.0f;
                blackHole[x][y][z] = (1 + ((x*5 + y*3 + z*2) % 7)) / 64.0f;
                for (int k=0; k<CLASSES; k++) {
                    memberships[x][y][z][k] =
                        (2 + ((x*7 + y*5 + z*3 + k*11) % 15)) / 32.0f;
                    available[x][y][z][k] = ((x + y + z + k) % 2) == 0;
                }
            }
        }
    }

    private static float lesionScore(float ordinary, float lesion, float blackHole,
            boolean useLesionMax) {
        if (useLesionMax) {
            float inner = lesion > blackHole ? lesion : blackHole;
            return ordinary > inner ? ordinary : inner;
        }
        return (ordinary + lesion) + blackHole;
    }

    private static int bestMembership(Inputs in, int x, int y, int z, boolean useLesionMax) {
        int id = 0;
        float max = 0.0f, score = 0.0f;
        for (int m=0; m<CLASSES; m++) {
            if (m==1 || m==2)
                score = lesionScore(in.memberships[x][y][z][m], in.lesion[x][y][z],
                    in.blackHole[x][y][z], useLesionMax);
            else score = in.memberships[x][y][z][m];
            if (score > max) {
                id = m;
                // Source quirk: compare adjusted score, but retain ordinary membership.
                max = in.memberships[x][y][z][m];
            }
        }
        return id;
    }

    private static float maxMembershipRatio(Inputs in, int x, int y, int z,
            int current, int best, boolean useLesionMax) {
        if (current == best) return 1.0f;
        float cur = current==1 || current==2
            ? lesionScore(in.memberships[x][y][z][current], in.lesion[x][y][z], in.blackHole[x][y][z], useLesionMax)
            : in.memberships[x][y][z][current];
        float max = best==1 || best==2
            ? lesionScore(in.memberships[x][y][z][best], in.lesion[x][y][z], in.blackHole[x][y][z], useLesionMax)
            : in.memberships[x][y][z][best];
        return (max + ZERO) / (cur + ZERO);
    }

    private static final class Result {
        final float[] volume;
        final int[] boundaryCounts;
        final byte[] skeleton;
        final boolean[] available;
        final float[] insertionPriorities;
        final int[] insertionPayloads;
        final float[] removalPriorities;
        final int[] removalPayloads;
        Result(float[] v, int[] c, byte[] s, boolean[] a, float[] ip, int[] ix,
                float[] rp, int[] rx) {
            volume=v; boundaryCounts=c; skeleton=s; available=a;
            insertionPriorities=ip; insertionPayloads=ix;
            removalPriorities=rp; removalPayloads=rx;
        }
    }

    private static Result execute(boolean useLesionMax) {
        Inputs in = new Inputs();
        float[] volume = new float[CLASSES];
        int[] boundaryCounts = new int[CLASSES];
        float[] insertedPriority = new float[NX*NY*NZ*CLASSES];
        int[] insertedPayload = new int[NX*NY*NZ*CLASSES*4];
        int inserted = 0;
        for (int x=1; x<NX-1; x++) for (int y=1; y<NY-1; y++) for (int z=1; z<NZ-1; z++) {
            int current = in.segmentation[x][y][z];
            in.skeleton[x][y][z] = (byte)current;
            for (int k=0; k<CLASSES; k++) {
                if (current==k) {
                    float value = k==1
                        ? lesionScore(in.memberships[x][y][z][k], in.lesion[x][y][z], in.blackHole[x][y][z], useLesionMax)
                        : in.memberships[x][y][z][k];
                    volume[k] += value;
                }
                in.available[x][y][z][k] = current==k;
            }
        }
        BinaryHeap4D tree = new BinaryHeap4D(NX*NY*NZ*CLASSES, BinaryHeap4D.MAXTREE);
        for (int x=1; x<NX-1; x++) for (int y=1; y<NY-1; y++) for (int z=1; z<NZ-1; z++) {
            int current = in.segmentation[x][y][z];
            boolean boundary = in.segmentation[x-1][y][z]!=current || in.segmentation[x+1][y][z]!=current
                || in.segmentation[x][y-1][z]!=current || in.segmentation[x][y+1][z]!=current
                || in.segmentation[x][y][z-1]!=current || in.segmentation[x][y][z+1]!=current;
            if (!boundary) continue;
            in.available[x][y][z][current] = false;
            int best = bestMembership(in,x,y,z,useLesionMax);
            float ratio = maxMembershipRatio(in,x,y,z,current,best,useLesionMax);
            float order = -1.0f / ratio;
            tree.addValue(order,x,y,z,current);
            insertedPriority[inserted] = order;
            int offset = inserted*4;
            insertedPayload[offset] = x; insertedPayload[offset+1] = y;
            insertedPayload[offset+2] = z; insertedPayload[offset+3] = current;
            inserted++;
            boundaryCounts[current]++;
        }
        float[] removalPriority = new float[inserted];
        int[] removalPayload = new int[inserted*4];
        for (int i=0; i<inserted; i++) {
            removalPriority[i] = tree.getFirst();
            int offset = i*4;
            removalPayload[offset] = tree.getFirstX();
            removalPayload[offset+1] = tree.getFirstY();
            removalPayload[offset+2] = tree.getFirstZ();
            removalPayload[offset+3] = tree.getFirstK();
            tree.removeFirst();
        }
        byte[] skeletonFlat = new byte[NX*NY*NZ];
        boolean[] availableFlat = new boolean[NX*NY*NZ*CLASSES];
        int voxel=0, ai=0;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            skeletonFlat[voxel++] = in.skeleton[x][y][z];
            for (int k=0;k<CLASSES;k++) availableFlat[ai++] = in.available[x][y][z][k];
        }
        float[] ip = new float[inserted];
        int[] ix = new int[inserted*4];
        System.arraycopy(insertedPriority,0,ip,0,inserted);
        System.arraycopy(insertedPayload,0,ix,0,inserted*4);
        return new Result(volume,boundaryCounts,skeletonFlat,availableFlat,ip,ix,removalPriority,removalPayload);
    }

    private static void emitFloats(StringBuilder out, float[] values) {
        out.append('[');
        for (int i=0;i<values.length;i++) { if (i>0) out.append(','); out.append(Float.floatToRawIntBits(values[i])); }
        out.append(']');
    }
    private static void emitInts(StringBuilder out, int[] values) {
        out.append('[');
        for (int i=0;i<values.length;i++) { if (i>0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    private static void emitBools(StringBuilder out, boolean[] values) {
        out.append('[');
        for (int i=0;i<values.length;i++) { if (i>0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    private static void emitBytes(StringBuilder out, byte[] values) {
        out.append('[');
        for (int i=0;i<values.length;i++) { if (i>0) out.append(','); out.append(values[i]); }
        out.append(']');
    }
    private static void emitCase(StringBuilder out, boolean useLesionMax) {
        Result r=execute(useLesionMax);
        out.append("{\"use_lesion_max\":").append(useLesionMax).append(",\"volume_bits\":"); emitFloats(out,r.volume);
        out.append(",\"boundary_counts\":"); emitInts(out,r.boundaryCounts);
        out.append(",\"skeleton\":"); emitBytes(out,r.skeleton);
        out.append(",\"available\":"); emitBools(out,r.available);
        out.append(",\"inserted_priority_bits\":"); emitFloats(out,r.insertionPriorities);
        out.append(",\"inserted_payloads\":"); emitInts(out,r.insertionPayloads);
        out.append(",\"removed_priority_bits\":"); emitFloats(out,r.removalPriorities);
        out.append(",\"removed_payloads\":"); emitInts(out,r.removalPayloads);
        out.append('}');
    }
    public static void main(String[] args) {
        StringBuilder out = new StringBuilder("{\"source_sha256\":\"").append(SOURCE_SHA)
            .append("\",\"heap_source_sha256\":\"").append(HEAP_SOURCE_SHA)
            .append("\",\"shape\":[4,4,4],\"cases\":[");
        emitCase(out,true); out.append(','); emitCase(out,false);
        System.out.println(out.append("]}").toString());
    }
}
