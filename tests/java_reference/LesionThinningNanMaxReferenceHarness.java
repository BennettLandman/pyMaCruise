import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.Numerics;

/**
 * Selected initThinningLabels/propagateThinning excerpt for the useLesionMax
 * NaN comparison path. Source identities:
 * LesionToadSegmentation.java 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99
 * BinaryHeap4D.java 9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a
 * Numerics.java 13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb
 */
public final class LesionThinningNanMaxReferenceHarness {
    private static final int NX=5, NY=5, NZ=5, CLASSES=13;
    private static final float ZERO=1.0e-10f;

    private static int label(int x) { return x < 2 ? 0 : 10; }

    private static float membership(int x, int y, int z, int k) {
        if ((k!=0 && k!=10) || (x==2 && y==2 && z==2 && k==10)) return Float.NaN;
        if (k==0) return (8 + (x + 2*y + 3*z) % 5) / 32.0f;
        return (18 + (2*x + y + z) % 7) / 32.0f;
    }

    private static float lesionScore(float ordinary) {
        // Selected Numerics.max(a,b,c) evaluates max(a,max(b,c)); for NaN a,
        // its strict greater-than comparison returns the finite inner result.
        return Numerics.max(ordinary, 0.0f, 0.0f);
    }

    private static int bits(float value) { return Float.floatToRawIntBits(value); }

    private static void append(StringBuilder out, String name, int[] values) {
        out.append("\"").append(name).append("\":[");
        for (int i=0;i<values.length;i++) {
            if (i>0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static String run() {
        float[] volume=new float[CLASSES];
        byte[] skeleton=new byte[NX*NY*NZ];
        boolean[] available=new boolean[NX*NY*NZ*CLASSES];
        int[] segmentation=new int[NX*NY*NZ];
        int[] boundaryCounts=new int[CLASSES];
        BinaryHeap4D tree=new BinaryHeap4D(NX*NY*NZ*CLASSES, BinaryHeap4D.MAXTREE);
        int[] insertedBits=new int[NX*NY*NZ*CLASSES];
        int[] insertedPayloads=new int[NX*NY*NZ*CLASSES*4];
        int inserted=0;
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            boolean border=x==0 || x==NX-1 || y==0 || y==NY-1 || z==0 || z==NZ-1;
            int offset=(x*NY+y)*NZ+z;
            segmentation[offset]=label(x);
            skeleton[offset]=(byte)(border?7:-1);
            for (int k=0;k<CLASSES;k++)
                available[offset*CLASSES+k]=((x+y+z+k)%2)==0;
        }
        int[] dx={-1,1,0,0,0,0}, dy={0,0,-1,1,0,0}, dz={0,0,0,0,-1,1};
        for (int x=1;x<NX-1;x++) for (int y=1;y<NY-1;y++) for (int z=1;z<NZ-1;z++) {
            int current=label(x);
            int offset=(x*NY+y)*NZ+z;
            skeleton[offset]=(byte)current;
            for (int k=0;k<CLASSES;k++) {
                available[offset*CLASSES+k]=current==k;
                if (current==k) volume[k]+=(k==10?lesionScore(membership(x,y,z,k)):membership(x,y,z,k));
            }
        }
        for (int x=1;x<NX-1;x++) for (int y=1;y<NY-1;y++) for (int z=1;z<NZ-1;z++) {
            int current=label(x);
            boolean boundary=false;
            for (int face=0;face<6;face++)
                if (label(x+dx[face])!=current) boundary=true;
            if (!boundary) continue;
            int offset=(x*NY+y)*NZ+z;
            available[offset*CLASSES+current]=false;
            int best=0; float max=0.0f;
            for (int k=0;k<CLASSES;k++) {
                float score=k==10?lesionScore(membership(x,y,z,k)):membership(x,y,z,k);
                if (score>max) { best=k; max=membership(x,y,z,k); }
            }
            float currentScore=current==10?lesionScore(membership(x,y,z,current)):membership(x,y,z,current);
            float bestScore=best==10?lesionScore(membership(x,y,z,best)):membership(x,y,z,best);
            float ratio=current==best?1.0f:(bestScore+ZERO)/(currentScore+ZERO);
            float priority=-1.0f/ratio;
            tree.addValue(priority,x,y,z,current);
            insertedBits[inserted]=bits(priority);
            int itemOffset=4*inserted;
            insertedPayloads[itemOffset]=x; insertedPayloads[itemOffset+1]=y;
            insertedPayloads[itemOffset+2]=z; insertedPayloads[itemOffset+3]=current;
            inserted++; boundaryCounts[current]++;
        }
        float[] minimum=new float[CLASSES];
        for (int k=0;k<CLASSES;k++) minimum[k]=Numerics.max(100.0f,0.4f*volume[k]);
        int[] volumeBits=new int[CLASSES], minimumBits=new int[CLASSES];
        for (int k=0;k<CLASSES;k++) { volumeBits[k]=bits(volume[k]); minimumBits[k]=bits(minimum[k]); }
        int[] skeletonInts=new int[skeleton.length];
        for (int i=0;i<skeleton.length;i++) skeletonInts[i]=skeleton[i];
        int[] availableInts=new int[available.length];
        for (int i=0;i<available.length;i++) availableInts[i]=available[i]?1:0;
        int[] popBits=new int[inserted], popPayloads=new int[4*inserted], processed=new int[inserted];
        int count=0; float val=0.0f;
        while (tree.isNotEmpty() && val>-1000.0f) {
            int x=tree.getFirstX(), y=tree.getFirstY(), z=tree.getFirstZ(), k=tree.getFirstK();
            val=tree.getFirst(); tree.removeFirst();
            popBits[count]=bits(val);
            int offset=4*count;
            popPayloads[offset]=x; popPayloads[offset+1]=y; popPayloads[offset+2]=z; popPayloads[offset+3]=k;
            processed[count]=volume[k]>minimum[k]?1:0;
            count++;
        }
        int[] remainingBits=new int[tree.getCurrentSize()];
        int[] remainingPayloads=new int[4*tree.getCurrentSize()];
        for (int i=0;i<remainingBits.length;i++) {
            remainingBits[i]=bits(tree.getFirst());
            int offset=4*i;
            remainingPayloads[offset]=tree.getFirstX(); remainingPayloads[offset+1]=tree.getFirstY();
            remainingPayloads[offset+2]=tree.getFirstZ(); remainingPayloads[offset+3]=tree.getFirstK();
            tree.removeFirst();
        }
        StringBuilder out=new StringBuilder("{");
        append(out,"volume_bits",volumeBits); out.append(',');
        append(out,"minimum_bits",minimumBits); out.append(',');
        append(out,"boundary_counts",boundaryCounts); out.append(',');
        append(out,"skeleton",skeletonInts); out.append(',');
        append(out,"available",availableInts); out.append(',');
        append(out,"inserted_bits",java.util.Arrays.copyOf(insertedBits,inserted)); out.append(',');
        append(out,"inserted_payloads",java.util.Arrays.copyOf(insertedPayloads,4*inserted)); out.append(',');
        append(out,"pop_bits",java.util.Arrays.copyOf(popBits,count)); out.append(',');
        append(out,"pop_payloads",java.util.Arrays.copyOf(popPayloads,4*count)); out.append(',');
        append(out,"processed",java.util.Arrays.copyOf(processed,count)); out.append(',');
        append(out,"remaining_bits",remainingBits); out.append(',');
        append(out,"remaining_payloads",remainingPayloads);
        return out.append('}').toString();
    }

    public static void main(String[] args) { System.out.println(run()); }
}
