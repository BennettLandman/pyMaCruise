import edu.jhmi.rad.medic.structures.BinaryHeap4D;
import edu.jhmi.rad.medic.utilities.Numerics;

/**
 * Selected-method excerpt of LesionToadSegmentation.initThinningLabels for
 * source-generated NaN memberships. The selected source hashes are pinned here:
 * LesionToadSegmentation.java 545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99
 * BinaryHeap4D.java 9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a
 * Numerics.java 13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb
 */
public final class LesionThinningSourceNanHarness {
    private static final int NX=5, NY=5, NZ=5, CLASSES=13;
    private static final float ZERO=1.0e-10f;

    private static int label(int x) { return x < 2 ? 0 : 10; }

    private static float membership(int x, int y, int z, int k) {
        if ((k!=0 && k!=10) || (x==2 && y==2 && z==2 && k==10)) return Float.NaN;
        if (k==0) return (8 + (x + 2*y + 3*z) % 5) / 32.0f;
        return (18 + (2*x + y + z) % 7) / 32.0f;
    }

    private static void ints(StringBuilder out, String name, int[] values) {
        out.append("\"").append(name).append("\":[");
        for (int i=0;i<values.length;i++) {
            if (i>0) out.append(',');
            out.append(values[i]);
        }
        out.append(']');
    }

    private static int bits(float value) { return Float.floatToRawIntBits(value); }

    private static String runSourcePass() {
        float[] volume=new float[CLASSES];
        byte[] skeleton=new byte[NX*NY*NZ];
        boolean[] available=new boolean[NX*NY*NZ*CLASSES];
        int[] segmentation=new int[NX*NY*NZ];
        BinaryHeap4D tree=new BinaryHeap4D(NX*NY*NZ*CLASSES, BinaryHeap4D.MAXTREE);
        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            boolean border=x==0 || x==NX-1 || y==0 || y==NY-1 || z==0 || z==NZ-1;
            skeleton[(x*NY+y)*NZ+z]=(byte)(border?7:-1);
            segmentation[(x*NY+y)*NZ+z]=label(x);
            for (int k=0;k<CLASSES;k++)
                available[((x*NY+y)*NZ+z)*CLASSES+k]=((x+y+z+k)%2)==0;
        }
        for (int x=1;x<NX-1;x++) for (int y=1;y<NY-1;y++) for (int z=1;z<NZ-1;z++) {
            int current=label(x);
            skeleton[(x*NY+y)*NZ+z]=(byte)current;
            for (int k=0;k<CLASSES;k++) {
                available[((x*NY+y)*NZ+z)*CLASSES+k]=current==k;
                if (current==k) volume[k]+=(k==10 ? (membership(x,y,z,k)+0.0f)+0.0f : membership(x,y,z,k));
            }
        }
        for (int x=1;x<NX-1;x++) for (int y=1;y<NY-1;y++) for (int z=1;z<NZ-1;z++) {
            int current=label(x);
            boolean boundary=false;
            for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++)
                if (dx*dx+dy*dy+dz*dz==1 && label(x+dx)!=current) boundary=true;
            if (!boundary) continue;
            available[((x*NY+y)*NZ+z)*CLASSES+current]=false;
            int best=0; float max=0.0f;
            for (int k=0;k<CLASSES;k++) {
                float score=membership(x,y,z,k);
                if (k==10) score=(membership(x,y,z,k)+0.0f)+0.0f;
                if (score>max) { best=k; max=membership(x,y,z,k); }
            }
            float cur=current==10?(membership(x,y,z,current)+0.0f)+0.0f:membership(x,y,z,current);
            float maximum=best==10?(membership(x,y,z,best)+0.0f)+0.0f:membership(x,y,z,best);
            float ratio=current==best?1.0f:(maximum+ZERO)/(cur+ZERO);
            tree.addValue(-1.0f/ratio,x,y,z,current);
        }
        float[] minimum=new float[CLASSES];
        for (int k=0;k<CLASSES;k++) minimum[k]=Numerics.max(100.0f,0.4f*volume[k]);
        StringBuilder poppedBits=new StringBuilder("[");
        StringBuilder poppedPayloads=new StringBuilder("[");
        StringBuilder processed=new StringBuilder("[");
        int count=0; float val=0.0f;
        while (tree.isNotEmpty() && val>-1000.0f) {
            int x=tree.getFirstX(), y=tree.getFirstY(), z=tree.getFirstZ(), k=tree.getFirstK();
            val=tree.getFirst(); tree.removeFirst();
            if (count>0) { poppedBits.append(','); poppedPayloads.append(','); processed.append(','); }
            poppedBits.append(bits(val));
            poppedPayloads.append(x).append(',').append(y).append(',').append(z).append(',').append(k);
            processed.append(volume[k]>minimum[k]);
            count++;
        }
        poppedBits.append(']'); poppedPayloads.append(']'); processed.append(']');
        int remaining=tree.getCurrentSize();
        StringBuilder remainingBits=new StringBuilder("[");
        StringBuilder remainingPayloads=new StringBuilder("[");
        for (int i=0;i<remaining;i++) {
            if (i>0) { remainingBits.append(','); remainingPayloads.append(','); }
            remainingBits.append(bits(tree.getFirst()));
            remainingPayloads.append(tree.getFirstX()).append(',').append(tree.getFirstY()).append(',')
                .append(tree.getFirstZ()).append(',').append(tree.getFirstK());
            tree.removeFirst();
        }
        remainingBits.append(']'); remainingPayloads.append(']');
        StringBuilder out=new StringBuilder("{\"minimum_bits\":[");
        for (int k=0;k<CLASSES;k++) { if (k>0) out.append(','); out.append(bits(minimum[k])); }
        out.append("],\"pop_bits\":").append(poppedBits).append(",\"pop_payloads\":").append(poppedPayloads)
            .append(",\"processed\":").append(processed).append(",\"remaining_bits\":").append(remainingBits)
            .append(",\"remaining_payloads\":").append(remainingPayloads).append('}');
        return out.toString();
    }

    public static void main(String[] args) {
        if (args.length>0 && args[0].equals("pass")) {
            System.out.println(runSourcePass());
            return;
        }
        float[] volume=new float[CLASSES];
        int[] counts=new int[CLASSES];
        byte[] skeleton=new byte[NX*NY*NZ];
        boolean[] available=new boolean[NX*NY*NZ*CLASSES];
        BinaryHeap4D tree=new BinaryHeap4D(NX*NY*NZ*CLASSES, BinaryHeap4D.MAXTREE);
        int[] insertionBits=new int[NX*NY*NZ*CLASSES];
        int[] insertionPayloads=new int[NX*NY*NZ*CLASSES*4];
        int inserted=0;

        for (int x=0;x<NX;x++) for (int y=0;y<NY;y++) for (int z=0;z<NZ;z++) {
            boolean border=x==0 || x==NX-1 || y==0 || y==NY-1 || z==0 || z==NZ-1;
            skeleton[(x*NY+y)*NZ+z]=(byte)(border?7:-1);
            for (int k=0;k<CLASSES;k++)
                available[((x*NY+y)*NZ+z)*CLASSES+k]=((x+y+z+k)%2)==0;
        }

        for (int x=1;x<NX-1;x++) for (int y=1;y<NY-1;y++) for (int z=1;z<NZ-1;z++) {
            int current=label(x);
            skeleton[(x*NY+y)*NZ+z]=(byte)current;
            for (int k=0;k<CLASSES;k++) {
                available[((x*NY+y)*NZ+z)*CLASSES+k]=current==k;
                if (current==k) volume[k]+=(k==10 ? (membership(x,y,z,k)+0.0f)+0.0f : membership(x,y,z,k));
            }
        }
        int[] dx={-1,1,0,0,0,0}, dy={0,0,-1,1,0,0}, dz={0,0,0,0,-1,1};
        for (int x=1;x<NX-1;x++) for (int y=1;y<NY-1;y++) for (int z=1;z<NZ-1;z++) {
            int current=label(x);
            boolean boundary=false;
            for (int f=0;f<6;f++) if (label(x+dx[f])!=current) boundary=true;
            if (!boundary) continue;
            available[((x*NY+y)*NZ+z)*CLASSES+current]=false;
            int best=0;
            float max=0.0f;
            for (int k=0;k<CLASSES;k++) {
                float score=membership(x,y,z,k);
                if (k==10) score=(membership(x,y,z,k)+0.0f)+0.0f;
                if (score>max) { best=k; max=membership(x,y,z,k); }
            }
            float cur=current==10?(membership(x,y,z,current)+0.0f)+0.0f:membership(x,y,z,current);
            float maximum=best==10?(membership(x,y,z,best)+0.0f)+0.0f:membership(x,y,z,best);
            float ratio;
            if (current==best) ratio=1.0f;
            else ratio=(maximum+ZERO)/(cur+ZERO);
            float priority=-1.0f/ratio;
            tree.addValue(priority,x,y,z,current);
            insertionBits[inserted]=bits(priority);
            int off=4*inserted;
            insertionPayloads[off]=x; insertionPayloads[off+1]=y;
            insertionPayloads[off+2]=z; insertionPayloads[off+3]=current;
            inserted++;
            counts[current]++;
        }

        int[] removalBits=new int[inserted];
        int[] removalPayloads=new int[4*inserted];
        for (int i=0;i<inserted;i++) {
            removalBits[i]=bits(tree.getFirst());
            int off=4*i;
            removalPayloads[off]=tree.getFirstX(); removalPayloads[off+1]=tree.getFirstY();
            removalPayloads[off+2]=tree.getFirstZ(); removalPayloads[off+3]=tree.getFirstK();
            tree.removeFirst();
        }
        int[] volumeBits=new int[CLASSES];
        for (int i=0;i<CLASSES;i++) volumeBits[i]=bits(volume[i]);
        StringBuilder out=new StringBuilder("{");
        ints(out,"volume_bits",volumeBits); out.append(',');
        ints(out,"boundary_counts",counts); out.append(',');
        int[] skeletonInts=new int[skeleton.length];
        for (int i=0;i<skeleton.length;i++) skeletonInts[i]=skeleton[i];
        ints(out,"skeleton",skeletonInts); out.append(',');
        int[] availableInts=new int[available.length];
        for (int i=0;i<available.length;i++) availableInts[i]=available[i]?1:0;
        ints(out,"available",availableInts); out.append(',');
        ints(out,"inserted_priority_bits",java.util.Arrays.copyOf(insertionBits,inserted)); out.append(',');
        ints(out,"inserted_payloads",java.util.Arrays.copyOf(insertionPayloads,4*inserted)); out.append(',');
        ints(out,"removed_priority_bits",removalBits); out.append(',');
        ints(out,"removed_payloads",removalPayloads);
        System.out.println(out.append('}'));
    }
}
