import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.zip.GZIPInputStream;

/**
 * Fixed two-pass excerpt of AlgorithmLesionToads.runAlgorithm's active loop.
 * Energy samples are from the P39v/P39w synthetic composition fixture; this
 * harness pins only Java float32 loop controls and source call order. The
 * numeric methods retain their own selected-source harnesses.
 * VandyCRUISE AlgorithmLesionToads.java SHA-256:
 * 556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb.
 */
public final class LesionToadsRepeatedIterationReferenceHarness {
    private static float ratio(float previous, float current) {
        return Math.abs(2.0f * (previous - current) / (previous + current));
    }
    private static String quote(String value) { return "\"" + value + "\""; }
    private static String sha256(byte[] bytes) throws Exception {
        byte[] digest=MessageDigest.getInstance("SHA-256").digest(bytes);
        StringBuilder result=new StringBuilder();
        for(byte value:digest)result.append(String.format("%02x",value&255));
        return result.toString();
    }
    private static String volumeBits(byte[] gzipped) throws Exception {
        ByteArrayOutputStream stream=new ByteArrayOutputStream();
        try(GZIPInputStream input=new GZIPInputStream(new ByteArrayInputStream(gzipped))) {
            input.transferTo(stream);
        }
        byte[] nifti=stream.toByteArray();
        ByteBuffer header=ByteBuffer.wrap(nifti).order(ByteOrder.LITTLE_ENDIAN);
        if(header.getInt(0)!=348||Short.toUnsignedInt(header.getShort(40))!=4
                ||Short.toUnsignedInt(header.getShort(70))!=16)
            throw new IllegalArgumentException("expected NIfTI-1 4-D float32 state fixture");
        int nx=Short.toUnsignedInt(header.getShort(42));
        int ny=Short.toUnsignedInt(header.getShort(44));
        int nz=Short.toUnsignedInt(header.getShort(46));
        int channels=Short.toUnsignedInt(header.getShort(48));
        int offset=(int)header.getFloat(108);
        if(channels!=14||offset<348||nifti.length<offset+nx*ny*nz*channels*4)
            throw new IllegalArgumentException("invalid P39x loop state dimensions or payload");
        ByteBuffer voxels=ByteBuffer.wrap(nifti).order(ByteOrder.LITTLE_ENDIAN);
        StringBuilder out=new StringBuilder("[");
        for(int pass=0;pass<2;pass++) {
            float[] volumes=new float[4];
            int base=pass*7;
            for(int x=1;x<nx-1;x++)for(int y=1;y<ny-1;y++)for(int z=1;z<nz-1;z++) {
                int voxel=((z*ny+y)*nx+x);
                int channelStride=nx*ny*nz;
                int current=(int)voxels.getFloat(offset+(base*channelStride+voxel)*4);
                if(current<0||current>=4)throw new IllegalArgumentException("invalid class label in loop state");
                float increment=voxels.getFloat(offset+((base+1+current)*channelStride+voxel)*4);
                if(current==1) {
                    increment+=voxels.getFloat(offset+((base+5)*channelStride+voxel)*4);
                    increment+=voxels.getFloat(offset+((base+6)*channelStride+voxel)*4);
                }
                volumes[current]+=increment;
            }
            if(pass>0)out.append(',');
            out.append('[');
            for(int k=0;k<4;k++){if(k>0)out.append(',');out.append(Float.floatToRawIntBits(volumes[k]));}
            out.append(']');
        }
        return out.append(']').toString();
    }
    public static void main(String[] args) throws Exception {
        if(args.length!=1)throw new IllegalArgumentException("supply the generated RAS loop-state NIfTI");
        byte[] loopStates;
        try { loopStates=Files.readAllBytes(Path.of(args[0])); }
        catch(Exception error) { throw new RuntimeException(error); }
        Locale.setDefault(Locale.ROOT);
        float[] energies = new float[] {49.18265151977539f, 26.26949119567871f, 25.682783126831055f};
        int nIterMax = 2;
        float maxDistance = 0.0f;
        List<Integer> ratios = new ArrayList<>();
        List<String> reasons = new ArrayList<>();
        List<String> stages = new ArrayList<>();
        for (int iter=1; iter<=Math.abs(nIterMax); iter++) {
            float ePrevious = energies[iter-1];
            float energy = energies[iter];
            float eRatio = ratio(ePrevious, energy);
            ratios.add(Float.floatToRawIntBits(eRatio));
            boolean stop = false;
            if (eRatio < maxDistance) stop = true;
            if (iter >= Math.abs(nIterMax)) stop = true;
            boolean doGrowing = !stop || nIterMax > 0;
            StringBuilder trace = new StringBuilder("[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\"");
            if (doGrowing) trace.append(",\"propagateGrowing\"");
            if (stop) trace.append(",\"updateClassificationFromLabels\",\"updateLesionClassifcation\"]");
            else trace.append(",\"computeRelations\"]");
            stages.add(trace.toString());
            reasons.add(stop ? (eRatio < maxDistance ? "convergence" : "iteration_limit") : "continue");
            if (stop) break;
        }
        String loopVolumeBits;
        try { loopVolumeBits=volumeBits(loopStates); }
        catch(Exception error) { throw new RuntimeException(error); }
        System.out.print("{\"source_sha256\":\"556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb\",\"loop_state_sha256\":\""+sha256(loopStates)+"\",\"thinning_initial_volume_bits\":"+loopVolumeBits+",\"energy_bits\":[");
        for (int i=0;i<energies.length;i++) { if (i>0) System.out.print(","); System.out.print(Float.floatToRawIntBits(energies[i])); }
        System.out.print("],\"ratio_bits\":[");
        for (int i=0;i<ratios.size();i++) { if (i>0) System.out.print(","); System.out.print(ratios.get(i)); }
        System.out.print("],\"stop_reasons\":[");
        for (int i=0;i<reasons.size();i++) { if (i>0) System.out.print(","); System.out.print(quote(reasons.get(i))); }
        System.out.print("],\"stage_order\":[");
        for (int i=0;i<stages.size();i++) { if (i>0) System.out.print(","); System.out.print(stages.get(i)); }
        float convergenceThreshold=0.7f;
        float firstRatio=ratio(energies[0],energies[1]);
        boolean convergenceStop=firstRatio<convergenceThreshold;
        boolean positiveGrowing=!convergenceStop||2>0;
        boolean negativeGrowing=!convergenceStop||-2>0;
        String positiveTrace="[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\""+(positiveGrowing?",\"propagateGrowing\"":"")+",\"updateClassificationFromLabels\",\"updateLesionClassifcation\"]";
        String negativeTrace="[\"propagateThinning\",\"computeCentroids\",\"computeVariances\",\"computeChannelWeights\",\"computeMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\""+(negativeGrowing?",\"propagateGrowing\"":"")+",\"updateClassificationFromLabels\",\"updateLesionClassifcation\"]";
        System.out.print("],\"negative_ceiling\":{\"n_iter_max\":-2,\"stop\":true,\"run_growing\":false,\"stop_reason\":\"iteration_limit\"},\"convergence_cases\":[");
        System.out.print("{\"n_iter_max\":2,\"max_distance_bits\":"+Float.floatToRawIntBits(convergenceThreshold)+",\"ratio_bits\":"+Float.floatToRawIntBits(firstRatio)+",\"passes\":1,\"stop\":"+convergenceStop+",\"run_growing\":"+positiveGrowing+",\"stop_reason\":\"convergence\",\"stage_order\":"+positiveTrace+"},");
        System.out.print("{\"n_iter_max\":-2,\"max_distance_bits\":"+Float.floatToRawIntBits(convergenceThreshold)+",\"ratio_bits\":"+Float.floatToRawIntBits(firstRatio)+",\"passes\":1,\"stop\":"+convergenceStop+",\"run_growing\":"+negativeGrowing+",\"stop_reason\":\"convergence\",\"stage_order\":"+negativeTrace+"}]}");
        System.out.println();
    }
}
