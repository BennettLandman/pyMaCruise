import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Selected wrapper call-order and atlas registration-selection excerpt. */
public final class ToadsLesionAtlasRegistrationHandoffReferenceHarness {
    private static final String WRAPPER_SHA="b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String ATLAS_SHA="e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2";

    private static List<String> descriptorTypes(Path path) throws Exception {
        List<String> lines=Files.readAllLines(path);
        int header=-1;
        for(int i=0;i<lines.size();i++) if(lines.get(i).startsWith("Structures:")){header=i;break;}
        if(header<0) throw new IllegalArgumentException("descriptor has no Structures header");
        int count=Integer.parseInt(lines.get(header).split(":",2)[1].trim());
        List<String> types=new ArrayList<>();
        for(int k=0;k<count;k++){
            String[] fields=lines.get(header+1+k).split("\\t");
            String tag=fields[2];
            if(tag.equals("mask")) types.add("MASK");
            else if(tag.equals("optimized_csf")) types.add("OPTIMIZED_CSF");
            else if(tag.equals("lesion")) types.add("LESION");
            else types.add("OBJECT");
        }
        return types;
    }

    public static void main(String[] args) throws Exception {
        List<String> types=descriptorTypes(Path.of(args[0]));
        List<Integer> registered=new ArrayList<>();
        for(int k=0;k<types.size();k++) if(!types.get(k).equals("MASK")) registered.add(k);
        System.out.println("{\"wrapper_source_sha256\":\""+WRAPPER_SHA+"\","
            +"\"atlas_source_sha256\":\""+ATLAS_SHA+"\","
            +"\"class_types\":"+quoted(types)+","
            +"\"registered_class_indices\":"+integers(registered)+","
            +"\"call_order\":[\"computeExactMemberships\",\"energyRatio\",\"registerShapes\",\"convergenceCheck\",\"iterationLimitCheck\"],"
            +"\"registration_dispatch\":\"SINGLE/registerAllShapes\","
            +"\"registration_resets\":[\"lambda=INIT_LAMBDA\",\"Nturn=0\",\"itPlus=0\"],"
            +"\"registration_stop\":[\"n < itMax\",\"diff > minEdiff\",\"lambda > minLambda\"],"
            +"\"registration_finalize\":[\"computeRotation(transform)\",\"computeTransformedShapeBoundingBox()\"]}");
    }

    private static String quoted(List<String> values){
        List<String> out=new ArrayList<>();for(String v:values)out.add("\""+v+"\"");return out.toString();
    }
    private static String integers(List<Integer> values){return values.toString();}
}
