/** Pinned AlgorithmMultipleToads post-membership call sequence excerpt. */
public final class ToadsDefaultRegistrationHandoffReferenceHarness {
    private static final String WRAPPER_SHA = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d";
    private static final String ATLAS_SHA = "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db";
    public static void main(String[] args) {
        // Source excerpt in AlgorithmMultipleToads: after the membership
        // energy update, ratio is computed, optional registration runs, then
        // convergence/iteration stops are checked, then growth/output dispatch.
        String[] callOrder = {
            "computeExactMemberships", "energyRatio", "registerShapes",
            "convergenceCheck", "iterationLimitCheck", "propagateGrowing",
            "classificationOrRelationRefresh"
        };
        StringBuilder order = new StringBuilder("[");
        for (int i=0; i<callOrder.length; i++) {
            if (i>0) order.append(',');
            order.append('"').append(callOrder[i]).append('"');
        }
        System.out.print("{\"algorithm_source_sha256\":\""+WRAPPER_SHA
            +"\",\"atlas_source_sha256\":\""+ATLAS_SHA
            +"\",\"call_order\":"+order+
            "],\"align_shapes_default\":true,\"registration_transform_type\":\"rigid\","
            +"\"registration_dispatch\":\"SINGLE/registerAllShapes\","
            +"\"registration_resets\":[\"lambda=INIT_LAMBDA\",\"Nturn=0\",\"itPlus=0\"],"
            +"\"registration_stop\":[\"n < itMax\",\"diff > minEdiff\",\"lambda > minLambda\"],"
            +"\"registration_finalize\":[\"computeRotation(transform)\",\"computeTransformedShapeBoundingBox()\"]}");
    }
}
