/*
 * Test-only reference excerpt for the final exact-registration coefficient
 * assembly in VandyCRUISE MultipleToadDeformableAtlas. It intentionally keeps
 * the source call order, including passing the already-normalized Hessian to
 * registrationMeasureGradient. This documents observed source behavior, not
 * a recommendation for a corrected optimizer.
 */
public final class ToadsRegistrationCoefficientReferenceHarness {
    private static float registrationMeasureHessian(float hess, float nhess,
            float norm, float scaling) {
        return hess / norm;
    }

    private static float registrationMeasureGradient(float grad, float ngrad,
            float norm, float scaling) {
        return grad / norm;
    }

    public static void main(String[] args) {
        float[] gradient = {0.75f, -0.125f, 3.5f};
        float[] hessian = {2.0f, 0.5f, 8.25f};
        float[] normGradient = {12.0f, 4.0f, -1.0f};
        float[] normHessian = {9.0f, 8.0f, 7.0f};
        float norm = 2.5f;
        float scaling = 60.0f;
        StringBuilder result = new StringBuilder("{\"gradient\":[");
        for (int i = 0; i < gradient.length; i++) {
            hessian[i] = registrationMeasureHessian(hessian[i], normHessian[i], norm, scaling);
            gradient[i] = registrationMeasureGradient(hessian[i], normHessian[i], norm, scaling);
            if (i > 0) result.append(',');
            result.append(Float.toString(gradient[i]));
        }
        result.append("],\"hessian\":[");
        for (int i = 0; i < hessian.length; i++) {
            if (i > 0) result.append(',');
            result.append(Float.toString(hessian[i]));
        }
        System.out.println(result.append("]}"));
    }
}
