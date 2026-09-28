package edu.jhu.ece.iacl.algorithms.gvf;

import edu.jhu.ece.iacl.algorithms.VersionUtil;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
/**
 * Gradient Vector Flow Calculation using either generalized or original method
 * @author Blake Lucas
 *
 */
public abstract class GradVecFlow extends AbstractCalculation{
	public static String getVersion() {
		return VersionUtil.parseRevisionNumber("$Revision: 1.1 $");
	}



	public GradVecFlow(AbstractCalculation parent){
		super(parent);
		setLabel("Gradient Vector Flow");
	}
	public GradVecFlow(){
		super();
		setLabel("Gradient Vector Flow");
	}
	/**
	 * Create Gradient Vector Flow (GVF) solver
	 * @param method {"Generalized","Original"}
	 * @return GVF solver
	 */
	public static GradVecFlow create(String method){
		if(method.equals("Original"))return new GradVecFlowOriginal();
		else if(method.equals("Generalized"))return new GradVecFlowGeneralized();
		else return null;
	}

	/**
	 * Solve GVF using matrix 3d representation of volume
	 * 
	 * @param vol
	 *            matrix 3d
	 * @param lambda
	 *            lambda
	 * @param normalize
	 *            normalize if true
	 * @param iter
	 *            maximum iterations
	 * @return vector matrix3d
	 */
	public abstract ImageData solve(ImageData vol,float lambda,boolean normalize,int iter);
}
