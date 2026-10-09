package edu.jhu.ece.iacl.jist.pipeline;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
/** Headless lifecycle stub for selected JIST algorithm source. */
public abstract class ProcessingAlgorithm {
    protected abstract void createInputParameters(ParamCollection inputParams);
    protected abstract void createOutputParameters(ParamCollection outputParams);
    protected abstract void execute(CalculationMonitor monitor);
}
