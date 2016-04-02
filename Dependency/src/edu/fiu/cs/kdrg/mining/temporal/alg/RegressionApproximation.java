package edu.fiu.cs.kdrg.mining.temporal.alg;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;

public interface RegressionApproximation {

	public void train(SparseVector vectors);

	public SparseVector getPosteriorMean(int dim);

	public SparseVector getPosteriorVariance(int dim);

}
