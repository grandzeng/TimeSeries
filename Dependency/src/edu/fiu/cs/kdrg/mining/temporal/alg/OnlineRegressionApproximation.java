package edu.fiu.cs.kdrg.mining.temporal.alg;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;

public interface OnlineRegressionApproximation extends RegressionApproximation {
	public void onlineTrain(SparseVector vector);
}
