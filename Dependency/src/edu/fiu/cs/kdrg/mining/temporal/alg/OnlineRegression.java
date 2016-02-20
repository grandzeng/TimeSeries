package edu.fiu.cs.kdrg.mining.temporal.alg;

import org.ejml.simple.SimpleMatrix;

public interface OnlineRegression extends Regression {
	public void onlineTrain(SimpleMatrix vector);
}
