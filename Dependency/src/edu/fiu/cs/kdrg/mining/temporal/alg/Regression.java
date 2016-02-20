package edu.fiu.cs.kdrg.mining.temporal.alg;

import org.ejml.simple.SimpleMatrix;

public interface Regression {

	public void train(SimpleMatrix vectors);

	public SimpleMatrix getPosteriorMean(int dim);
	
	public SimpleMatrix getPosteriorVariance(int dim);

}
