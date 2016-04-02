package edu.fiu.cs.kdrg.mining.temporal.alg;

import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;
import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordWriter;

/**
 * It can incremental train a vector auto-regression model from a Bayesian
 * perspective. This algorithm is approximation version which only considers the
 * diagonal elements.
 * 
 * @author Chunqiu Zeng
 * @date Feb.19, 2016
 *
 */
public class BayesianLinearRegressionApproximation implements OnlineRegression {

	private int lag;

	private int dimension;

	private SparseVector[] meanWeights;

	private SparseVector precisionWeights;

	private double iga;

	private double[] igbs;

	private TimeFrame tf;

	private RecordWriter rw;

	public BayesianLinearRegressionApproximation(int lag, int dimension) {
		this.lag = lag;
		this.dimension = dimension;
		init();
	}

	private int numFeatures() {
		return lag * dimension;
	}

	private void init() {
		iga = 1.0;
		igbs = new double[dimension];
		tf = new TimeFrame(lag, dimension);

		meanWeights = new SparseVector[dimension];
		for (int i = 0; i < dimension; i++) {
			meanWeights[i] = new SparseVector(numFeatures(), 1);
			igbs[i] = 1.0;
		}

		precisionWeights = new SparseVector(numFeatures(), 1.0);
	}

	@Override
	public void train(SimpleMatrix vectors) {
		for (int i = 0; i < vectors.numCols(); i++) {
			SimpleMatrix vec = vectors.extractVector(false, i);
			onlineTrain(vec);
		}

	}

	@Override
	public SimpleMatrix getPosteriorMean(int dim) {
		SimpleMatrix ret = new SimpleMatrix(numFeatures(), 1);
		for (int i = 0; i < meanWeights[dim].getDimension(); i++) {
			ret.set(i, meanWeights[dim].get(i));
		}
		return ret;
	}

	@Override
	public SimpleMatrix getPosteriorVariance(int dim) {
		double v = 2 * iga;
		SparseVector sigma = precisionWeights.copyNew().inverse().mul(igbs[dim] / iga).mul(v / (v - 2));
		SimpleMatrix ret = new SimpleMatrix(numFeatures(), numFeatures());
		for (int i = 0; i < sigma.getDimension(); i++) {
			ret.set(i, i, sigma.get(i));
		}
		return ret;
	}

	@Override
	public void onlineTrain(SimpleMatrix vector) {
		SparseVector x = tf.getSparseVector();
		SparseVector temp = x.elementwiseMultiply(x);
		// SimpleMatrix temp = x.mult(x.transpose());
		SparseVector newPrecisionWeights = precisionWeights.copyNew().add(temp);
		SparseVector newMeanWeights;
		iga += 0.5;
		double y;
		SparseVector xy;
		for (int i = 0; i < dimension; i++) {
			y = vector.get(i);
			xy = x.copyNew().mul(y);
			SparseVector temp1 = precisionWeights.elementwiseMultiply(meanWeights[i]).add(xy);
			newMeanWeights = newPrecisionWeights.copyNew().inverse().elementwiseMultiply(temp1);
			igbs[i] += 0.5 * (meanWeights[i].elementwiseMultiply(precisionWeights).innerProduct(meanWeights[i]) + y * y
					- newMeanWeights.elementwiseMultiply(newPrecisionWeights).innerProduct(newMeanWeights));
			meanWeights[i] = newMeanWeights;
		}

		precisionWeights = newPrecisionWeights;
		if (rw != null)
			rw.outLine(meanWeights);
		tf.shift(vector);

	}
	
	public void setTimeFrame(TimeFrame tf) {
		this.tf = tf.copy();
	}

	public int getLag() {
		return lag;
	}

	public int getDimension() {
		return dimension;
	}

	public void setOutput(RecordWriter rw) {
		this.rw = rw;
	}

}
