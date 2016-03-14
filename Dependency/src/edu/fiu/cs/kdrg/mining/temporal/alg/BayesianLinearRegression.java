package edu.fiu.cs.kdrg.mining.temporal.alg;

import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordWriter;

/**
 * 
 * It can incremental train a vector auto-regression model from a Bayesian
 * perspective.
 * 
 * @author Chunqiu Zeng
 * @date Feb.19, 2016
 *
 */
public class BayesianLinearRegression implements OnlineRegression {

	private int lag;

	private int dimension;

	private SimpleMatrix[] meanWeights;

	private SimpleMatrix precisionWeights;

	private double iga;

	private double[] igbs;

	private TimeFrame tf;

	private RecordWriter rw;

	public BayesianLinearRegression(int lag, int dimension) {
		this.lag = lag;
		this.dimension = dimension;
		init();
	}

	private void init() {
		meanWeights = new SimpleMatrix[dimension];
		iga = 1.0;
		igbs = new double[dimension];
		tf = new TimeFrame(lag, dimension);

		meanWeights = new SimpleMatrix[dimension];
		for (int i = 0; i < dimension; i++) {
			meanWeights[i] = new SimpleMatrix(lag * dimension, 1);
			igbs[i] = 1.0;
		}
		precisionWeights = SimpleMatrix.identity(lag * dimension);
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
		return meanWeights[dim];
	}

	@Override
	public SimpleMatrix getPosteriorVariance(int dim) {
		// return precisionWeights.invert();
		double v = 2 * iga;
		SimpleMatrix sigma = precisionWeights.invert().scale(igbs[dim] / iga);
		return sigma.scale(v / (v - 2));
	}

	@Override
	public void onlineTrain(SimpleMatrix vector) {
		SimpleMatrix x = tf.getVector();
		SimpleMatrix temp = x.mult(x.transpose());
		SimpleMatrix newPrecisionWeights = precisionWeights.plus(temp);
		SimpleMatrix newMeanWeights;
		iga += 0.5;
		double y;
		SimpleMatrix xy;
		for (int i = 0; i < dimension; i++) {
			y = vector.get(i);
			xy = x.scale(y);
			SimpleMatrix temp1 = precisionWeights.mult(meanWeights[i]).plus(xy);
			newMeanWeights = newPrecisionWeights.invert().mult(temp1);
			igbs[i] += 0.5 * (meanWeights[i].transpose().mult(precisionWeights).dot(meanWeights[i]) + y * y
					- newMeanWeights.transpose().mult(newPrecisionWeights).dot(newMeanWeights));
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
