package edu.fiu.cs.kdrg.mining.temporal.alg;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordWriter;
import edu.fiu.cs.kdrg.mining.temporal.util.Utility;

public class BayesianLasso implements OnlineRegression {

	/**
	 * The number of dimensions.
	 */
	private int dimension = 0;

	/**
	 * The number of particles.
	 */
	private int numOfParticles = 0;

	/**
	 * The time lag for time frame.
	 */
	private int lag = 1;

	/**
	 * The time series within a time lag.
	 */
	private TimeFrame tf = null;

	/**
	 * The parameter for L1 regularization about the constant part of
	 * coefficient.
	 */
	private double lambda = 1.0;

	/**
	 * The coefficients for the linear model. Since it's a vector auto
	 * regression, coefficients contains dimension*numOfParticles SimpleMatrix.
	 * The dimension of each SimpleMatrix is (dimension X 1). Each SimpleMatrix
	 * refers c_t1.
	 * 
	 * y_t = x_t * c_t1 + noise.
	 */
	private SimpleMatrix[][] coefficients = null;

	/**
	 * It is the variance of noise, where noise ~ N(0,sigma^2). It contains
	 * dimension*numOfParticles variances.
	 * 
	 * y_t = x_t * c_t1 + x_t*state_t* c_t2 + noise.
	 */
	private double[][] sigmaSquares = null;

	/**
	 * IG(a, b): a prior distribution for sigma^2. All share a single prior a.
	 */
	private double aPriorForSigmaSquare = 1.0;

	/**
	 * IG(a, b): b prior distribution for sigma^2. It contains
	 * dimension*numOfParticles prior b.
	 */
	private double[][] bPriorForSigmaSquares = null;

	/**
	 * coefficient ~ N(B_t, sigma^2 * M_t): M_t is the covariates. It contains
	 * dimension*numOfParticles SimpleMatrix. The dimension of covariate is :
	 * (2* dimension) X (2*dimension).
	 */
	private SimpleMatrix[][] covariatesForCoefficient = null;

	/**
	 * coefficient ~ N(B_t, sigma^2 * M_t): B_t is the means. It contains
	 * dimension*numOfParticles SimpleMatrix. The dimension of each SimpleMatrix
	 * is dimension*2.
	 */
	private SimpleMatrix[][] meansForCoefficient = null;

	/**
	 * The random variable for coefficients follows exponential distribution. It
	 * contains dimension*numOfParticles SimpleMatrix. The dimension of each
	 * SimpleMatrix is (dimension X dimension).
	 */
	private SimpleMatrix[][] randForCoefficient = null;

	/**
	 * It denotes dimension*numOfParticles weights for all particles of all
	 * dimensions.
	 */
	private double[][] particleWeights = null;

	/**
	 * Generate sample from exponential distribution with parameter lambda1.
	 */
	private ExponentialDistribution exponent = null;

	/**
	 * The result writer.
	 */
	private RecordWriter rw;

	public BayesianLasso() {
	}

	public BayesianLasso(int dim, int numOfParticles) {
		this(dim, numOfParticles, 1, 1.0);

	}

	public BayesianLasso(int dim, int numOfParticles, int lag, double lambda) {
		this.dimension = dim;
		this.numOfParticles = numOfParticles;
		this.lag = lag;
		this.lambda = lambda;

		exponent = new ExponentialDistribution(2.0 / (lambda * lambda));

		this.sigmaSquares = new double[dimension][numOfParticles];
		this.aPriorForSigmaSquare = 1.0;
		this.bPriorForSigmaSquares = new double[dimension][numOfParticles];

		this.coefficients = new SimpleMatrix[dimension][numOfParticles];
		this.meansForCoefficient = new SimpleMatrix[dimension][numOfParticles];
		this.covariatesForCoefficient = new SimpleMatrix[dimension][numOfParticles];

		this.particleWeights = new double[dimension][numOfParticles];
		this.randForCoefficient = new SimpleMatrix[dimension][numOfParticles];

		for (int d = 0; d < dimension; d++) {
			for (int p = 0; p < numOfParticles; p++) {
				bPriorForSigmaSquares[d][p] = 1.0;
				sigmaSquares[d][p] = Utility.sampleInverseGammaDistribution(this.aPriorForSigmaSquare,
						this.bPriorForSigmaSquares[d][p]);

				meansForCoefficient[d][p] = new SimpleMatrix(getNumOfFeatures(), 1);
				covariatesForCoefficient[d][p] = SimpleMatrix.identity(getNumOfFeatures());

				particleWeights[d][p] = 1.0 / numOfParticles;
				randForCoefficient[d][p] = Utility.sampleExpDistribution(getNumOfFeatures(), exponent);

				// sample coefficients
				SimpleMatrix var = (randForCoefficient[d][p].elementPower(2)
						.mult(covariatesForCoefficient[d][p].scale(sigmaSquares[d][p])));
				coefficients[d][p] = Utility.sampleMultivariateGaussion(meansForCoefficient[d][p], var);
			}
		}

	}

	/**
	 * When a new instance comes, the weights for each particles are
	 * recalculated.
	 * 
	 * @param instance
	 *            the new coming instance
	 */
	protected void computeWeights(SimpleMatrix instance) {
		SimpleMatrix x = tf.getVector();
		double defaultWeight = 1.0 / numOfParticles;
		TDistribution student = Utility.createTDistribution(2 * aPriorForSigmaSquare);
		for (int d = 0; d < dimension; d++) {
			double y = instance.get(d);
			double mean = 0.0;
			double var = 0.0;
			double total = 0.0;
			for (int p = 0; p < numOfParticles; p++) {
				SimpleMatrix weights = coefficients[d][p];
				mean = x.dot(weights);

				SimpleMatrix xtemp = randForCoefficient[d][p].mult(x);
				var = aPriorForSigmaSquare / bPriorForSigmaSquares[d][p]
						* (xtemp.dot(covariatesForCoefficient[d][p].mult(xtemp)) + 1.0);

				double w = student.density((y - mean) / Math.sqrt(var));
				total += w;
				particleWeights[d][p] = w;
			}

			// normalize
			for (int indexP = 0; indexP < numOfParticles; indexP++) {
				if (total == 0.0) {
					particleWeights[d][indexP] = defaultWeight;
				} else {
					particleWeights[d][indexP] = particleWeights[d][indexP] / total;
				}
			}
		}
	}

	/**
	 * This algorithm sample the element in place and the time cost if O(2n).
	 */
	protected void resample() {
		for (int d = 0; d < dimension; d++) {
			int[] frequency = Utility.sample(numOfParticles, particleWeights[d], numOfParticles);
			int src = 0;
			int dst = 0;
			int len = frequency.length;
			int index = 0;
			while (index < numOfParticles) {
				// find the elements with frequency > 0
				while (src < len && frequency[src] == 0) {
					src++;
				}
				index += frequency[src];

				// move frequency -1 to other empty positions
				int i = 1;
				while (i < frequency[src]) {
					// find the next empty position
					while (dst < len && frequency[dst] != 0) {
						dst++;
					}
					// copy data[src] to data[dst]
					copyParticleFromSrcToDst(d, src, dst);
					// the dst postion is filled
					dst++;
					// the number of elements at the src position filled.
					i++;
				}
				// the element at src is done.
				src++;
			}
		}
	}

	protected void copyParticleFromSrcToDst(int dim, int src, int dst) {
		coefficients[dim][dst] = coefficients[dim][src].copy();
		sigmaSquares[dim][dst] = sigmaSquares[dim][dst];
		bPriorForSigmaSquares[dim][dst] = bPriorForSigmaSquares[dim][dst];
		covariatesForCoefficient[dim][dst] = covariatesForCoefficient[dim][src].copy();
		meansForCoefficient[dim][dst] = meansForCoefficient[dim][src].copy();
		randForCoefficient[dim][dst] = randForCoefficient[dim][src].copy();
		particleWeights[dim][dst] = particleWeights[dim][src];
	}

	/**
	 * Update the statistics of parameter.
	 * 
	 * @param instance
	 *            the new coming instance
	 */
	protected void updateParameters(SimpleMatrix instance) {
		SimpleMatrix x = tf.getVector();
		// update the parameter a of sigma^2 ~ ig(a, b)
		aPriorForSigmaSquare += 0.5;

		for (int d = 0; d < dimension; d++) {
			double y = instance.get(d);
			for (int p = 0; p < numOfParticles; p++) {
				// update the variance of coefficient
				SimpleMatrix obs = x;
				SimpleMatrix old_covariatesForCoefficient = covariatesForCoefficient[d][p].copy();

				covariatesForCoefficient[d][p] = old_covariatesForCoefficient.invert()
						.plus(randForCoefficient[d][p].mult(obs.mult(obs.transpose())).mult(randForCoefficient[d][p]))
						.invert();

				// update the mean of coefficient
				SimpleMatrix temp1 = randForCoefficient[d][p].mult(covariatesForCoefficient[d][p])
						.mult(randForCoefficient[d][p]);
				SimpleMatrix temp2 = covariatesForCoefficient[d][p].mult(old_covariatesForCoefficient.invert());
				temp2 = randForCoefficient[d][p].mult(temp2).mult(randForCoefficient[d][p].invert());
				SimpleMatrix old_meansForCoefficient = meansForCoefficient[d][p].copy();
				meansForCoefficient[d][p] = temp1.mult(obs).scale(y).plus(temp2.mult(old_meansForCoefficient));

				// update random variable for coefficients
				SimpleMatrix old_randForCoefficient = randForCoefficient[d][p].copy();
				randForCoefficient[d][p] = Utility.sampleExpDistribution(getNumOfFeatures(), exponent);

				// update the parameter b of sigma^2 ~ ig(a, b)
				SimpleMatrix bTemp = old_randForCoefficient
						.mult(old_covariatesForCoefficient.mult(old_randForCoefficient));
				double bTemp1 = 1.0 + obs.dot(bTemp.mult(obs));
				double val = y - obs.dot(meansForCoefficient[d][p]);
				val = val * val / bTemp1;
				bPriorForSigmaSquares[d][p] += 0.5 * val;

				// sample sigma^2
				sigmaSquares[d][p] = Utility.sampleInverseGammaDistribution(aPriorForSigmaSquare,
						bPriorForSigmaSquares[d][p]);

				// sample coefficients
				SimpleMatrix var = (randForCoefficient[d][p]
						.mult(covariatesForCoefficient[d][p].scale(sigmaSquares[d][p])));
				coefficients[d][p] = Utility.sampleMultivariateGaussion(meansForCoefficient[d][p], var);
			}
		}
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

		int numFeatures = this.getNumOfFeatures();
		SimpleMatrix ret = new SimpleMatrix(numFeatures, 1);
		for (int f = 0; f < numFeatures; f++) {
			double temp = 0.0;
			for (int p = 0; p < numOfParticles; p++) {
				temp += (meansForCoefficient[dim][p].get(f));
			}
			ret.set(f, temp / numOfParticles);
		}

		return ret;
	}

	@Override
	public SimpleMatrix getPosteriorVariance(int dim) {
		// TODO Auto-generated method stub
		int numFeatures = this.getNumOfFeatures();
		SimpleMatrix ret = new SimpleMatrix(numFeatures, 1);
		return ret;
	}

	public SimpleMatrix[] getMeans() {
		SimpleMatrix[] ret = new SimpleMatrix[dimension];
		for (int i = 0; i < dimension; i++) {
			ret[i] = getPosteriorMean(i);
		}
		return ret;
	}

	@Override
	public void onlineTrain(SimpleMatrix instance) {
		computeWeights(instance);

		resample();

		updateParameters(instance);

		if (rw != null)
			rw.outLine(this.getMeans());
		tf.shift(instance);

	}

	private int getNumOfFeatures() {
		return dimension * lag;
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
