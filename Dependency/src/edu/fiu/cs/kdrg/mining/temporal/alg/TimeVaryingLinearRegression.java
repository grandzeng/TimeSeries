package edu.fiu.cs.kdrg.mining.temporal.alg;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordWriter;
import edu.fiu.cs.kdrg.mining.temporal.util.Utility;

public class TimeVaryingLinearRegression implements OnlineRegression {

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
	private double lambda1 = 1.0;

	/**
	 * The parameter for L1 regularization about the varying part of
	 * coefficient.
	 */
	private double lambda2 = 1.0;

	/**
	 * The coefficients for the linear model. Since it's a vector auto
	 * regression, coefficients contains dimension*numOfParticles SimpleMatrix.
	 * The dimension of each SimpleMatrix is (2*dimension X 1). Each
	 * SimpleMatrix refers (c_t1,c_t2).
	 * 
	 * y_t = x_t * c_t1 + x_t*state_t* c_t2 + noise.
	 */
	private SimpleMatrix[][] coefficients = null;

	/**
	 * The states change over time. The states contains dimension*numOfParticles
	 * SimpleMatrix. The dimension of each SimpleMatrix is (dimension X 1). Each
	 * SimpleMatrix refers state_t.
	 * 
	 * y_t = x_t * c_t1 + x_t*state_t* c_t2 + noise.
	 */
	private SimpleMatrix[][] states = null;

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
	 * SimpleMatrix is (2*dimension X 2*dimension).
	 */
	private SimpleMatrix[][] randForCoefficient = null;

	/* Kalman filtering sufficient statistics of state */

	/**
	 * the means for state at time t. i.e., m_t
	 */
	private SimpleMatrix[][] meansForState = null;

	/**
	 * The variance for state at time t. i.e., C_t (dim X dim). It contains
	 * dimension*numOfParticles SimpleMatrix. The dimension of each SimpleMatrix
	 * is (dimension x dimension).
	 */
	private SimpleMatrix[][] variancesForState = null;

	/**
	 * The variance for predictive error. i.e., Q_t. It contains
	 * dimension*numOfParticles values.
	 */
	private double[][] variancesForPredictiveError = null;

	/**
	 * The kalman gain. i.e., A_t (kx1). It contains dimension*numOfParticles
	 * SimpleMatrix. The dimension of each SimpleMatrix is (dimension x 1).
	 */
	private SimpleMatrix[][] kalmanGain = null;

	/**
	 * It denotes dimension*numOfParticles weights for all particles of all
	 * dimensions.
	 */
	private double[][] particleWeights = null;

	/**
	 * Generate sample from exponential distribution with parameter lambda1.
	 */
	private ExponentialDistribution exp1 = null;

	/**
	 * Generate sample from exponential distribution with parameter lambda2.
	 */
	private ExponentialDistribution exp2 = null;

	/**
	 * The result writer.
	 */
	private RecordWriter rw;

	/**
	 * Default constructor.
	 */
	public TimeVaryingLinearRegression() {
	}

	/**
	 * Constructor
	 * 
	 * @param dim
	 * @param numOfParticles
	 */
	public TimeVaryingLinearRegression(int dim, int numOfParticles) {
		this(dim, numOfParticles, 1, 1.0, 1.0);
	}

	/**
	 * Constructor
	 * 
	 * @param dim
	 * @param numOfParticles
	 * @param lambda1
	 * @param lambda2
	 */
	public TimeVaryingLinearRegression(int dim, int numOfParticles, int lag, double lambda1, double lambda2) {
		this.dimension = dim;
		this.lag = lag;
		this.numOfParticles = numOfParticles;
		this.lambda1 = lambda1;
		this.lambda2 = lambda2;

		exp1 = new ExponentialDistribution(2.0 / (lambda1 * lambda1));
		exp2 = new ExponentialDistribution(2.0 / (lambda2 * lambda2));

		this.sigmaSquares = new double[dimension][numOfParticles];
		this.aPriorForSigmaSquare = 1.0;
		this.bPriorForSigmaSquares = new double[dimension][numOfParticles];

		this.coefficients = new SimpleMatrix[dimension][numOfParticles];
		this.meansForCoefficient = new SimpleMatrix[dimension][numOfParticles];
		this.covariatesForCoefficient = new SimpleMatrix[dimension][numOfParticles];

		this.states = new SimpleMatrix[dimension][numOfParticles];
		this.meansForState = new SimpleMatrix[dimension][numOfParticles];
		this.variancesForState = new SimpleMatrix[dimension][numOfParticles];

		this.kalmanGain = new SimpleMatrix[dimension][numOfParticles];
		this.variancesForPredictiveError = new double[dimension][numOfParticles];

		this.particleWeights = new double[dimension][numOfParticles];
		this.randForCoefficient = new SimpleMatrix[dimension][numOfParticles];
		
		tf = new TimeFrame(lag, dimension);

		for (int d = 0; d < dimension; d++) {
			for (int p = 0; p < numOfParticles; p++) {
				bPriorForSigmaSquares[d][p] = 1.0;
				sigmaSquares[d][p] = Utility.sampleInverseGammaDistribution(this.aPriorForSigmaSquare,
						this.bPriorForSigmaSquares[d][p]);

				meansForCoefficient[d][p] = new SimpleMatrix(2 * getNumOfFeatures(), 1);
				covariatesForCoefficient[d][p] = SimpleMatrix.identity(2 * getNumOfFeatures());

				meansForState[d][p] = new SimpleMatrix(getNumOfFeatures(), 1);
				variancesForState[d][p] = new SimpleMatrix(getNumOfFeatures(), getNumOfFeatures());

				kalmanGain[d][p] = new SimpleMatrix(getNumOfFeatures(), 1);
				variancesForPredictiveError[d][p] = 0.0;

				particleWeights[d][p] = 1.0 / numOfParticles;
				randForCoefficient[d][p] = Utility.sampleExpDistribution(getNumOfFeatures(), exp1, exp2);

				// System.out.println(randForCoeffi);

				SimpleMatrix var = (covariatesForCoefficient[d][p].mult(randForCoefficient[d][p]))
						.scale(sigmaSquares[d][p]);
				coefficients[d][p] = Utility.sampleMultivariateGaussion(meansForCoefficient[d][p], var);

				states[d][p] = new SimpleMatrix(getNumOfFeatures(), 1);// Utility.sampleMultivariateGaussion(meansForState[d][p],
																		// variancesForState[d][p]);
			}
		}

	}

	protected SimpleMatrix combineEfficientAndMeansOfState(int indexParticle) {
		int numOfFeatures = getNumOfFeatures();
		SimpleMatrix ret = new SimpleMatrix(numOfFeatures, 1);
		double temp = 0.0;
		for (int d = 0; d < dimension; d++)
			for (int i = 0; i < numOfFeatures; i++) {
				temp = coefficients[d][indexParticle].get(i) + coefficients[d][indexParticle].get(i + numOfFeatures)
						* meansForState[d][indexParticle].get(i);
				ret.set(i, temp);
			}
		return ret;
	}

	protected SimpleMatrix combineEfficientAndState(int indexParticle) {
		int numOfFeatures = getNumOfFeatures();
		SimpleMatrix ret = new SimpleMatrix(numOfFeatures, 1);
		double temp = 0.0;
		for (int d = 0; d < dimension; d++)
			for (int i = 0; i < numOfFeatures; i++) {
				temp = coefficients[d][indexParticle].get(i)
						+ coefficients[d][indexParticle].get(i + numOfFeatures) * states[d][indexParticle].get(i);
				ret.set(i, temp);
			}
		return ret;
	}

	protected SimpleMatrix getStateCoefficient(int d, int indexP) {
		int numFeatures = getNumOfFeatures();
		SimpleMatrix stateCoefficient = new SimpleMatrix(numFeatures, 1);
		double value = 0.0;
		for (int i = 0; i < numFeatures; i++) {
			value = coefficients[d][indexP].get(i + numFeatures);
			stateCoefficient.set(i, value);
		}
		return stateCoefficient;
	}

	protected SimpleMatrix getObservationCoefficient(int d, int indexP) {
		int numFeatures = getNumOfFeatures();
		SimpleMatrix observationCoefficient = new SimpleMatrix(numFeatures, 1);
		double value = 0.0;
		for (int i = 0; i < numFeatures; i++) {
			value = coefficients[d][indexP].get(i);
			observationCoefficient.set(i, value);
		}
		return observationCoefficient;
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
		SimpleMatrix identity = SimpleMatrix.identity(getNumOfFeatures());
		double defaultWeight = 1.0 / numOfParticles;
		for (int d = 0; d < dimension; d++) {
			double y = instance.get(d);
			double mean = 0.0;
			double var = 0.0;
			double total = 0.0;
			for (int p = 0; p < numOfParticles; p++) {
				SimpleMatrix weights = combineEfficientAndMeansOfState(p);
				mean = x.dot(weights);
				SimpleMatrix I_plus_C = identity.plus(variancesForState[d][p]);
				SimpleMatrix temp = x.elementMult(getStateCoefficient(d, p));
//				System.out.println(temp.numRows() +"X"+temp.numCols());
//				System.out.println(I_plus_C.numRows() + "X" + I_plus_C.numCols());
				var = sigmaSquares[d][p] + temp.transpose().mult(I_plus_C).dot(temp);

				double w = Utility.standardNormDensity((y - mean) / Math.sqrt(var));
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
		states[dim][dst] = states[dim][src].copy();
		sigmaSquares[dim][dst] = sigmaSquares[dim][dst];
		bPriorForSigmaSquares[dim][dst] = bPriorForSigmaSquares[dim][dst];
		covariatesForCoefficient[dim][dst] = covariatesForCoefficient[dim][src].copy();
		meansForCoefficient[dim][dst] = meansForCoefficient[dim][src].copy();
		randForCoefficient[dim][dst] = randForCoefficient[dim][src].copy();
		meansForState[dim][dst] = meansForState[dim][src].copy();
		variancesForState[dim][dst] = variancesForState[dim][src].copy();
		variancesForPredictiveError[dim][dst] = variancesForPredictiveError[dim][src];
		kalmanGain[dim][dst] = kalmanGain[dim][src].copy();
		particleWeights[dim][dst] = particleWeights[dim][src];
	}

	/**
	 * Update the statistics of state.
	 * 
	 * @param instance
	 *            the new coming instance
	 */
	protected void kalmanUpdateStates(SimpleMatrix instance) {

		SimpleMatrix x = tf.getVector();
		SimpleMatrix identity = SimpleMatrix.identity(getNumOfFeatures());
		for (int d = 0; d < dimension; d++) {
			double y = instance.get(d);
			for (int p = 0; p < numOfParticles; p++) {
				// update Q_t
				SimpleMatrix I_plus_C = identity.plus(variancesForState[d][p]);
				SimpleMatrix stateCoef = getStateCoefficient(d, p);
				SimpleMatrix temp = x.elementMult(stateCoef);
				SimpleMatrix temp1 = I_plus_C.mult(temp);
				variancesForPredictiveError[d][p] = sigmaSquares[d][p] + temp.dot(temp1);
				// update A_t
				kalmanGain[d][p] = temp1.divide(variancesForPredictiveError[d][p]);
				// update M_t
				SimpleMatrix weights = combineEfficientAndMeansOfState(p);
				double b = y - x.dot(weights);
				meansForState[d][p] = meansForState[d][p].plus(kalmanGain[d][p].scale(b));

				// update C_t
				SimpleMatrix help = kalmanGain[d][p]
						.mult(kalmanGain[d][p].transpose().scale(variancesForPredictiveError[d][p]));
				variancesForState[d][p] = I_plus_C.minus(help);

				// sample states
				states[d][p] = Utility.sampleMultivariateGaussion(meansForState[d][p], variancesForState[d][p]);
			}
		}
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
				SimpleMatrix obs = x.combine(SimpleMatrix.END, 0, x.elementMult(states[d][p]));
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
				randForCoefficient[d][p] = Utility.sampleExpDistribution(getNumOfFeatures(), exp1, exp2);

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
				SimpleMatrix var = (randForCoefficient[d][p].elementPower(2)
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
				temp += (meansForCoefficient[dim][p].get(f)
						+ meansForCoefficient[dim][p].get(f + numFeatures) * meansForState[dim][p].get(f));
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

		kalmanUpdateStates(instance);

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
