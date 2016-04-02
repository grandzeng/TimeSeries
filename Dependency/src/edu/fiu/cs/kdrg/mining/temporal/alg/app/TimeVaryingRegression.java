package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.io.IOException;
import java.io.Writer;
import java.util.Collection;
import java.util.Random;

import org.apache.commons.math3.distribution.ExponentialDistribution;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;
import edu.fiu.cs.kdrg.mining.temporal.tv.util.Utility;
import edu.fiu.cs.kdrg.mining.temporal.tv.util.MathFunctions;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TLongList;

public abstract class TimeVaryingRegression implements OnlineBayesianRegression {

	protected int numOfParticles = 1;

	protected int dimension = 0;

	protected long impression = 0;

	protected long click = 0;

	/**
	 * The variance of noise: i.e., sigma^2
	 */
	protected double[] squareOfSigma = null;

	/**
	 * IG(a, b): a prior distribution for sigma^2
	 */
	protected double[] aPriorOfSquareSigma = null;

	/**
	 * IG(a, b): b prior distribution for sigma^2
	 */
	protected double[] bPriorOfSquareSigma = null;

	/**
	 * The coefficients for the linear model. The number of coefficients is
	 * twice of the dimension.
	 * 
	 * y_t = x_t * c_t1 + x_t*state_t* c_t2 + noise.
	 */
	protected SparseVector[] coefficients = null;

	/**
	 * coefficient ~ N(B_t, sigma^2 * M_t): M_t is the covariates. the dimension
	 * of covariate is : (2* dimension) X (2*dimension)
	 */
	protected SparseVector[] varOfCoefficient = null;

	/**
	 * coefficient ~ N(B_t, sigma^2 * M_t): B_t is the means. the dimension of
	 * means is dimension*2.
	 */
	protected SparseVector[] meanOfCoefficient = null;

	/**
	 * The states changed over time. the dimension is : dimension.
	 */
	protected SparseVector[] states = null;

	/**
	 * the means for state at time t. i.e., m_t
	 */
	protected SparseVector[] meanOfState = null;

	/**
	 * the variance for state at time t. i.e., C_t (dim X dim)
	 */
	protected SparseVector[] varOfState = null;

	/** the weight for each particle */
	protected double[] particleWeights = null;

	/**
	 * cache some intermediate result for efficiency.
	 */

	SparseVector[] sigma2XB_2 = null;

	SparseVector[] sigma2XB_2XM = null;

	SparseVector[] miuB1 = null;

	SparseVector[] miuB2 = null;

	SparseVector[] vB1 = null;

	SparseVector[] vB2 = null;

	public boolean withoutSigma = true;

	/**
	 * Constructor
	 * 
	 * @param numOfParticles
	 */
	public TimeVaryingRegression(int numOfParticles, int dimension) {
		this.numOfParticles = numOfParticles;
		this.dimension = dimension;
	}

	protected abstract void computeWeights(Instance instance);

	protected abstract void computeWeights(long timestamp, Instance instance);

	protected abstract void cpFromSrcToDst(int src, int dst);

	protected abstract void updateStateStatistics(Instance instance);

	protected abstract void updateStateStatistics(long timestamp, Instance instance);

	protected abstract void sampleStates(long timestamp, Instance instance);

	protected abstract void sampleStates(Instance instance);

	protected abstract void updateParameterStatistics(Instance instance);

	protected abstract void updateParameterStatistics(long timestamp, Instance instance);

	protected abstract void sampleParameters(long timestamp, Instance instance);

	protected abstract void sampleParameters(Instance instance);

	public abstract double predict(long timestamp, SparseVector instance, double alpha);

	public abstract double predict(SparseVector instance, double alpha);

	protected void resample() {
		int[] frequency = Utility.sample(numOfParticles, particleWeights, numOfParticles);
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
				cpFromSrcToDst(src, dst);
				// the dst postion is filled
				dst++;
				// the number of elements at the src position filled.
				i++;
			}
			// the element at src is done.
			src++;
		}
	}

	@Override
	public void train(Collection<Instance> instances) {
		trainOnline(instances);
	}

	@Override
	public void train(TLongList timestamps, Collection<Instance> instances) {
		trainOnline(timestamps, instances);
	}

	public double predict(SparseVector instance, SparseVector weights) {
		return instance.innerProduct(weights);
	}

	public SparseVector getPosteriorMeans() {
		SparseVector ret = new SparseVector(dimension);
		for (int ip = 0; ip < numOfParticles; ip++) {

			ret.add(getPosteriorMean(ip, 1));
		}
		return ret.div(numOfParticles);
	}

	public SparseVector getPosteriorMeansFromParticles() {
		SparseVector ret = new SparseVector(dimension);
		for (int ip = 0; ip < numOfParticles; ip++) {
			ret.add(combineMeanStateAndCoefficient(ip));
		}
		return ret.div(numOfParticles);
	}

	@Override
	public SparseVector getPosteriorPrecisions() {

		SparseVector ret = new SparseVector(dimension);
		for (int ip = 0; ip < numOfParticles; ip++) {
			ret.add(getPosteriorVar(ip, 1));
		}
		return ret.div(numOfParticles * numOfParticles).inverse();
	}

	public SparseVector getPosteriorMean1(int indexP, long interval) {
		SparseVector C = plusOnes(varOfState[indexP], interval);
		SparseVector B2 = varOfCoefficient[indexP].subvector(dimension, 2 * dimension);
		SparseVector uB1 = meanOfCoefficient[indexP].subvector(0, dimension);
		SparseVector uB2 = meanOfCoefficient[indexP].subvector(dimension, dimension * 2);
		SparseVector term1 = B2.copyNew().mul(squareOfSigma[indexP]).add(C).inverse();
		SparseVector term2 = C.elementwiseMultiply(uB2)
				.add(B2.mul(squareOfSigma[indexP]).elementwiseMultiply(meanOfState[indexP]));
		term1 = uB1.add(term2.elementwiseMultiply(term1));
		return term1;

	}

	public SparseVector getPosteriorMean(int indexP, long interval) {
		interval = 0;
		SparseVector C = plusOnes(varOfState[indexP], interval);
		SparseVector term1 = C.copyNew().add(this.sigma2XB_2[indexP]).inverse();
		SparseVector term2 = C.elementwiseMultiply(this.miuB2[indexP]).add(this.sigma2XB_2XM[indexP]);
		term1 = term2.elementwiseMultiply(term1).add(this.miuB1[indexP]);
		return term1;

	}

	public SparseVector getPosteriorVar1(int indexP, long interval) {
		SparseVector C = plusOnes(varOfState[indexP], interval);
		SparseVector B1 = varOfCoefficient[indexP].subvector(0, dimension);
		SparseVector B2 = varOfCoefficient[indexP].subvector(dimension, 2 * dimension);
		SparseVector term1 = B2.copyNew().mul(squareOfSigma[indexP]).add(C).inverse();
		SparseVector term2 = B2.elementwiseMultiply(C).elementwiseMultiply(term1);
		term1 = B1.add(term2).mul(squareOfSigma[indexP]);
		return term1;
	}

	public SparseVector getPosteriorVar(int indexP, long interval) {
		interval = 0;
		SparseVector C = plusOnes(varOfState[indexP], interval);
		SparseVector term1 = C.copyNew().add(this.sigma2XB_2[indexP]).inverse();
		SparseVector term2 = C.elementwiseMultiply(this.vB2[indexP]).elementwiseMultiply(term1);
		term1 = term2.add(this.vB1[indexP]).mul(squareOfSigma[indexP]);
		return term1;
	}

	public double avgSigma2() {
		double sum = 0.0;
		for (int i = 0; i < numOfParticles; i++) {
			sum += squareOfSigma[i];
		}
		return sum / numOfParticles;
	}

	@Override
	public void outputDebugInfo(Writer writer) throws IOException {
		

	}

	public void sampleSquareOfSigma() {
		for (int i = 0; i < numOfParticles; i++) {
			if (!withoutSigma){
				this.squareOfSigma[i] = Utility.sampleInverseGammaDistribution(this.aPriorOfSquareSigma[i],
						this.bPriorOfSquareSigma[i]);
			}
			else{
				this.squareOfSigma[i] = 1.0;
			}
		}
	}

	protected SparseVector extendFeatureWithMeanState(SparseVector x, int indexP) {
		int dim = x.getDimension();

		SparseVector temp = x.elementwiseMultiply(meanOfState[indexP]);
		assert (x.getDefaultValue() == temp.getDefaultValue());
		SparseVector ret = new SparseVector(dim * 2, x.getDefaultValue());
		TIntIterator dimIter = temp.nonDefaultDimensions().iterator();
		int dimIndex;
		double val;
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();
			val = temp.get(dimIndex);
			ret.set(dimIndex + dim, val);
		}
		dimIter = x.nonDefaultDimensions().iterator();
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();
			val = x.get(dimIndex);
			ret.set(dimIndex, val);
		}
		return ret;
	}

	protected SparseVector extendFeatureWithState(SparseVector x, int indexP) {
		int dim = x.getDimension();
		SparseVector temp = x.elementwiseMultiply(states[indexP]);
		assert (x.getDefaultValue() == temp.getDefaultValue());
		SparseVector ret = new SparseVector(dim * 2, x.getDefaultValue());
		TIntIterator dimIter = temp.nonDefaultDimensions().iterator();
		int dimIndex;
		double val;
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();
			val = temp.get(dimIndex);
			ret.set(dimIndex + dim, val);
		}

		dimIter = x.nonDefaultDimensions().iterator();
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();
			val = x.get(dimIndex);
			ret.set(dimIndex, val);
		}
		return ret;

	}

	protected SparseVector combineMeanStateAndCoefficient(int indexP) {
		SparseVector front = coefficients[indexP].subvector(0, dimension);
		SparseVector rear = coefficients[indexP].subvector(dimension, 2 * dimension);
		front.add(rear.elementwiseMultiply(meanOfState[indexP]));
		return front;
	}

	protected SparseVector combineStateAndCoefficient(int indexP) {
		SparseVector front = coefficients[indexP].subvector(0, dimension);
		SparseVector rear = coefficients[indexP].subvector(dimension, 2 * dimension);
		front.add(rear.elementwiseMultiply(states[indexP]));
		return front;
	}

	protected SparseVector plusOnes(SparseVector x, long interval) {
		SparseVector ret = x.copyNew();
		return ret.add(interval);
	}

	protected SparseVector multFeatureWithStateCoef(SparseVector x, int ip) {
		SparseVector rear = coefficients[ip].subvector(dimension, 2 * dimension);
		return rear.elementwiseMultiply(x);
	}

	/**
	 * xAy
	 * 
	 * @param x
	 * @param A
	 * @param y
	 * @return
	 */
	protected double mult(SparseVector x, SparseVector A, SparseVector y) {
		return x.elementwiseMultiply(A).innerProduct(y);
	}

	protected SparseVector multCovAndCoef(long interval, int ip) {
		SparseVector I_plus_C = plusOnes(this.varOfState[ip], interval);
		SparseVector rear = coefficients[ip].subvector(dimension, 2 * dimension);
		return rear.square().elementwiseMultiply(I_plus_C);
	}

}
