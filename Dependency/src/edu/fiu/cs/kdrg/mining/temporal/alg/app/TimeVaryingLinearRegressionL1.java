/**
 * 
 */
package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.util.Collection;

import org.apache.commons.math3.distribution.ExponentialDistribution;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;
import edu.fiu.cs.kdrg.mining.temporal.tv.util.Utility;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TLongList;

/**
 * Time Varying Linear Regression with L1 regularity.
 * 
 * @author Chunqiu Zeng
 * @date Jan. 22nd, 2016
 *
 */
public class TimeVaryingLinearRegressionL1 extends TimeVaryingLinearRegression {

	private double lambda1 = 1.0;

	private double lambda2 = 1.0;

	private SparseVector[] exponentialRandom = null;

	private ExponentialDistribution exp1 = null;

	private ExponentialDistribution exp2 = null;

	public TimeVaryingLinearRegressionL1(int numOfParticles, int dimension, double lambda1, double lambda2) {
		super(numOfParticles, dimension);
		this.lambda1 = lambda1;
		this.lambda2 = lambda2;
		this.exp1 = new ExponentialDistribution(2.0 / (lambda1 * lambda1));
		this.exp2 = new ExponentialDistribution(2.0 / (lambda2 * lambda2));
		exponentialRandom = new SparseVector[numOfParticles];
		for (int ip = 0; ip < numOfParticles; ip++) {
			exponentialRandom[ip] = new SparseVector(dimension * 2);
		}
	}

	public TimeVaryingLinearRegressionL1(int numOfParticles, int dimension, double priorMean, double priorPrecision,
			double lambda1, double lambda2) {
		super(numOfParticles, dimension, priorMean, priorPrecision);
		this.lambda1 = lambda1;
		this.lambda2 = lambda2;
		this.exp1 = new ExponentialDistribution(2.0 / (lambda1 * lambda1));
		this.exp2 = new ExponentialDistribution(2.0 / (lambda2 * lambda2));
		exponentialRandom = new SparseVector[numOfParticles];
		for (int ip = 0; ip < numOfParticles; ip++) {
			exponentialRandom[ip] = new SparseVector(dimension * 2);
		}
	}

	@Override
	public OnlineRegression copyNew() {
		TimeVaryingLinearRegressionL1 newObj = new TimeVaryingLinearRegressionL1(numOfParticles, dimension, lambda1,
				lambda2);
		copyTo(newObj);
		for (int ip = 0; ip < numOfParticles; ip++) {
			newObj.exponentialRandom[ip] = this.exponentialRandom[ip].copyNew();
		}
		return newObj;
	}

	@Override
	public String getModelName() {
		return "tvlrL1_" + lambda1 + "_" + lambda2 + "_" + priorMeanOfCoef + "_" + priorPrecisionOfCoef;
	}

	protected void sampleAugumentedVariables(int ip, Instance instance) {
		SparseVector x = instance.getFeatures();
		int dim = x.nonDefaultDimensions().size();
		double[] e1 = exp1.sample(dim);
		double[] e2 = exp2.sample(dim);
		TIntIterator dimIter = x.nonDefaultDimensions().iterator();
		int i = 0;
		int dimIndex;
		exponentialRandom[ip] = new SparseVector(dimension * 2);
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();
			exponentialRandom[ip].set(dimIndex, e1[i]);
			exponentialRandom[ip].set(dimIndex + dimension, e2[i]);
			i++;
		}
	}

	@Override
	protected void cpFromSrcToDst(int src, int dst) {
		super.cpFromSrcToDst(src, dst);
		exponentialRandom[dst] = exponentialRandom[src].copyNew();
	}

	@Override
	protected void updateParameterStatistics(Instance instance) {
		double y = instance.getResponse();
		SparseVector x = instance.getFeatures();
		for (int ip = 0; ip < numOfParticles; ip++) {
			SparseVector z = extendFeatureWithState(x, ip);
			SparseVector old_varOfCoef = varOfCoefficient[ip].copyNew();
			SparseVector temp = z.copyNew().square().elementwiseMultiply(this.exponentialRandom[ip]);
			varOfCoefficient[ip] = (varOfCoefficient[ip].inverse().add(temp)).inverse();

			SparseVector old_meanOfCoef = meanOfCoefficient[ip].copyNew();
			temp = old_varOfCoef.copyNew().inverse().elementwiseMultiply(old_meanOfCoef);
			meanOfCoefficient[ip] = varOfCoefficient[ip]
					.elementwiseMultiply(z.elementwiseMultiply(exponentialRandom[ip]).mul(y).add(temp));

			aPriorOfSquareSigma[ip] += 0.5;

			double term1 = y - z.innerProduct(old_meanOfCoef);
			term1 *= term1;
			temp = old_varOfCoef.elementwiseMultiply(this.exponentialRandom[ip]);
			double term2 = 1.0 + mult(z, temp, z);

			double term = (term1 / term2);
			bPriorOfSquareSigma[ip] += 0.5 * term;
		}

	}

	@Override
	protected void sampleParameters(long time, Instance instance) {

		sampleParameters(instance);

	}

	@Override
	protected void sampleParameters(Instance instance) {

		for (int ip = 0; ip < numOfParticles; ip++) {
			sampleAugumentedVariables(ip, instance);
			// SparseVector var =
			// varOfCoefficient[ip].elementwiseMultiply(exponentialRandom[ip]).mul(squareOfSigma[ip]);
			coefficients[ip] = Utility.sampleMultivariateNormalCoef(meanOfCoefficient[ip],
					varOfCoefficient[ip].elementwiseMultiply(exponentialRandom[ip]), squareOfSigma[ip],
					instance.getFeatures());
		}

	}
}
