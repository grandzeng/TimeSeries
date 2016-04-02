package edu.fiu.cs.kdrg.mining.temporal.tv.util;

import java.util.Random;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseMatrix;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;

import gnu.trove.iterator.TIntIterator;

public class Utility {

	public static double epsilon = 0.000000001;

	public static SparseVector expOnes = null;

	private static final NormalDistribution standardNormDistribution = new NormalDistribution();
	private static final Random rand = new Random();

	/**
	 * This function is used to compute the sum of all the values in a given
	 * array.
	 * 
	 * @param array
	 * @return
	 */
	public static double sum(double[] array) {
		double ret = 0.0;
		for (double v : array) {
			ret += v;
		}
		return ret;
	}

	public static double squareSum(double[] array) {
		double ret = 0.0;
		for (double v : array) {
			ret += v * v;
		}
		return ret;
	}

	/**
	 * This function is used to compute the average value, given an array of
	 * values.
	 * 
	 * @param array
	 * @return
	 */
	public static double avg(double[] array) {
		double ret = sum(array);
		return ret / array.length;
	}

	/**
	 * This function is used to sample a set of particles from a given list of
	 * particles according to a given list of probabilities.
	 * 
	 * @param numOfParticles
	 * @param probabilities
	 * @param count
	 * @return
	 */
	public static int[] sample(int numOfParticles, double probabilities[], int count) {
		int frequency[] = new int[numOfParticles];
		int singletons[] = new int[numOfParticles];
		double total = 0.0;
		for (int i = 0; i < numOfParticles; i++) {
			frequency[i] = 0;
			singletons[i] = i;
			// if (Double.isNaN(probabilities[i])) {
			// probabilities[i] = 0.0;
			// }
			total += probabilities[i];
		}
		if (total == 0.0) {
			for (int i = 0; i < numOfParticles; i++)
				probabilities[i] = 1.0 / numOfParticles;
		}
		EnumeratedIntegerDistribution distribution = new EnumeratedIntegerDistribution(singletons, probabilities);
		for (int i = 0; i < count; i++) {
			int x = distribution.sample();
			frequency[x] += 1;
		}

		return frequency;
	}

	/**
	 * This function is used to compute the density of standard normal
	 * distribution, given a real number x.
	 * 
	 * @param x
	 * @return
	 */
	public static double standardNormDensity(double x) {
		// return standardNormDistribution.density(x);
		return Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
	}

	public static double[] sampleStandardNormal(int size) {
		return standardNormDistribution.sample(size);
	}

	public static SparseVector sampleMultivariateNormal(SparseVector mean, SparseVector cov) {
		int dim = mean.getDimension();
		double[] ret = standardNormDistribution.sample(dim);
		double tempMean, tempVar;
		for (int i = 0; i < dim; i++) {
			tempMean = mean.get(i);
			tempVar = cov.get(i) == 0.0 ? epsilon : cov.get(i);
			ret[i] = tempMean + Math.sqrt(tempVar) * ret[i];
		}
		return new SparseVector(ret);
	}

	public static SparseVector sampleMultivariateNormalState(SparseVector mean, SparseVector cov, double SquareOfSigma,
			SparseVector features) {
		int dim = features.nonDefaultDimensions().size();
		// int dimension = features.getDimension();
		double[] probs = standardNormDistribution.sample(dim);
		assert (MathFunctions.almostEqual(features.getDefaultValue(), 0));
		SparseVector ret = new SparseVector(mean.getDimension());
		TIntIterator dimIter = features.nonDefaultDimensions().iterator();
		int dimIndex;
		int i = 0;
		double help, tempVar;
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();

			help = cov.get(dimIndex) * SquareOfSigma;
			tempVar = (help == 0.0 ? epsilon : help);
			help = mean.get(dimIndex) + Math.sqrt(tempVar) * probs[i];
			ret.set(dimIndex, help);
			i++;
		}
		return ret;
	}

	public static SparseVector sampleMultivariateNormalCoef(SparseVector mean, SparseVector cov, double SquareOfSigma,
			SparseVector features) {
		int dim = features.nonDefaultDimensions().size();
		int dimension = features.getDimension();
		double[] probs = standardNormDistribution.sample(2 * dim);
		assert (MathFunctions.almostEqual(features.getDefaultValue(), 0));
		SparseVector ret = new SparseVector(mean.getDimension());
		TIntIterator dimIter = features.nonDefaultDimensions().iterator();
		int dimIndex;
		int i = 0;
		double help, tempVar;
		while (dimIter.hasNext()) {
			dimIndex = dimIter.next();

			help = cov.get(dimIndex) * SquareOfSigma;
			tempVar = (help == 0.0 ? epsilon : help);
			help = mean.get(dimIndex) + Math.sqrt(tempVar) * probs[i];
			ret.set(dimIndex, help);

			help = cov.get(dimIndex + dimension) * SquareOfSigma;
			tempVar = (help == 0.0 ? epsilon : help);
			help = mean.get(dimIndex + dimension) + Math.sqrt(tempVar) * probs[i + dim];
			ret.set(dimIndex + dimension, help);
			i++;
		}
		return ret;
	}

	public static SparseVector sampleExpDistribution(int size, ExponentialDistribution exp1,
			ExponentialDistribution exp2) {
		double[] temp1 = exp1.sample(size);
		double[] temp2 = exp2.sample(size);
		SparseVector ret = new SparseVector(2 * size);
		for (int i = 0; i < size; i++) {
			ret.set(i, temp1[i]);
			ret.set(i + size, temp2[i]);
		}
		return new SparseVector(ret);
	}

	public static double sampleInverseGammaDistribution(double a, double b) {
		GammaDistribution gamma = new GammaDistribution(a, 1.0 / b);
		return 1.0 / gamma.sample();
		// return nextGamma(a, 1.0 / b, 0);
	}

	public static double nextGamma(double alpha, double beta, double lambda) {
		double gamma = 0;
		if (alpha <= 0 || beta <= 0) {
			throw new IllegalArgumentException("alpha and beta must be strictly positive.");
		}
		if (alpha < 1) {
			double b, p;
			boolean flag = false;
			b = 1 + alpha * Math.exp(-1);
			while (!flag) {
				p = b * rand.nextDouble();
				if (p > 1) {
					gamma = -Math.log((b - p) / alpha);
					if (rand.nextDouble() <= Math.pow(gamma, alpha - 1))
						flag = true;
				} else {
					gamma = Math.pow(p, 1 / alpha);
					if (rand.nextDouble() <= Math.exp(-gamma))
						flag = true;
				}
			}
		} else if (alpha == 1) {
			gamma = -Math.log(rand.nextDouble());
		} else {
			double y = -Math.log(rand.nextDouble());
			while (rand.nextDouble() > Math.pow(y * Math.exp(1 - y), alpha - 1))
				y = -Math.log(rand.nextDouble());
			gamma = alpha * y;
		}
		return beta * gamma + lambda;
	}

}
