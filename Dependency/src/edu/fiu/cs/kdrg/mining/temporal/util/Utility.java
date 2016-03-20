package edu.fiu.cs.kdrg.mining.temporal.util;

import java.util.Random;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.ejml.simple.SimpleMatrix;

import gnu.trove.iterator.TIntIterator;

public class Utility {

	public static double epsilon = 0.000000001;

	public static NormalDistribution standardNormDistribution = new NormalDistribution(0.0, 1.0);

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
			if (Double.isNaN(probabilities[i])) {
				probabilities[i] = 0.0;
			}
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

	public static double sampleInverseGammaDistribution(double a, double b) {
		GammaDistribution gamma = new GammaDistribution(a, 1.0 / b);
		return 1.0 / gamma.sample();
		// return nextGamma(a, 1.0 / b, 0);
	}

	public static SimpleMatrix sampleExpDistribution(int size, ExponentialDistribution exp1,
			ExponentialDistribution exp2) {
		double[] temp1 = exp1.sample(size);
		double[] temp2 = exp2.sample(size);
		SimpleMatrix ret = SimpleMatrix.identity(2 * size);
		for (int i = 0; i < size; i++) {
			ret.set(i, i, Math.sqrt(temp1[i]));
			ret.set(i + size, i + size, Math.sqrt(temp2[i]));
		}
		return ret;
	}

	public static SimpleMatrix sampleExpDistribution(int size, ExponentialDistribution exp) {
		double[] temp1 = exp.sample(size);
		SimpleMatrix ret = SimpleMatrix.identity(size);
		for (int i = 0; i < size; i++) {
			ret.set(i, i, Math.sqrt(temp1[i]));
		}
		return ret;
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

	public static SimpleMatrix sampleMultivariateNormal(SimpleMatrix mean, SimpleMatrix cov) {
		int dim = mean.getNumElements();
		SimpleMatrix retMatrix = new SimpleMatrix(dim, 1);
		double[] ret = standardNormDistribution.sample(dim);
		double tempMean, tempVar;
		for (int i = 0; i < dim; i++) {
			tempMean = mean.get(i);
			tempVar = cov.get(i) == 0.0 ? epsilon : cov.get(i);
			retMatrix.set(i, tempMean + Math.sqrt(tempVar) * ret[i]);

		}
		return retMatrix;
	}

	public static SimpleMatrix sampleMultivariateGaussion(SimpleMatrix mean, SimpleMatrix cov) {
		double[] meanArray = new double[mean.getNumElements()];
		double[][] covArray = new double[cov.numRows()][cov.numCols()];
		for (int i = 0; i < meanArray.length; i++) {
			meanArray[i] = mean.get(i);
		}

		for (int i = 0; i < cov.numRows(); i++) {
			for (int j = 0; j < cov.numCols(); j++) {
				covArray[i][j] = cov.get(i, j);
			}
		}
		MultivariateNormalDistribution norm = new MultivariateNormalDistribution(meanArray, covArray);
		return new SimpleMatrix(meanArray.length, 1, true, norm.sample());
	}

	public static TDistribution createTDistribution(double degree) {
		return new TDistribution(degree);
	}
}
