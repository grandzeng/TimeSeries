package edu.fiu.cs.kdrg.mining.temporal.util;

import java.util.Random;

import org.apache.commons.math3.distribution.GammaDistribution;

public class MathFunctions {

	public static final double EPS = 1.0E-13;

	public static double crossEntropy(double target, double actual) {
		double adjustedTarget = (target == 0 ? 0.000001 : target);
		adjustedTarget = (target == 1.0 ? 0.999999 : target);
		double adjustedActual = (actual == 0 ? 0.000001 : actual);
		adjustedActual = (actual == 1 ? 0.999999 : actual);
		return -adjustedTarget * Math.log(adjustedActual) - (1 - adjustedTarget) * Math.log(1 - adjustedActual);
	}

	public static double crossEntropyDerivative(double target, double actual) {
		double adjustedTarget = target;
		double adjustedActual = actual;
		if (adjustedActual == 1) {
			adjustedActual = 0.999;
		} else if (actual == 0) {
			adjustedActual = 0.001;
		}
		if (adjustedTarget == 1) {
			adjustedTarget = 0.999;
		} else if (adjustedTarget == 0) {
			adjustedTarget = 0.001;
		}
		return -adjustedTarget / adjustedActual + (1 - adjustedTarget) / (1 - adjustedActual);
	}

	public static double sigmoid(double v) {
		return 1.0 / (1 + Math.pow(Math.E, -v));
	}

	public static double sigmoidDerivative(double v) {
		double sigmoid = sigmoid(v);
		return sigmoid * (1 - sigmoid);
	}

	public static double norm2(double[] v) {
		double sum = 0;
		for (double elem : v) {
			sum += elem * elem;
		}
		return Math.sqrt(sum);
	}

	public static double cosineSim(double[] v1, double[] v2) {
		if (v1.length != v2.length) {
			throw new IllegalArgumentException("The dimensions are not identifical !");
		} else {
			double s = 0;
			for (int i = 0; i < v1.length; i++) {
				s += v1[i] * v2[i];
			}
			double n2 = norm2(v2);
			double n1 = norm2(v1);
			if (n1 == 0) {
				throw new IllegalArgumentException("The norm of v1 is 0");
			}
			if (n2 == 0) {
				throw new IllegalArgumentException("The norm of v2 is 0");
			}
			return s / n1 / n2;
		}
	}

	public static boolean almostEqual(double x, double y) {
		double s = x - y;
		return (-EPS < s) && (s < EPS);
	}
	
	public static void normalizeProbabilites(double[] probs) {
		// Normalize the probabilities (make sure the summation to be 1)
		double sum = 0;
		for (int i=0; i<probs.length; i++) {
			sum += probs[i];
		}
		if (almostEqual(sum, 0)) {
			throw new IllegalArgumentException("The sum of probabilities is 0");
		}
		
		for (int i=0; i<probs.length; i++) {
			probs[i] /= sum;
		}
	}

	public static int randomSample(double[] probs, Random rand) {
		// Convert to cumulative probs
		double[] cumProbs = new double[probs.length];
		for (int i = 0; i < probs.length; i++) {
			if (i == 0) {
				cumProbs[i] = probs[i];
			} else {
				cumProbs[i] = cumProbs[i - 1] + probs[i];
			}
		}
		double v = rand.nextDouble();
		for (int i = 0; i < cumProbs.length; i++) {
			if (i == 0) {
				if (v <= cumProbs[i]) {
					return i;
				}
			} else {
				if (v > cumProbs[i - 1] && v <= cumProbs[i]) {
					return i;
				}
			}
		}
		return cumProbs.length - 1;
	}

	/**
	 * Generate the beta distribution sample via gamma distribution.
	 * 
	 * @param alpha
	 * @param beta
	 * @return
	 */
	public static double randomBeta(double alpha, double beta) {
		GammaDistribution gamma1 = new GammaDistribution(alpha, 1);
		GammaDistribution gamma2 = new GammaDistribution(beta, 1);
		double x = gamma1.sample();
		double y = gamma2.sample();
		return x / (x + y);
	}

}
