package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;
import edu.fiu.cs.kdrg.mining.temporal.tv.util.Utility;
import gnu.trove.list.TLongList;

public class TimeVaryingLinearRegression extends TimeVaryingRegression {

	private Random rand = new Random(System.nanoTime());

	double priorMeanOfCoef;

	double priorPrecisionOfCoef;

	long startTime = 0;

	long previousTime = 0;

	long currentTime = 0;

	public TimeVaryingLinearRegression(int numOfParticles, int dim) {
		this(numOfParticles, dim, 0.0, 1.0);
	}

	public TimeVaryingLinearRegression(int numOfParticles, int dim, double priorMean, double priorPrecision) {
		super(numOfParticles, dim);
		priorMeanOfCoef = priorMean;
		priorPrecisionOfCoef = priorPrecision;

		aPriorOfSquareSigma = new double[numOfParticles];
		Arrays.fill(aPriorOfSquareSigma, 1.0);
		bPriorOfSquareSigma = new double[numOfParticles];
		Arrays.fill(bPriorOfSquareSigma, 1.0);
		squareOfSigma = new double[numOfParticles];
		sampleSquareOfSigma();

		meanOfCoefficient = new SparseVector[numOfParticles];
		varOfCoefficient = new SparseVector[numOfParticles];
		coefficients = new SparseVector[numOfParticles];
		meanOfState = new SparseVector[numOfParticles];
		varOfState = new SparseVector[numOfParticles];
		states = new SparseVector[numOfParticles];
		particleWeights = new double[numOfParticles];

		// cache
		sigma2XB_2 = new SparseVector[numOfParticles];

		sigma2XB_2XM = new SparseVector[numOfParticles];

		miuB1 = new SparseVector[numOfParticles];

		miuB2 = new SparseVector[numOfParticles];

		vB1 = new SparseVector[numOfParticles];

		vB2 = new SparseVector[numOfParticles];

		for (int index = 0; index < numOfParticles; index++) {
			this.meanOfCoefficient[index] = new SparseVector(2 * dim, priorMean);
			this.varOfCoefficient[index] = new SparseVector(2 * dim, 1.0 / priorPrecision);
			// this.coefficients[index] =
			// Utility.sampleMultivariateNormal(meanOfCoefficient[index],
			// varOfCoefficient[index], squareOfSigma[index]);
			this.coefficients[index] = meanOfCoefficient[index].copyNew();
			this.meanOfState[index] = new SparseVector(dim);
			this.varOfState[index] = new SparseVector(dim, 1);
			this.states[index] = new SparseVector(dim);

			this.particleWeights[index] = 1.0 / numOfParticles;

			sigma2XB_2[index] = new SparseVector(dim, squareOfSigma[index] / priorPrecision);
			sigma2XB_2XM[index] = new SparseVector(dim);
			miuB1[index] = new SparseVector(dim);
			miuB2[index] = new SparseVector(dim);
			vB1[index] = new SparseVector(dim, 1.0 / priorPrecision);
			vB2[index] = new SparseVector(dim, 1.0 / priorPrecision);
		}

	}

	@Override
	public String getModelName() {
		return "tvlr_" + numOfParticles + "_" + priorMeanOfCoef + "_" + priorPrecisionOfCoef + this.withoutSigma;
	}

	@Override
	public void train(TLongList timestamps, Collection<Instance> instances) {
		trainOnline(timestamps, instances);

	}

	@Override
	protected void computeWeights(Instance instance) {
		this.computeWeightsWithInterval(instance, 1);
	}

	@Override
	protected void computeWeights(long timestamp, Instance instance) {
		computeWeightsWithInterval(instance, currentTime - previousTime);

	}

	protected void computeWeightsWithInterval(Instance instance, long interval) {
		double y = instance.getResponse();
		SparseVector x = instance.getFeatures();
		double mean;
		double var;
		SparseVector tempVec;
		double total = 0.0;
		for (int ip = 0; ip < numOfParticles; ip++) {
			tempVec = combineMeanStateAndCoefficient(ip);
			mean = tempVec.innerProduct(x);
			SparseVector I_plus_C = plusOnes(varOfState[ip], interval);
			SparseVector help = multFeatureWithStateCoef(x, ip);
			var = mult(help, I_plus_C, help) + squareOfSigma[ip];

			double std = Math.sqrt(var);
			particleWeights[ip] = Utility.standardNormDensity((y - mean) / std) / std;
			total += particleWeights[ip];
		}

		// normalize
		double defaultWeight = 1.0 / numOfParticles;
		for (int ip = 0; ip < numOfParticles; ip++) {
			if (total == 0.0) {
				particleWeights[ip] = defaultWeight;
			} else {
				particleWeights[ip] = particleWeights[ip] / total;
			}
		}
	}

	@Override
	protected void updateStateStatistics(Instance instance) {
		long interval = 1;
		kalmanFilterUpdate(instance, interval);
	}

	protected void kalmanFilterUpdate(Instance instance, long interval) {
		double y = instance.getResponse();
		SparseVector x = instance.getFeatures();
		for (int ip = 0; ip < numOfParticles; ip++) {
			// update Q_t
			// SparseVector I_plus_C = plusOnes(varOfState[ip], interval + 1);
			SparseVector I_plus_C = plusOnes(varOfState[ip], interval);
			SparseVector help = multFeatureWithStateCoef(x, ip);
			double Q_t = mult(help, I_plus_C, help) + squareOfSigma[ip];

			// update kalman gain A_t
			SparseVector A_t = I_plus_C.elementwiseMultiply(help).div(Q_t);

			// update meanOfState
			SparseVector temp = A_t.copyNew();
			double error = y - combineMeanStateAndCoefficient(ip).innerProduct(x);
			temp = temp.mul(error);
			meanOfState[ip] = meanOfState[ip].add(temp);

			// update varOfState
			varOfState[ip] = I_plus_C.sub(A_t.square().mul(Q_t));
		}
	}

	@Override
	protected void updateStateStatistics(long timestamp, Instance instance) {
		long interval = currentTime - previousTime;
		kalmanFilterUpdate(instance, interval);
	}

	@Override
	protected void sampleStates(long timestamp, Instance instance) {
		sampleStates(instance);
	}

	@Override
	protected void sampleStates(Instance instance) {
		SparseVector x = instance.getFeatures();
		for (int ip = 0; ip < numOfParticles; ip++) {
			states[ip] = Utility.sampleMultivariateNormalState(meanOfState[ip], varOfState[ip], 1.0, x);
		}
	}

	@Override
	protected void updateParameterStatistics(Instance instance) {
		double y = instance.getResponse();
		SparseVector x = instance.getFeatures();
		for (int ip = 0; ip < numOfParticles; ip++) {
			SparseVector z = extendFeatureWithState(x, ip);
			SparseVector old_varOfCoef = varOfCoefficient[ip].copyNew();
			varOfCoefficient[ip] = (varOfCoefficient[ip].inverse().add(z.copyNew().square())).inverse();

			SparseVector old_meanOfCoef = meanOfCoefficient[ip].copyNew();
			meanOfCoefficient[ip] = varOfCoefficient[ip].elementwiseMultiply(
					z.copyNew().mul(y).add(old_varOfCoef.copyNew().inverse().elementwiseMultiply(old_meanOfCoef)));

			aPriorOfSquareSigma[ip] += 0.5;

			double term1 = y - z.innerProduct(old_meanOfCoef);
			term1 *= term1;

			double term2 = 1.0 + mult(z, old_varOfCoef, z);

			double term = (term1 / term2);
			bPriorOfSquareSigma[ip] += 0.5 * term;
		}

	}

	@Override
	protected void updateParameterStatistics(long timestamp, Instance instance) {
		updateParameterStatistics(instance);
	}

	@Override
	protected void sampleParameters(Instance instance) {
		if (!this.withoutSigma)
			sampleSquareOfSigma();
		for (int ip = 0; ip < numOfParticles; ip++) {
			// SparseVector var =
			// varOfCoefficient[ip].copyNew().mul(squareOfSigma[ip]);
			coefficients[ip] = Utility.sampleMultivariateNormalCoef(meanOfCoefficient[ip], varOfCoefficient[ip],
					squareOfSigma[ip], instance.getFeatures());
		}

	}

	@Override
	protected void sampleParameters(long timestamp, Instance instance) {
		sampleParameters(instance);
	}

	@Override
	public void trainOnline(TLongList timestamps, Collection<Instance> instances) {
		// long start = System.currentTimeMillis();
		long tempTime = 0;
		int index = 0;
		for (Instance instance : instances) {
			tempTime = timestamps.get(index);
			if (impression == 0) {
				startTime = tempTime;
				previousTime = tempTime - 1;
				currentTime = tempTime;
			} else {
				previousTime = currentTime;
				currentTime = tempTime;
			}

			trainOnline(tempTime, instance);
			index++;
		}
	}

	protected void trainOnline(long timestamp, Instance instance) {
		impression++;
		if (instance.getResponse() > 0.0)
			click++;
		this.computeWeights(timestamp, instance);
		this.resample();
		this.updateStateStatistics(timestamp, instance);
		this.sampleStates(timestamp, instance);
		this.updateParameterStatistics(timestamp, instance);
		this.sampleParameters(timestamp, instance);
		// cache for efficiency
		this.updateCache();
	}

	@Override
	public void trainOnline(Collection<Instance> instances) {
		for (Instance instance : instances) {
			trainOnline(instance);
		}

	}

	protected void trainOnline(Instance instance) {
		impression++;
		if (instance.getResponse() > 0.0)
			click++;
		this.computeWeights(instance);
		this.resample();
		this.updateStateStatistics(instance);
		this.sampleStates(instance);
		this.updateParameterStatistics(instance);
		this.sampleParameters(instance);
		this.updateCache();
	}

	protected void updateCache() {
		for (int i = 0; i < numOfParticles; i++) {
			this.miuB1[i] = this.meanOfCoefficient[i].subvector(0, dimension);
			this.miuB2[i] = this.meanOfCoefficient[i].subvector(dimension, 2 * dimension);
			this.vB1[i] = this.varOfCoefficient[i].subvector(0, dimension);
			this.vB2[i] = this.varOfCoefficient[i].subvector(dimension, 2 * dimension);
			this.sigma2XB_2[i] = this.vB2[i].copyNew().mul(squareOfSigma[i]);
			this.sigma2XB_2XM[i] = this.sigma2XB_2[i].elementwiseMultiply(this.meanOfState[i]);

		}
	}

	@Override
	protected void cpFromSrcToDst(int src, int dst) {
		this.coefficients[dst] = this.coefficients[src].copyNew();
		this.states[dst] = this.states[src].copyNew();
		this.squareOfSigma[dst] = this.squareOfSigma[src];
		this.aPriorOfSquareSigma[dst] = this.aPriorOfSquareSigma[src];
		this.bPriorOfSquareSigma[dst] = this.bPriorOfSquareSigma[src];
		this.meanOfCoefficient[dst] = this.meanOfCoefficient[src].copyNew();
		this.varOfCoefficient[dst] = this.varOfCoefficient[src].copyNew();
		this.varOfState[dst] = this.varOfState[src].copyNew();
		this.meanOfState[dst] = this.meanOfState[src].copyNew();
		this.particleWeights[dst] = this.particleWeights[src];
		this.sigma2XB_2[dst] = this.sigma2XB_2[src].copyNew();
		this.sigma2XB_2XM[dst] = this.sigma2XB_2XM[src].copyNew();
		this.miuB1[dst] = this.miuB1[src].copyNew();
		this.miuB2[dst] = this.miuB1[src].copyNew();
		this.vB1[dst] = this.vB1[src].copyNew();
	}

	@Override
	public OnlineRegression copyNew() {
		TimeVaryingLinearRegression newObj = new TimeVaryingLinearRegression(numOfParticles, dimension);
		copyTo(newObj);
		return newObj;
	}

	protected void copyTo(TimeVaryingLinearRegression newObj) {
		newObj.aPriorOfSquareSigma = this.aPriorOfSquareSigma.clone();
		newObj.bPriorOfSquareSigma = this.bPriorOfSquareSigma.clone();
		newObj.click = this.click;
		newObj.coefficients = new SparseVector[numOfParticles];
		newObj.currentTime = this.currentTime;
		// newObj.dimension = this.dimension;
		newObj.impression = this.impression;
		newObj.meanOfCoefficient = new SparseVector[numOfParticles];
		newObj.meanOfState = new SparseVector[numOfParticles];
		// newObj.numOfParticles = this.numOfParticles;
		newObj.particleWeights = this.particleWeights.clone();
		newObj.previousTime = this.previousTime;
		newObj.priorMeanOfCoef = this.priorMeanOfCoef;
		newObj.priorPrecisionOfCoef = this.priorPrecisionOfCoef;
		newObj.squareOfSigma = this.squareOfSigma.clone();
		newObj.startTime = this.startTime;
		newObj.states = new SparseVector[numOfParticles];
		newObj.varOfCoefficient = new SparseVector[numOfParticles];
		newObj.varOfState = new SparseVector[numOfParticles];

		// cache
		newObj.sigma2XB_2 = new SparseVector[numOfParticles];

		newObj.sigma2XB_2XM = new SparseVector[numOfParticles];

		newObj.miuB1 = new SparseVector[numOfParticles];

		newObj.miuB2 = new SparseVector[numOfParticles];
		for (int i = 0; i < numOfParticles; i++) {
			newObj.coefficients[i] = this.coefficients[i].copyNew();
			newObj.meanOfCoefficient[i] = this.meanOfCoefficient[i].copyNew();
			newObj.meanOfState[i] = this.meanOfState[i].copyNew();
			newObj.states[i] = this.states[i].copyNew();
			newObj.varOfState[i] = this.varOfState[i].copyNew();
			newObj.varOfCoefficient[i] = this.varOfCoefficient[i].copyNew();

			newObj.sigma2XB_2[i] = this.sigma2XB_2[i].copyNew();
			newObj.sigma2XB_2XM[i] = this.sigma2XB_2XM[i].copyNew();
			newObj.miuB1[i] = this.miuB1[i].copyNew();
			newObj.miuB2[i] = this.miuB2[i].copyNew();
			newObj.vB1[i] = this.vB1[i].copyNew();
			newObj.vB2[i] = this.vB2[i].copyNew();
		}
	}

	public double predict(long timestamp, SparseVector instance, double alpha) {
		long interval = timestamp - currentTime;
		if (impression == 0)
			interval = 1;
		double mean = this.getPosteriorMeansFromParticles().innerProduct(instance);
		double sumVar = 0.0;
		for (int ip = 0; ip < numOfParticles; ip++) {
			SparseVector help = multCovAndCoef(interval, ip);
			sumVar += mult(instance, help, instance) + squareOfSigma[ip];
		}
		double std = Math.sqrt(sumVar) / numOfParticles;
		return mean + alpha * std;
	}

	public double predict(SparseVector instance, double alpha) {
		long interval = 1;
		double mean = getPosteriorMeansFromParticles().innerProduct(instance);
		double sumVar = 0.0;
		for (int ip = 0; ip < numOfParticles; ip++) {
			SparseVector help = multCovAndCoef(interval, ip);
			sumVar += mult(instance, help, instance) + squareOfSigma[ip];
		}
		double std = Math.sqrt(sumVar) / numOfParticles;
		return mean + alpha * std;
	}

	public double predict(long timestamp, SparseVector instance) {
		long interval = timestamp - currentTime;
		if (impression == 0)
			interval = 1;
		double[] ret = new double[numOfParticles];
		double[] sampleY = Utility.sampleStandardNormal(numOfParticles);
		double sum = 0.0;
		for (int ip = 0; ip < numOfParticles; ip++) {
			SparseVector tempVec = combineMeanStateAndCoefficient(ip);
			double mean = tempVec.innerProduct(instance);
			SparseVector I_plus_C = plusOnes(this.varOfState[ip], interval);
			SparseVector help = multFeatureWithStateCoef(instance, ip);
			double var = mult(help, I_plus_C, help) + squareOfSigma[ip];
			ret[ip] = sampleY[ip] * Math.sqrt(var) + mean;
			sum += ret[ip];
		}
		return sum / numOfParticles;
	}

	public double predict(SparseVector instance) {
		long interval = 1;
		double[] ret = new double[numOfParticles];
		double[] sampleY = Utility.sampleStandardNormal(numOfParticles);
		double sum = 0.0;
		for (int p = 0; p < numOfParticles; p++) {
			SparseVector x = combineMeanStateAndCoefficient(p);
			double mean = x.innerProduct(instance);
			SparseVector I_plus_C = plusOnes(this.varOfState[p], interval);
			SparseVector help = multFeatureWithStateCoef(instance, p);
			double var = mult(help, I_plus_C, help) + squareOfSigma[p];
			ret[p] = sampleY[p] * Math.sqrt(var) + mean;
			sum += ret[p];
		}
		return sum / numOfParticles;
	}

	@Override
	public SparseVector getPosteriorMeans(long timestamp) {
		long interval = timestamp - currentTime;
		if (impression == 0)
			interval = 0;
		interval = 0;
		SparseVector ret = new SparseVector(dimension);
		for (int ip = 0; ip < numOfParticles; ip++) {
			ret.add(getPosteriorMean(ip, interval));
		}
		return ret.div(numOfParticles);
	}

	@Override
	public SparseVector getPosteriorPrecisions(long timestamp) {
		long interval = timestamp - currentTime;
		if (impression == 0)
			interval = 0;
		interval = 0;
		SparseVector ret = new SparseVector(dimension);
		for (int ip = 0; ip < numOfParticles; ip++) {
			ret.add(getPosteriorVar(ip, interval));
		}
		return ret.div(numOfParticles * numOfParticles).inverse();
	}

}
