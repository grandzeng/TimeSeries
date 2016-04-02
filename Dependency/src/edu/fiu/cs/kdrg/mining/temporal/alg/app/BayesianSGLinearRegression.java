package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.io.Writer;
import java.util.Collection;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TLongList;

/**
 * Batched Bayesian linear regression.
 * 
 * @author Liang Tang
 * @update Nov 25, 2013 10:58:18 AM
 * 
 */
public class BayesianSGLinearRegression implements OnlineBayesianRegression {

	private int dimension; // dimension of the feature

	private SparseVector weights; 

	private SparseVector weightPrecisions;

	private double regularizationWeight = -1;

	private SparseVector prioriMeans;

	private SparseVector prioriPrecisions;

	private long numTrainedInstances = 0;

	private final static int MAX_LOOP = 100;

	private final static double MIN_LEARNING_RATE = Double.MIN_VALUE * 10;

	protected double noise_var = 1.0;

	public BayesianSGLinearRegression(int dimension) {
		this(dimension, new SparseVector(dimension), new SparseVector(dimension));
	}

	public BayesianSGLinearRegression(int dimension, SparseVector prioriMeans, SparseVector prioriPrecisions) {
		this.dimension = dimension;
		this.weights = new SparseVector(dimension);
		this.weightPrecisions = new SparseVector(dimension);
		this.prioriMeans = new SparseVector(prioriMeans);
		this.prioriPrecisions = new SparseVector(prioriPrecisions);
		this.numTrainedInstances = 0;
	}

	public void train(Collection<Instance> instances) {
		double oldObjective = negativeAvgSquareError(instances, this.weights);
		double newObjective = -1;
		SparseVector newWeights = null;
		boolean bConverged = false;

		// Maximize the posterior estimation
		for (int iterCount = 0; iterCount < MAX_LOOP; iterCount++) {
			double learningRate = 0.5;
			// Find the right learning rate
			while (true) {
				newWeights = ascent(instances, learningRate, this.weights);
				newObjective = negativeAvgSquareError(instances, newWeights);
				if (newObjective < oldObjective) {
					learningRate /= 2.0;
					if (learningRate <= MIN_LEARNING_RATE) {
						bConverged = true;
						break;
					}
				} else if (learningRate <= MIN_LEARNING_RATE) {
					bConverged = true;
					break;
				} else {
					break;
				}
			}
			if (bConverged) {
				break;
			}

			// Update the weights
			this.weights = newWeights;
			oldObjective = newObjective;
			// System.out.println(DateUtil.getDate() + " objective : " +
			// newObjective + ", at iteration : "+ iterCount);
		}

		// Compute the precision(inverse of variance) of weights
		this.weightPrecisions = computePrecisions(weightPrecisions, instances, this.weights);
		numTrainedInstances = instances.size();
	}



	@Override
	public void trainOnline(Collection<Instance> instances) {

		// This learning-rate is decayed by the number of trained instances
		double learningRate = 10.0 / Math.sqrt(numTrainedInstances + 1);

		weights = ascent(instances, learningRate, weights);

		// compute the precisions of posterior distributions
		if (numTrainedInstances == 0) {
			weightPrecisions = new SparseVector(prioriPrecisions);
		}
		this.weightPrecisions = computePrecisions(weightPrecisions, instances, this.weights);
		numTrainedInstances += instances.size();
	}


	public void setRegluarizationWeight(double weight) {
		this.regularizationWeight = weight;
	}

	@Override
	public SparseVector getPosteriorMeans() {
		// TODO Auto-generated method stub
		if (numTrainedInstances == 0) {
			return this.prioriMeans;
		} else {
			return this.weights;
		}
	}

	@Override
	public SparseVector getPosteriorPrecisions() {
		// TODO Auto-generated method stub
		if (numTrainedInstances == 0) {
			return this.prioriPrecisions;
		} else {
			return weightPrecisions;
		}
	}

	public void setPrioriMeans(SparseVector means) {
		this.prioriMeans = means;
	}

	public void setPrioriPrecisions(SparseVector precisions) {
		this.prioriPrecisions = precisions;
	}

	/**
	 * Compute the precision from MAP estimated weights
	 * 
	 * @param instances
	 * @param weights
	 * @param w0
	 * @return
	 */
	protected SparseVector computePrecisions(SparseVector oldPrecision, Collection<Instance> instances,
			SparseVector weights) {
		SparseVector precisions = new SparseVector(dimension);
		for (Instance inst : instances) {
			SparseVector features = inst.getFeatures();
			TIntIterator dimIter = features.nonDefaultDimensions().iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				precisions.add(dimIndex, 1);
			}
		}
		SparseVector tmp = oldPrecision.copyNew();
		tmp.mul(noise_var);
		precisions.add(tmp);
		precisions.mul(1.0 / noise_var);
		return precisions;
	}

	/**
	 * Ascent the weights and w0 by using instances at one time
	 * 
	 * @param instances
	 * @param learningRate
	 * @param weights
	 * @param w0
	 * @return
	 */
	protected SparseVector ascent(Collection<Instance> instances, double learningRate, SparseVector weights) {
		SparseVector delta = new SparseVector(dimension);
		// ascent for instances
		for (Instance instance : instances) {
			SparseVector features = instance.getFeatures();
			double actual = instance.getResponse();
			double predicted = predict(features, weights);
			double diff = actual - predicted;
			TIntIterator dimIter = features.nonDefaultDimensions().iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				double weight = weights.get(dimIndex);
				double prioriPrecision = prioriPrecisions.get(dimIndex);
				double prioriMean = prioriMeans.get(dimIndex);
				double regularization = -prioriPrecision * (weight - prioriMean);
				if (regularizationWeight >= 0) {
					regularization = -regularizationWeight * (weight - prioriMean);
				}
				delta.add(dimIndex, learningRate * (diff + regularization));
			}
		}
		int numInsts = instances.size();
		if (numInsts > 0) {
			delta.div(instances.size());
		}
		SparseVector newWeights = new SparseVector(weights);
		newWeights.add(delta);
		return newWeights;
	}

	protected SparseVector ascent(SparseVector delta, SparseVector weights) {
		SparseVector newWeights = new SparseVector(weights);
		newWeights.add(delta);
		return newWeights;
	}

	protected double negativeAvgSquareError(Collection<Instance> instances, SparseVector weights) {
		double totalError = 0;
		for (Instance instance : instances) {
			SparseVector features = instance.getFeatures();
			double actual = instance.getResponse();
			double predicted = predict(features, weights);
			double error = actual - predicted;
			totalError += error * error;
		}
		return -totalError / instances.size();
	}

	@Override
	public double predict(SparseVector instance) {
		// TODO Auto-generated method stub
		return predict(instance, this.weights);
	}

	@Override
	public double predict(SparseVector features, SparseVector weights) {
		return features.innerProduct(weights);
	}

	

	@Override
	public String getModelName() {
		return "SGLinearRegression_"+prioriMeans.getDefaultValue()+"_"+prioriPrecisions.getDefaultValue();
	}

	@Override
	public OnlineRegression copyNew() {
		// TODO Auto-generated method stub
		BayesianSGLinearRegression copy = new BayesianSGLinearRegression(dimension, prioriMeans, prioriPrecisions);
		copy.noise_var = noise_var;
		copy.numTrainedInstances = numTrainedInstances;
		copy.regularizationWeight = regularizationWeight;
		copy.weights = weights.copyNew();
		copy.weightPrecisions = weightPrecisions.copyNew();
		return copy;
	}

	@Override
	public void outputDebugInfo(Writer writer) {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public void train(TLongList timestamps, Collection<Instance> instances) {
		train(instances);

	}

	@Override
	public double predict(long timestamp, SparseVector instance) {
		// TODO Auto-generated method stub
		return predict(instance);
	}

	@Override
	public void trainOnline(TLongList timestamps, Collection<Instance> instances) {
		// TODO Auto-generated method stub
		trainOnline(instances);
	}
	
	@Override
	public SparseVector getPosteriorMeans(long timestamp) {
		// TODO Auto-generated method stub
		return getPosteriorMeans();
	}

	@Override
	public SparseVector getPosteriorPrecisions(long timestamp) {
		// TODO Auto-generated method stub
		return getPosteriorPrecisions();
	}

}
