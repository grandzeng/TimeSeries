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
public class BayesianLinearRegression implements OnlineBayesianRegression {

	private int dimension; // dimension of the feature

	private SparseVector weights; // Dimension -1 is the bias

	private SparseVector covariance; // X^TX

	private SparseVector b; // \sum_{i=1}{n}rx_n
	
	private SparseVector prioriMeans;
	
	private SparseVector prioriPrecision;

	private long numTrainedInstances = 0;

	public BayesianLinearRegression(int dimension) {
		this.dimension = dimension;
		this.weights = new SparseVector(dimension);
		this.covariance = new SparseVector(dimension, 1.0);
		this.b = new SparseVector(dimension);
		this.numTrainedInstances = 0;
	}

	public BayesianLinearRegression(int dimension, SparseVector prioriMeans, SparseVector prioriPrecision) {
		this.dimension = dimension;
		this.prioriMeans = prioriMeans.copyNew();
		this.prioriPrecision = prioriPrecision.copyNew();
		this.weights = prioriMeans.copyNew();
		this.covariance = prioriPrecision.copyNew().inverse();
		this.b = new SparseVector(dimension);
		this.numTrainedInstances = 0;
	}

	public void train(Collection<Instance> instances) {
		trainOnline(instances);
	}

	@Override
	public void trainOnline(Collection<Instance> instances) {
		for (Instance inst : instances) {
			SparseVector features = inst.getFeatures();
			covariance.add(features.copyNew().elementwiseMultiply(features));
			b.add(features.copyNew().mul(inst.getResponse()));
		}
		this.numTrainedInstances += instances.size();
		weights = covariance.copyNew().inverse().elementwiseMultiply(b);
	}

	@Override
	public double predict(SparseVector instance) {
		// TODO Auto-generated method stub
		return predict(instance, this.weights);
	}

	public double predict(SparseVector instance, SparseVector weights) {
		return instance.innerProduct(weights);
	}

	@Override
	public String getModelName() {
		return "LinearRegression_"+prioriMeans.getDefaultValue()+"_"+prioriPrecision.getDefaultValue();
	}

	@Override
	public SparseVector getPosteriorMeans() {
		// TODO Auto-generated method stub
		return weights;
	}

	@Override
	public SparseVector getPosteriorPrecisions() {
		// TODO Auto-generated method stub
		return covariance.copyNew();
	}

	@Override
	public OnlineRegression copyNew() {
		// TODO Auto-generated method stub
		BayesianLinearRegression copy = new BayesianLinearRegression(dimension, this.prioriMeans, prioriPrecision);
		copy.weights = this.weights.copyNew();
		copy.covariance = this.covariance.copyNew();
		copy.b = this.b.copyNew();
		copy.numTrainedInstances = this.numTrainedInstances;
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
