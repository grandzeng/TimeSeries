package edu.fiu.cs.kdrg.mining.temporal.core;

/**
 * Generic instance object for regression and classification algorithms
 * 
 * @author Liang Tang
 * @date Dec 19, 2013 12:31:35 PM
 */
public class Instance {

	SparseVector features = null;

	double response;

	public Instance() {

	}

	public Instance(SparseVector features, double response) {
		this.features = features;
		this.response = response;
	}

	public void setFeatures(SparseVector features) {
		this.features = features;
	}

	public void setResponse(double response) {
		this.response = response;
	}

	public SparseVector getFeatures() {
		return features;
	}

	public double getResponse() {
		return response;
	}
	
	public String toString() {
		return features.toString()+" : "+response;
	}
}
