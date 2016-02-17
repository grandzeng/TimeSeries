package edu.fiu.cs.kdrg.mining.temporal.core;

import java.util.List;

/**
 * Log entry for recommendation history
 * 
 * @author Liang Tang
 * @date Jan 16, 2014 4:13:08 AM
 */
public class RecommendLog {
	long timestamp;

	SparseVector contextFeatures = null;

	String servedItemId = null;

	double reward;

	List<String> itemIdPool = null;

	public RecommendLog() {

	}

	public RecommendLog(SparseVector contextFeatures, String servedItemId, double reward) {
		this.contextFeatures = contextFeatures;
		this.servedItemId = servedItemId;
		this.reward = reward;
	}

	public RecommendLog(SparseVector contextFeatures, String servedItemId, double reward, List<String> itemIdPool) {
		this(contextFeatures, servedItemId, reward);
		this.itemIdPool = itemIdPool;
	}

	public RecommendLog(SparseVector contextFeatures, String servedItemId, double reward, List<String> itemIdPool,
			long timestamp) {
		this(contextFeatures, servedItemId, reward, itemIdPool);
		this.timestamp = timestamp;
	}

	public long getTimestamp() {
		return timestamp;
	}

	public void setTimestamp(long timestamp) {
		this.timestamp = timestamp;
	}

	public SparseVector getContextFeatures() {
		return contextFeatures;
	}

	public String getServedItemId() {
		return servedItemId;
	}

	public double getReward() {
		return reward;
	}

	public void setReward(double reward) {
		this.reward = reward;
	}

	public List<String> getItemIdPool() {
		return itemIdPool;
	}

}
