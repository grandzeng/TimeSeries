package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.io.IOException;
import java.io.Writer;
import java.util.Collection;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;
import gnu.trove.list.TLongList;

/**
 * @author Chunqiu Zeng
 *
 */
public interface TemporalRegression {
	
	public void train(TLongList timestamps, Collection<Instance> instances);

	public double predict(long timestamp, SparseVector instance);

}
