package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.util.Collection;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import gnu.trove.list.TLongList;

/**
 * Online training interface
 * 
 * @author Liang Tang
 * @date Jan 8, 2014 2:33:25 PM
 */
public interface OnlineRegression extends Regression, TemporalRegression {

	public void trainOnline(Collection<Instance> instances);

	public void trainOnline(TLongList timestamps, Collection<Instance> instances);

	public OnlineRegression copyNew();

}
