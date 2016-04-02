package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import java.io.IOException;
import java.io.Writer;
import java.util.Collection;

import edu.fiu.cs.kdrg.mining.temporal.core.Instance;
import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;

/**
 * 
 * @author Liang Tang
 * @date Nov 16, 2013 10:52:48 PM
 */
public interface Regression {

	void train(Collection<Instance> instances);

	double predict(SparseVector instance);

	String getModelName();

	void outputDebugInfo(Writer writer) throws IOException;

}
