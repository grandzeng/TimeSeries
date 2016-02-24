package edu.fiu.cs.kdrg.mining.temporal.evaluation;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.util.MathFunctions;

/**
 * Compute AUC score: The area under curve.
 * 
 * @author Chunqiu Zeng
 * @date 22nd Feb, 2016
 *
 */
public class AUC {

	/**
	 * If the value is not zero, it gets a positive value.
	 * 
	 * @param value
	 * @return positive or negative
	 */
	public static boolean positive(double value) {
		return !MathFunctions.almostEqual(value, 0.0);
	}

	public static double compute(SimpleMatrix expected, SimpleMatrix actual) {
		int dim = expected.getNumElements();
		int numPositives = 0;
		int numNegatives = 0;
		for (int i = 0; i < dim; i++) {
			if (!positive(expected.get(i))) {
				numNegatives++;
			} else {
				numPositives++;
			}
		}
		// sort actual by the absolute value of its elements.
		List<Integer> indices = sortInDecendingOrder(actual);
		int TP = 0;
		int area = 0;
		for (Integer index : indices) {
			if (positive(expected.get(index))) {
				TP++;
			} else {
				area += TP;
			}
		}
		double prod = numPositives * numNegatives;
		if (prod == 0.0)
			return 0.0;
		return area / prod;
	}

	public static List<Integer> sortInDecendingOrder(SimpleMatrix vec) {
		List<Integer> indices = new ArrayList<Integer>();
		for (int i = 0; i < vec.getNumElements(); i++) {
			indices.add(i);
		}
		Collections.sort(indices, new IndexComparator(vec));
		Collections.reverse(indices);
		return indices;
	}

	static class IndexComparator implements Comparator<Integer> {

		private SimpleMatrix vec;

		public IndexComparator(SimpleMatrix v) {
			vec = v;
		}

		@Override
		public int compare(Integer o1, Integer o2) {
			double ret = Math.abs(vec.get(o1)) - Math.abs(vec.get(o2));
			if (ret == 0.0)
				return 0;
			else if (ret > 0.0)
				return 1;
			return -1;
		}
	}
}
