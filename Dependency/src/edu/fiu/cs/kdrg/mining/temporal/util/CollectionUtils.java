package edu.fiu.cs.kdrg.mining.temporal.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * 
 * @author Liang Tang
 * @update Nov 8, 2013 5:07:11 PM
 * 
 */
public class CollectionUtils {

	public static int[] asIntArray(Collection<Integer> l) {
		if (l.size() == 0) {
			return new int[0];
		} else {
			int[] arr = new int[l.size()];
			int i = 0;
			for (Integer val : l) {
				arr[i] = val;
				i++;
			}
			return arr;
		}
	}

	public static double[] asDoubleArray(Collection<Double> l) {
		if (l.size() == 0) {
			return new double[0];
		} else {
			double[] arr = new double[l.size()];
			int i = 0;
			for (Double val : l) {
				arr[i] = val;
				i++;
			}
			return arr;
		}
	}

	public static boolean[] asBoolArray(Collection<Boolean> l) {
		if (l.size() == 0) {
			return new boolean[0];
		} else {
			boolean[] arr = new boolean[l.size()];
			int i = 0;
			for (Boolean val : l) {
				arr[i] = val;
				i++;
			}
			return arr;
		}
	}

	public static List<Double> asList(double[] arr) {
		List<Double> l = new ArrayList<Double>(arr.length);
		for (double val : arr) {
			l.add(val);
		}
		return l;
	}
	
	public static int findLargestNumberIndex(long[] arr) {
		int largestIndex = -1;
		for (int i=0; i<arr.length; i++) {
			if (largestIndex == -1 || arr[largestIndex] < arr[i]) {
				largestIndex = i;				
			}
		}
		return largestIndex;
	}
}
