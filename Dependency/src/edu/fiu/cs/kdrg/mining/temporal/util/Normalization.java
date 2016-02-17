package edu.fiu.cs.kdrg.mining.temporal.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Normalization {

  public static List<double[]> zeroOneNormalization(List<double[]> instances) {
    List<double[]> normalized = new ArrayList<double[]>();
    int dimension = instances.get(0).length;
    double[] mins = new double[dimension];
    double[] maxs = new double[dimension];

    Arrays.fill(mins, Double.MAX_VALUE);
    Arrays.fill(maxs, Double.MIN_VALUE);

    for (double[] instance : instances) {
      for (int d = 0; d < instance.length; ++d) {
        if (instance[d] > maxs[d]) {
          maxs[d] = instance[d];
        }
        if (instance[d] < mins[d]) {
          mins[d] = instance[d];
        }
      }
    }

    for (double[] instance : instances) {
      double[] normalizedInstance = new double[dimension];
      for (int d = 0; d < dimension; ++d) {
        normalizedInstance[d] = (instance[d] - mins[d]) / (maxs[d] - mins[d]);
      }
      normalized.add(normalizedInstance);
    }

    return normalized;
  }

}
