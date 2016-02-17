package edu.fiu.cs.kdrg.mining.temporal.util;

/**
 * Basic vector operations.
 * 
 * @author yjian004
 * 
 */
public class DoubleVectorUtil {

  /**
   * Append a 1 in front of the 1st element of given vector.
   * 
   * @param vec
   * @return
   */
  public static double[] appendBias(double[] vec) {
    double[] res = new double[vec.length + 1];
    res[0] = 1;
    for (int i = 0; i < vec.length; ++i) {
      res[i + 1] = vec[i];
    }
    return res;
  }

  public static double innerProduct(double[] vec1, double[] vec2) {
    if (vec1.length != vec2.length) {
      throw new IllegalArgumentException("Dimension does not match.");
    }
    double res = 0;
    for (int i = 0; i < vec1.length; ++i) {
      res += vec1[i] * vec2[i];
    }
    return res;
  }

  public static double normN(double[] vec, int N) {
    double res = 0;
    for (int i = 0; i < vec.length; ++i) {
      res += Math.pow(Math.abs(vec[i]), N);
    }
    return Math.pow(res, 1.0 / N);
  }

  public static double norm2(double[] vec) {
    return normN(vec, 2);
  }

}
