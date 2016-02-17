package edu.fiu.cs.kdrg.mining.temporal.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;



/**
 * K-means clustering
 * 
 * @author Liang Tang
 * @date Nov 14, 2013 9:12:34 PM
 */
public class KMeans {

  private final static int MAX_ITER = 500;

  private int k = 0;

  private SparseVector[] insts = null;

  public KMeans(SparseVector[] insts, int k) {
    this.insts = insts;
    this.k = k;
  }

  public int[] build() {
    SparseVector[] centers = initializeCenters();
    int[] instLabels = new int[insts.length];
    Arrays.fill(instLabels, 0);
    SparseVector[] lastCenters = Arrays.copyOf(centers, centers.length);
    Arrays.fill(lastCenters, null);
    for (int iter = 0; iter < MAX_ITER; iter++) {
      if (isCenterChanged(lastCenters, centers) == false) {
        break;
      }
      System.arraycopy(centers, 0, lastCenters, 0, centers.length);
      System.out.println("iter = " + iter + ", overall cluster sim = "
          + this.computeOverallClusterSimilarity(lastCenters, instLabels));
      updateInstLabels(centers, instLabels);
      updateCenters(instLabels, centers);
    }

    return instLabels;
  }

  private SparseVector[] initializeCenters() {
    Random rand = new Random(System.currentTimeMillis());
    Set<SparseVector> seedSet = new HashSet<SparseVector>();
    while (seedSet.size() < k) {
      int index = rand.nextInt(insts.length);
      SparseVector seed = insts[index];
      if (seedSet.contains(seed)) {
        continue;
      }
      seedSet.add(seed);
    }
    SparseVector[] seeds = new SparseVector[k];
    int i = 0;
    for (SparseVector seed : seedSet) {
      seeds[i] = seed;
      i++;
    }
    return seeds;
  }

  private void updateInstLabels(SparseVector[] centers, int[] instLabels) {
    for (int i = 0; i < insts.length; i++) {
      SparseVector inst = insts[i];
      double maxSim = -Double.MAX_VALUE;
      int bestCluLabel = -1;
      for (int j = 0; j < centers.length; j++) {
        double sim = SparseVector.cosine(centers[j], inst);
        if (sim > maxSim) {
          maxSim = sim;
          bestCluLabel = j;
        }
      }
      instLabels[i] = bestCluLabel;
    }
  }

  private void updateCenters(int[] instLabels, SparseVector[] centers) {
    int[] cluSizes = new int[k];
    for (int cluIdx = 0; cluIdx < centers.length; cluIdx++) {
      List<SparseVector> instsOfClu = new ArrayList<SparseVector>(this.insts.length / k);
      for (int j = 0; j < instLabels.length; j++) {
        if (instLabels[j] == cluIdx) {
          instsOfClu.add(insts[j]);
        }
      }
      cluSizes[cluIdx] = instsOfClu.size();

      SparseVector[] instArr = new SparseVector[instsOfClu.size()];
      instsOfClu.toArray(instArr);
      if (instArr.length > 0) {
        centers[cluIdx] = getMedian(instArr);
      } else {
        // Empty cluster
        centers[cluIdx] = null;
      }
    }
    for (int cluIdx = 0; cluIdx < centers.length; cluIdx++) {
      if (centers[cluIdx] == null) {
        // Find the farthest point to be the new center
        double minSim = Double.MAX_VALUE;
        int farthestInstIdx = -1;
        for (int instIdx = 0; instIdx < this.insts.length; instIdx++) {
          double sim = 0;
          SparseVector inst = this.insts[instIdx];
          for (int j = 0; j < centers.length; j++) {
            if (centers[j] != null) {
              sim += SparseVector.cosine(inst, centers[j]);
            }
          }
          if (sim < minSim) {
            minSim = sim;
            farthestInstIdx = instIdx;
          }
        }
        // Set the one to be the new center
        centers[cluIdx] = this.insts[farthestInstIdx];
        instLabels[farthestInstIdx] = cluIdx;
      }
    }
  }

  private double computeOverallClusterSimilarity(SparseVector[] centers, int[] instLabels) {
    double objVal = 0;
    for (int i = 0; i < this.insts.length; i++) {
      SparseVector inst = this.insts[i];
      SparseVector center = centers[instLabels[i]];
      objVal += SparseVector.cosine(inst, center);
    }
    return objVal;
  }

  private boolean isCenterChanged(SparseVector[] centers1, SparseVector[] centers2) {
    if (centers1.length != centers2.length) {
      return true;
    }
    for (int i = 0; i < centers1.length; i++) {
      if (centers1[i] == null && centers2[i] != null) {
        return true;
      } else if (centers1[i] != null && centers2[i] == null) {
        return true;
      }
    }

    double totalSim = 0;
    for (int i = 0; i < centers1.length; i++) {
      SparseVector center1 = centers1[i];
      double maxSim = -Double.MAX_VALUE;
      for (int j = 0; j < centers2.length; j++) {
        SparseVector center2 = centers2[j];
        double sim = SparseVector.cosine(center1, center2);
        if (sim >= maxSim) {
          maxSim = sim;
        }
      }
      totalSim += maxSim;
    }

    if (totalSim >= k - 0.0001) {
      return false;
    } else {
      return true;
    }
  }

  private static SparseVector getMedian(SparseVector[] insts) {
    if (insts == null || insts.length == 0) {
      throw new IllegalArgumentException("insts is null or empty");
    }
    int n = insts.length;
    int dimension = insts[0].getDimension();
    SparseVector center = new SparseVector(dimension);
    for (SparseVector inst : insts) {
      center.add(inst);
    }
    center.div(n);
    return center;
  }

}
