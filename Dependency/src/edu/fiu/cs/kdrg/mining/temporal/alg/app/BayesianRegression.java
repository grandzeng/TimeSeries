package edu.fiu.cs.kdrg.mining.temporal.alg.app;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;

public interface BayesianRegression {

  double predict(SparseVector instance, SparseVector weights);

  SparseVector getPosteriorMeans();

  SparseVector getPosteriorPrecisions();
  
  SparseVector getPosteriorMeans(long timestamp);

  SparseVector getPosteriorPrecisions(long timestamp);

}
