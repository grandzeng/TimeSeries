package edu.fiu.cs.kdrg.mining.temporal.util;

import java.util.Comparator;

/**
 * To implement the comparator for pairs (comparing by the right element)
 * 
 * @author Liang Tang
 * @date Dec 20, 2013 12:04:35 AM
 * @param <T>
 */
public class PairComparator<T> implements Comparator<Pair<T, Double>> {

  private int reverse; // 1 means sorting in natural order, -1 means sorting in
                       // reverse order

  public PairComparator() {
    this.reverse = 1;
  }

  public PairComparator(boolean reverse) {
    this.reverse = reverse ? -1 : 1;
  }

  @Override
  public int compare(Pair<T, Double> o1, Pair<T, Double> o2) {
    if (o1.getRight() > o2.getRight()) {
      return this.reverse * 1;
    } else if (o1.getRight() < o2.getRight()) {
      return this.reverse * -1;
    } else {
      return 0;
    }
  }

}
