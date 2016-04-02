package edu.fiu.cs.kdrg.mining.temporal.core;

import org.ejml.simple.SimpleMatrix;

public class TimeFrame {

	private int org = 0;

	private int lag;

	private int dimension;

	private SimpleMatrix vector;

	public TimeFrame() {
		vector = new SimpleMatrix();
	}

	public TimeFrame(int lag, int dimension) {
		this.lag = lag;
		this.dimension = dimension;
		this.vector = new SimpleMatrix(lag * dimension, 1);
	}

	public TimeFrame copy() {
		TimeFrame copy = new TimeFrame(lag, dimension);
		copy.org = org;
		copy.vector = vector.copy();
		return copy;
	}

	/**
	 * Get the corresponding value with given lag index and dimension index.
	 * 
	 * @param l
	 *            the lag index
	 * @param d
	 *            the dimension index
	 * @return the corresponding value
	 */
	public double getValue(int l, int d) {
		assert (0 <= l && l < lag);
		assert (0 <= d && d < dimension);
		int lagIndex = (l + org) % lag;
		return vector.get(lagIndex * dimension + d);
	}

	/**
	 * Set the corresponding value with given lag index and dimension index.
	 * 
	 * @param l
	 *            the lag index
	 * @param d
	 *            the dimension index
	 * @param v
	 *            the value to be set
	 * @return the corresponding value
	 */
	public void setValue(int l, int d, double v) {
		assert (0 <= l && l < lag);
		assert (0 <= d && d < dimension);
		int lagIndex = (l + org) % lag;
		vector.set(lagIndex * dimension + d, v);
	}

	/**
	 * Append the latest values to the time frame and set the lag index of the
	 * new values to 0. Increase the lag indices of all the old values. lag *
	 * 
	 * @param values
	 *            the new appended values
	 */
	public void shift(SimpleMatrix values) {
		assert (values.getNumElements() == dimension);
		int l = lag - 1;
		double v;
		for (int d = 0; d < dimension; d++) {
			v = values.get(d);
			setValue(l, d, v);
		}
		org = (org + lag - 1) % lag;
	}

	/**
	 * Append the latest values to the time frame and set the lag index of the
	 * new values to 0. Increase the lag indices of all the old values. lag *
	 * 
	 * @param values
	 *            the new appended values
	 */
	public void shift(double[] values) {
		assert (values.length == dimension);
		int l = lag - 1;
		double v;
		for (int d = 0; d < dimension; d++) {
			v = values[d];
			setValue(l, d, v);
		}
		org = (org + lag - 1) % lag;
	}

	/**
	 * Get the vector of values.
	 * 
	 * @return vector
	 */
	public SimpleMatrix getVector() {
		SimpleMatrix ret = new SimpleMatrix(lag * dimension, 1);
		double value;
		for (int l = 0; l < lag; l++) {
			for (int d = 0; d < dimension; d++) {
				value = getValue(l, d);
				ret.set(l * dimension + d, value);
			}
		}
		return ret;
	}
	
	public SparseVector getSparseVector(){
		SparseVector sv = new SparseVector(lag*dimension,0);
		double value;
		for(int l = 0;l<lag;l++){
			for(int d=0;d<dimension;d++){
				value= getValue(l,d);
				sv.set(l*dimension+d, value);
			}
		}
		return sv;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + dimension;
		result = prime * result + lag;
		double val;
		for (int l = 0; l < lag; l++) {
			for (int d = 0; d < dimension; d++) {
				val = this.getValue(l, d);
				result = prime * result + Double.valueOf(val).hashCode();
			}
		}
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		TimeFrame other = (TimeFrame) obj;
		if (dimension != other.dimension)
			return false;
		if (lag != other.lag)
			return false;
		for (int l = 0; l < lag; l++) {
			for (int d = 0; d < dimension; d++) {
				if (this.getValue(l, d) != other.getValue(l, d))
					return false;
			}
		}
		return true;
	}

	/**
	 * @return the lag
	 */
	public int getLag() {
		return lag;
	}

	/**
	 * @param lag
	 *            the lag to set
	 */
	public void setLag(int lag) {
		this.lag = lag;
	}

	/**
	 * @return the dimension
	 */
	public int getDimension() {
		return dimension;
	}

	/**
	 * @param dimension
	 *            the dimension to set
	 */
	public void setDimension(int dimension) {
		this.dimension = dimension;
	}

}
