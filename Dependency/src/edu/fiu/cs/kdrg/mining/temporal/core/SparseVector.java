package edu.fiu.cs.kdrg.mining.temporal.core;

import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.ArrayList;
import java.util.List;

import edu.fiu.cs.kdrg.mining.temporal.util.MathFunctions;

/**
 * 
 * @author Liang Tang
 * @date Nov 14, 2013 10:10:47 PM
 */
public class SparseVector {

	TIntDoubleMap nonDefValues = new TIntDoubleHashMap();

	double defaultVal = 0;

	int dimension = -1;

	public SparseVector(int dimension) {
		this.dimension = dimension;
	}

	public SparseVector(int dimension, double defaultValue) {
		this.dimension = dimension;
		this.defaultVal = defaultValue;
	}

	public SparseVector(SparseVector copy) {
		this(copy, copy.defaultVal);
	}

	public SparseVector(double[] array) {
		this(array.length);
		for (int i = 0; i < array.length; i++) {
			this.set(i, array[i]);
		}
	}

	public SparseVector(SparseVector copy, double forceDefaultValue) {
		if (MathFunctions.almostEqual(copy.defaultVal, forceDefaultValue)) {
			this.defaultVal = copy.defaultVal;
			this.dimension = copy.dimension;
			this.nonDefValues.putAll(copy.nonDefValues);
		} else {
			this.defaultVal = forceDefaultValue;
			this.dimension = copy.dimension;
			for (int d = 0; d < copy.dimension; d++) {
				this.set(d, copy.get(d));
			}
		}
	}

	public SparseVector copyNew() {
		return new SparseVector(this);
	}

	/**
	 * Create a sparse vector from a string description, like {"1:10",
	 * "2:11",...}
	 * 
	 * @param descriptions
	 * @param dimension
	 * @return
	 */
	public static SparseVector create(String[] descriptions, int dimension) {
		SparseVector v = new SparseVector(dimension);
		for (String entry : descriptions) {
			entry = entry.trim();
			String[] token = entry.split(":");
			v.set(Integer.parseInt(token[0]), Double.parseDouble(token[1]));
		}
		return v;
	}

	public double getDefaultValue() {
		return defaultVal;
	}

	public int getDimension() {
		return this.dimension;
	}

	public void setDimension(int dimension) {
		this.dimension = dimension;
	}

	public void set(int dimIndex, double val) {
		boolean bContains = nonDefValues.containsKey(dimIndex);
		boolean bEqualDef = MathFunctions.almostEqual(val, defaultVal);
		if (bEqualDef) {
			if (bContains == false) {
				return;
			} else {
				nonDefValues.remove(dimIndex);
			}
		} else {
			nonDefValues.put(dimIndex, val);
		}
	}

	public void set(int dimIndex) {
		set(dimIndex, 1);
	}

	public void clear(int fromDimIndex, int endDimIndex) {
		List<Integer> dimToRemove = new ArrayList<Integer>();
		TIntIterator dimIter = nonDefValues.keySet().iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			if (dimIndex >= fromDimIndex && dimIndex < endDimIndex) {
				dimToRemove.add(dimIndex);
			}
		}
		for (Integer dimIndex : dimToRemove) {
			this.set(dimIndex, 0);
		}
	}

	public double get(int dimIndex) {
		if (nonDefValues.containsKey(dimIndex)) {
			return nonDefValues.get(dimIndex);
		} else {
			return defaultVal;
		}
	}

	public TIntSet nonDefaultDimensions() {
		return nonDefValues.keySet();
	}

	public SparseVector subvector(int start, int end) {
		int dim = end - start;
		if (start < 0 && end > dimension && dim <= 0) {
			throw new IllegalArgumentException("index out of range!");
		}
		SparseVector ret = new SparseVector(dim, this.defaultVal);
		TIntSet nonDefDim = new TIntHashSet(this.nonDefValues.keySet());
		TIntIterator dimIter = nonDefDim.iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			if (dimIndex < start || dimIndex >= end)
				continue;
			ret.set(dimIndex - start, get(dimIndex));
		}
		return ret;
	}

	public SparseVector add(SparseVector v) {
		if (MathFunctions.almostEqual(v.defaultVal, 0)) {
			TIntIterator dimIter = v.nonDefValues.keySet().iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				this.nonDefValues.put(dimIndex, this.get(dimIndex) + v.nonDefValues.get(dimIndex));
			}
		} else {
			TIntSet nonDefDimUnion = new TIntHashSet(this.nonDefValues.keySet());
			nonDefDimUnion.addAll(v.nonDefValues.keySet());
			TIntIterator dimIter = nonDefDimUnion.iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				double val = this.get(dimIndex) + v.get(dimIndex);
				this.nonDefValues.put(dimIndex, val);
			}
			this.defaultVal += v.defaultVal;
		}
		return this;
	}

	public SparseVector add(int dimIndex, double scalar) {
		this.set(dimIndex, this.get(dimIndex) + scalar);
		return this;
	}

	public SparseVector neg() {
		this.defaultVal = -this.defaultVal;
		TIntIterator dimIter = nonDefValues.keySet().iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			this.set(dimIndex, -this.get(dimIndex));
		}
		return this;
	}

	public SparseVector subtract(SparseVector v) {
		SparseVector tmp = v.copyNew();
		tmp.neg();
		return this.add(tmp);
	}

	public SparseVector plus(double scalar) {
		if (scalar == 0.0)
			return this;
		this.defaultVal += scalar;
		TIntSet dimIndices = new TIntHashSet(nonDefValues.keySet());
		TIntIterator dimIter = dimIndices.iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			this.set(dimIndex, get(dimIndex) + scalar);
		}
		return this;
	}

	public SparseVector mul(double scalar) {
		if (MathFunctions.almostEqual(scalar, 0)) {
			this.defaultVal = 0;
			this.nonDefValues = new TIntDoubleHashMap();
		} else {
			this.defaultVal *= scalar;
			TIntSet dimIndices = new TIntHashSet(nonDefValues.keySet());
			TIntIterator dimIter = dimIndices.iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				this.set(dimIndex, this.get(dimIndex) * scalar);
			}
		}
		return this;
	}

	public SparseVector div(double scalar) {
		return mul(1.0 / scalar);
	}

	public SparseVector pow(double p) {
		TIntSet dimIndices = new TIntHashSet(nonDefValues.keySet());
		this.defaultVal = Math.pow(defaultVal, p);
		TIntIterator dimIter = dimIndices.iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			this.set(dimIndex, Math.pow(this.get(dimIndex), p));
		}
		return this;
	}

	public SparseVector square() {
		TIntSet dimIndices = new TIntHashSet(nonDefValues.keySet());
		this.defaultVal = this.defaultVal * this.defaultVal;
		TIntIterator dimIter = dimIndices.iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			double val = this.get(dimIndex);
			this.set(dimIndex, val * val);
		}
		return this;
	}

	public SparseVector inverse() {
		TIntSet dimIndices = new TIntHashSet(nonDefValues.keySet());
		if (MathFunctions.almostEqual(this.defaultVal, 0)) {
			this.defaultVal = 0;
		} else {
			this.defaultVal = 1.0 / this.defaultVal;
		}
		TIntIterator dimIter = dimIndices.iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			double val = this.nonDefValues.get(dimIndex);
			if (MathFunctions.almostEqual(val, 0) == false) {
				this.set(dimIndex, 1.0 / val);
			}
		}
		return this;
	}

	// public SparseVector inverse() {
	// pow(-1);
	// return this;
	// }

	private double innerProductForBothZeroDefVal(SparseVector v) {
		SparseVector sparseVec = this;
		SparseVector denseVec = v;
		if (this.nonDefaultDimensions().size() > v.nonDefaultDimensions().size()) {
			sparseVec = v;
			denseVec = this;
		}
		TIntIterator dimIter = sparseVec.nonDefValues.keySet().iterator();
		double s = 0;
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			if (denseVec.nonDefValues.containsKey(dimIndex)) {
				s += sparseVec.nonDefValues.get(dimIndex) * denseVec.nonDefValues.get(dimIndex);
			}
		}
		return s;
	}

	private double innerProductForSingleZeroDefVal(SparseVector nonZeroDefVector) {
		TIntIterator dimIter = this.nonDefValues.keySet().iterator();
		double s = 0;
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			s += this.nonDefValues.get(dimIndex) * nonZeroDefVector.get(dimIndex);
		}
		return s;
	}

	public double innerProduct(SparseVector v) {
		if (this.dimension != v.dimension) {
			throw new IllegalArgumentException("Two vectors' dimensions are not identical !");
		}

		// For zero default values' sparse vectors, use this fast computation
		// method
		boolean bThisZeroDef = MathFunctions.almostEqual(this.defaultVal, 0);
		boolean bAnotherZeroDef = MathFunctions.almostEqual(v.defaultVal, 0);
		if (bThisZeroDef && bAnotherZeroDef) {
			return innerProductForBothZeroDefVal(v);
		} else if (bThisZeroDef && !bAnotherZeroDef) {
			return innerProductForSingleZeroDefVal(v);
		} else if (!bThisZeroDef && bAnotherZeroDef) {
			return v.innerProductForSingleZeroDefVal(this);
		} else {
			TIntSet dims = new TIntHashSet(nonDefValues.keySet());
			dims.addAll(v.nonDefValues.keySet());
			double s = 0;
			TIntIterator dimIter = dims.iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				double v1 = this.get(dimIndex);
				double v2 = v.get(dimIndex);
				s += v1 * v2;
			}
			int remainingDims = this.dimension - dims.size();
			s += remainingDims * (this.defaultVal * v.defaultVal);
			return s;
		}
	}

	/**
	 * do the self inverse and then compute the inner product with v
	 * 
	 * @param v
	 * @return
	 */
	public double inverseInnerProduct(SparseVector v) {
		if (this.dimension != v.dimension) {
			throw new IllegalArgumentException("Two vectors' dimensions are not identical !");
		}

		// For zero default values' sparse vectors, use this fast computation
		// method
		if (MathFunctions.almostEqual(v.defaultVal, 0)) {
			TIntIterator dimIter = v.nonDefValues.keySet().iterator();
			double s = 0;
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				s += v.nonDefValues.get(dimIndex) / this.get(dimIndex);
			}
			return s;
		} else {
			return this.copyNew().inverse().innerProduct(v);
		}
	}

	public double norm2() {
		return Math.sqrt(innerProduct(this));
	}

	private SparseVector elementwiseMultiplyForZeroDefVal(SparseVector v) {
		SparseVector sparseVec = this;
		SparseVector denseVec = v;
		if (this.nonDefaultDimensions().size() > v.nonDefaultDimensions().size()) {
			sparseVec = v;
			denseVec = this;
		}
		TIntIterator dimIter = sparseVec.nonDefValues.keySet().iterator();
		SparseVector ret = new SparseVector(v.dimension);
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			if (denseVec.nonDefValues.containsKey(dimIndex)) {
				double s = sparseVec.nonDefValues.get(dimIndex) * denseVec.nonDefValues.get(dimIndex);
				ret.set(dimIndex, s);
			}
		}
		return ret;
	}

	public SparseVector elementwiseMultiply(SparseVector v) {
		if (this.dimension != v.dimension) {
			throw new IllegalArgumentException("Two vectors' dimensions are not identical !");
		}
		// For zero default values' sparse vectors, use this fast computation
		// method
		if (MathFunctions.almostEqual(this.defaultVal, 0) && MathFunctions.almostEqual(v.defaultVal, 0)) {
			return elementwiseMultiplyForZeroDefVal(v);
		}

		SparseVector ret = new SparseVector(dimension, this.defaultVal * v.defaultVal);
		TIntSet dims = new TIntHashSet(nonDefValues.keySet());
		dims.addAll(v.nonDefValues.keySet());
		TIntIterator dimIter = dims.iterator();
		while (dimIter.hasNext()) {
			int dimIndex = dimIter.next();
			ret.set(dimIndex, this.get(dimIndex) * v.get(dimIndex));
		}
		return ret;
	}

	public static double cosine(SparseVector v1, SparseVector v2) {
		double inner = v1.innerProduct(v2);
		return inner / v1.norm2() / v2.norm2();
	}

	public double[] toDoubleArray() {
		double[] ret = new double[dimension];
		for (int index = 0; index < dimension; index++) {
			ret[index] = defaultVal;
		}

		TIntIterator iter = nonDefValues.keySet().iterator();
		while (iter.hasNext()) {
			int index = iter.next();
			ret[index] = nonDefValues.get(index);
		}

		return ret;
	}

	@Override
	public String toString() {
		return "def:" + defaultVal + ", dim:" + dimension + ", " + nonDefValues.toString();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(defaultVal);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + dimension;
		result = prime * result + ((nonDefValues == null) ? 0 : nonDefValues.hashCode());
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
		SparseVector other = (SparseVector) obj;
		if (dimension != other.dimension)
			return false;
		if (MathFunctions.almostEqual(this.defaultVal, other.defaultVal) == false) {
			other = new SparseVector(other, this.defaultVal);
		}

		if (nonDefValues == null) {
			if (other.nonDefValues != null)
				return false;
		} else if (!nonDefValues.equals(other.nonDefValues)) {
			return false;
		}
		return true;
	}

}
