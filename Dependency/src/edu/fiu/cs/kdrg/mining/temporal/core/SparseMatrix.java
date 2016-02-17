package edu.fiu.cs.kdrg.mining.temporal.core;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import edu.fiu.cs.kdrg.mining.temporal.util.MathFunctions;
import edu.fiu.cs.kdrg.mining.temporal.util.Pair;
import gnu.trove.iterator.TIntIterator;

/**
 * Sparse matrix
 * 
 * @author Liang Tang
 * @date Jan 19, 2014 8:11:17 PM
 */
public class SparseMatrix {

	int nRows = -1;

	int nCols = -1;

	double defaultVal = 0;

	Map<Pair<Integer, Integer>, Double> values = new HashMap<Pair<Integer, Integer>, Double>();

	public SparseMatrix(int nRows, int nCols) {
		this.nRows = nRows;
		this.nCols = nCols;
	}

	public SparseMatrix(int nRows, int nCols, double defaultVal) {
		this.nRows = nRows;
		this.nCols = nCols;
		this.defaultVal = defaultVal;
	}

	public SparseMatrix(SparseMatrix copy) {
		this.nRows = copy.nRows;
		this.nCols = copy.nCols;
		this.defaultVal = copy.defaultVal;
		this.values.putAll(copy.values);
	}

	public SparseMatrix(SparseVector vector, boolean isColumnVector) {
		if (isColumnVector) {
			TIntIterator dimIter = vector.nonDefaultDimensions().iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				this.set(dimIndex, 0, vector.get(dimIndex));
			}
			nRows = vector.getDimension();
			nCols = 1;
		} else {
			TIntIterator dimIter = vector.nonDefaultDimensions().iterator();
			while (dimIter.hasNext()) {
				int dimIndex = dimIter.next();
				this.set(0, dimIndex, vector.get(dimIndex));
			}
			nCols = vector.getDimension();
			nRows = 1;
		}
		this.defaultVal = vector.defaultVal;
	}

	public int getNumRows() {
		return this.nRows;
	}

	public int getNumCols() {
		return this.nCols;
	}

	public double getDefaultValue() {
		return this.defaultVal;
	}

	public Set<Pair<Integer, Integer>> getNonDefaultDimensions() {
		return values.keySet();
	}

	public double get(Pair<Integer, Integer> loc) {
		Double val = values.get(loc);
		if (val == null) {
			return defaultVal;
		} else {
			return val;
		}
	}

	public double get(int row, int col) {
		return get(new Pair<Integer, Integer>(row, col));
	}

	public void set(Pair<Integer, Integer> loc, double val) {
		if (values.containsKey(loc) && MathFunctions.almostEqual(val, defaultVal)) {
			values.remove(loc);
		} else {
			values.put(loc, val);
		}
	}

	public void set(int row, int col, double val) {
		set(new Pair<Integer, Integer>(row, col), val);
	}

	public SparseMatrix mul(double scalar) {
		if (MathFunctions.almostEqual(scalar, 0)) {
			this.defaultVal = 0;
			values.clear();
		} else {
			defaultVal *= scalar;
			for (Pair<Integer, Integer> key : values.keySet()) {
				values.put(key, values.get(key) * scalar);
			}
		}
		return this;
	}

	public SparseMatrix add(SparseMatrix other) {
		for (Pair<Integer, Integer> key : other.values.keySet()) {
			this.set(key, this.get(key) + other.get(key));
		}
		this.defaultVal += other.defaultVal;
		return this;
	}

	public SparseVector diagonal() {
		if (nCols != nRows) {
			throw new IllegalAccessError("This matrix is not diagonal matrix!");
		}
		SparseVector ret = new SparseVector(nCols);
		for (Pair<Integer, Integer> key : values.keySet()) {
			if (key.getLeft() == key.getRight()) {
				ret.set(key.getLeft(), values.get(key));
			}
		}
		ret.defaultVal = this.defaultVal;
		return ret;
	}

	/**
	 * The key is the row index, the value is the column index with the entry
	 * value
	 * 
	 * @return
	 */
	private Map<Integer, Map<Integer, Double>> getRowColMap() {
		Map<Integer, Map<Integer, Double>> rowColMap = new HashMap<Integer, Map<Integer, Double>>();
		for (Pair<Integer, Integer> key : this.getNonDefaultDimensions()) {
			int row = key.getLeft();
			int col = key.getRight();
			Map<Integer, Double> colMap = rowColMap.get(row);
			if (colMap == null) {
				colMap = new HashMap<Integer, Double>();
				rowColMap.put(row, colMap);
			}
			colMap.put(col, this.get(key));
		}
		return rowColMap;
	}

	/**
	 * The key is the column index, the value is the row index with the entry
	 * value
	 * 
	 * @return
	 */
	private Map<Integer, Map<Integer, Double>> getColRowMap() {
		Map<Integer, Map<Integer, Double>> colRowMap = new HashMap<Integer, Map<Integer, Double>>();
		for (Pair<Integer, Integer> key : this.getNonDefaultDimensions()) {
			int row = key.getLeft();
			int col = key.getRight();
			Map<Integer, Double> rowMap = colRowMap.get(col);
			if (rowMap == null) {
				rowMap = new HashMap<Integer, Double>();
				colRowMap.put(col, rowMap);
			}
			rowMap.put(row, this.get(key));
		}
		return colRowMap;
	}

	/**
	 * the product of two matrices, where the other matrix is on the right side
	 * 
	 * @param other
	 * @return
	 */
	public SparseMatrix rightProduct(SparseMatrix other) {
		if (this.nCols != other.nRows) {
			throw new IllegalArgumentException("This matrix's column is not equal to the other matrix's row");
		}
		if (MathFunctions.almostEqual(this.defaultVal, 0) && MathFunctions.almostEqual(other.defaultVal, 0)) {
			return rightProductOnZeroDefVal(this, other);
		} else {
			Map<Integer, Map<Integer, Double>> rowColMap = this.getRowColMap();
			Map<Integer, Map<Integer, Double>> colRowMap = other.getColRowMap();
			SparseMatrix ret = new SparseMatrix(this.nRows, other.nCols,
					this.nCols * this.defaultVal * other.defaultVal);
			for (int row = 0; row < this.nRows; row++) {
				if (rowColMap.containsKey(row)) { // The row of this matrix is
													// not
													// default
					Map<Integer, Double> colMap = rowColMap.get(row);
					Double valOfMulDefCol = null;
					for (int col = 0; col < other.nCols; col++) {
						if (colRowMap.containsKey(col)) { // The column of other
															// matrix is
															// not default
							Map<Integer, Double> rowMap = colRowMap.get(col);
							// Compute the inner product of the two vectors
							Set<Integer> allIndices = new HashSet<Integer>(colMap.keySet());
							allIndices.addAll(rowMap.keySet());
							double val = 0;
							for (int index : allIndices) {
								val += colMap.get(index) * rowMap.get(index);
							}
							val += (this.nCols * -allIndices.size()) * (this.defaultVal * other.defaultVal);
							ret.set(row, col, val);
						} else { // The column of other matrix is default
							if (valOfMulDefCol == null) {
								double val = 0;
								for (int index : colMap.keySet()) {
									val += colMap.get(index) * other.defaultVal;
								}
								val += (this.nCols - colMap.keySet().size()) * this.defaultVal * other.defaultVal;
								ret.set(row, col, val);
								valOfMulDefCol = val;
							} else {
								ret.set(row, col, valOfMulDefCol);
							}
						}
					}
				} else { // The row of this matrix is default
					for (int col = 0; col < other.nCols; col++) {
						if (colRowMap.containsKey(col)) { // The column of other
															// matrix is
															// not default
							Map<Integer, Double> rowMap = colRowMap.get(col);
							// Compute the inner product of the two vectors
							double val = 0;
							for (int index : rowMap.keySet()) {
								val += rowMap.get(index) * this.defaultVal;
							}
							val += (this.nCols * -rowMap.keySet().size()) * (this.defaultVal * other.defaultVal);
							ret.set(row, col, val);
						} else {
							// The column of other matrix is default, do nothing
						}
					}
				}
			}
			return ret;
		}
	}

	/**
	 * Fast compute the product of two matrices if their default values are both
	 * zero
	 * 
	 * @param m1
	 * @param m2
	 * @return
	 */
	private static SparseMatrix rightProductOnZeroDefVal(SparseMatrix m1, SparseMatrix m2) {
		if (m1.nCols != m2.nRows) {
			throw new IllegalArgumentException("This matrix's column is not equal to the other matrix's row");
		}
		Map<Integer, Map<Integer, Double>> rowColMap = m1.getRowColMap();
		Map<Integer, Map<Integer, Double>> colRowMap = m2.getColRowMap();
		SparseMatrix ret = new SparseMatrix(m1.nRows, m2.nCols);
		for (int row : rowColMap.keySet()) {
			// Find the row vector of this matrix
			Map<Integer, Double> colMap = rowColMap.get(row);
			for (int col : colRowMap.keySet()) {
				// Find the column vector of the other matrix
				Map<Integer, Double> rowMap = colRowMap.get(col);
				// Compute the inner product of the two vectors
				double val = 0;
				for (int index : colMap.keySet()) {
					if (rowMap.containsKey(index)) {
						val += colMap.get(index) * rowMap.get(index);
					}
				}
				ret.set(row, col, val);
			}
		}
		return ret;
	}

	/**
	 * Compute the product of a column vector and a row vector
	 * 
	 * @param columnVector
	 * @param rowVector
	 * @return
	 */
	public static SparseMatrix columnRowVectorProduct(SparseVector columnVector, SparseVector rowVector) {
		SparseMatrix columnMatrix = new SparseMatrix(columnVector, true);
		SparseMatrix rowMatrix = new SparseMatrix(rowVector, false);
		return columnMatrix.rightProduct(rowMatrix);
	}

	public SparseVector leftRowVectorProduct(SparseVector v) {
		SparseMatrix leftRow = new SparseMatrix(v, false);
		SparseMatrix retMatrix = leftRow.rightProduct(this);
		return retMatrix.toVector();
	}

	public SparseVector toVector() {
		if (nCols == 1) {
			SparseVector ret = new SparseVector(nRows, this.defaultVal);
			for (Pair<Integer, Integer> key : values.keySet()) {
				ret.set(key.getLeft(), values.get(key));
			}
			return ret;
		} else if (nRows == 1) {
			SparseVector ret = new SparseVector(nCols, this.defaultVal);
			for (Pair<Integer, Integer> key : values.keySet()) {
				ret.set(key.getRight(), values.get(key));
			}
			return ret;
		} else {
			throw new IllegalStateException("This matrix is not a single column or row matrix");
		}
	}

	public static SparseMatrix createIdentityMatrix(int dim) {
		SparseMatrix matrix = new SparseMatrix(dim, dim);
		for (int i = 0; i < dim; i++) {
			matrix.set(i, i, 1.0);
		}
		return matrix;
	}

	public double[][] toDoubleArray() {
		double[][] ret = new double[nRows][nCols];
		for (int i = 0; i < nRows; i++) {
			for (int j = 0; j < nCols; j++)
				ret[i][j] = 0.0;
		}
		Iterator<Pair<Integer, Integer>> it = values.keySet().iterator();
		while (it.hasNext()) {
			Pair<Integer, Integer> pair = it.next();
			ret[pair.getLeft()][pair.getRight()] = values.get(pair);
		}

		return ret;

	}

	@Override
	public String toString() {
		return values.toString();
	}

}
