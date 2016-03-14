package edu.fiu.cs.kdrg.mining.temporal.util;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

public class SimpleMatrixTest {

	@Test
	public void testSimpleMatrix(){
		double [][]testmat = {{2.0,0},{0,1}};
		SimpleMatrix mat = new SimpleMatrix(testmat);
		mat.print();
		mat.invert().print();
		mat.print();
	}
}
