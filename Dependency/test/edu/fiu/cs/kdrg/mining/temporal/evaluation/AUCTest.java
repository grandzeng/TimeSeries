package edu.fiu.cs.kdrg.mining.temporal.evaluation;

import static org.junit.Assert.*;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

public class AUCTest {

	@Test
	public void testCompute() {
		SimpleMatrix expected = new SimpleMatrix(10, 1);
		SimpleMatrix actual = new SimpleMatrix(10, 1);
		expected.set(0, 1);
		expected.set(1, 1);
		expected.set(2, 0);
		expected.set(3, 0);
		expected.set(4, 1);
		expected.set(5, 1);
		expected.set(6, 0);
		expected.set(7, 0);
		expected.set(8, 0);
		expected.set(9, 1);
		
		actual.set(0, 0.3);
		actual.set(1, -0.7);
		actual.set(2, 0.001);
		actual.set(3, 0.002);
		actual.set(4, 0.4);
		actual.set(5, 0.6);
		actual.set(6, 0.03);
		actual.set(7, 0.05);
		actual.set(8, -0.99);
		actual.set(9, -0.9);
		assertEquals(AUC.compute(expected, actual), 0.8,0.0);
//		System.out.println(AUC.compute(expected, actual));
	}

}
