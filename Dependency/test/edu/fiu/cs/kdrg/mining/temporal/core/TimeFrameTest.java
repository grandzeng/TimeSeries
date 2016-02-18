package edu.fiu.cs.kdrg.mining.temporal.core;

import static org.junit.Assert.*;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

public class TimeFrameTest {
	TimeFrame tf = null;

	@Test
	public void testCopy() {
		tf = new TimeFrame(3, 2);
		TimeFrame copy = tf.copy();
		assertEquals("Copy not equal", tf, copy);
		// tf.getVector().print();
	}

	@Test
	public void testGetValue() {
		tf = new TimeFrame(3, 2);
		tf.setValue(0, 0, 1);
		assertEquals("Not equal", 1, tf.getValue(0, 0), 0.0);
	}

	@Test
	public void testSetValue() {
		tf = new TimeFrame(3, 2);
		tf.setValue(0, 1, 5);
		assertEquals("Not equal", 5, tf.getValue(0, 1), 0.0);
	}

	@Test
	public void testShiftSimpleMatrix() {
		tf = new TimeFrame(3, 2);
		tf.setValue(0, 0, 1);
		tf.setValue(0, 1, 2);
		tf.setValue(1, 0, 3);
		tf.setValue(1, 1, 4);
		tf.setValue(2, 0, 5);
		tf.setValue(2, 1, 6);
		double[] values = { 9.0, 10.0 };
		SimpleMatrix matrix = new SimpleMatrix(2, 1, true, values);
		tf.shift(matrix);
		// tf.getVector().print();

		TimeFrame expected = new TimeFrame(3, 2);
		expected.setValue(0, 0, 9.0);
		expected.setValue(0, 1, 10.0);
		expected.setValue(1, 0, 1);
		expected.setValue(1, 1, 2);
		expected.setValue(2, 0, 3);
		expected.setValue(2, 1, 4);
		assertEquals("Two TimeFrame not equal", expected, tf);
	}

	@Test
	public void testShiftDoubleArray() {
		tf = new TimeFrame(3, 2);
		tf.setValue(0, 0, 1);
		tf.setValue(0, 1, 2);
		tf.setValue(1, 0, 3);
		tf.setValue(1, 1, 4);
		tf.setValue(2, 0, 5);
		tf.setValue(2, 1, 6);
		double[] values = { 7.0, 8.0 };
		tf.shift(values);
		// tf.getVector().print();

		TimeFrame expected = new TimeFrame(3, 2);
		expected.setValue(0, 0, 7.0);
		expected.setValue(0, 1, 8.0);
		expected.setValue(1, 0, 1);
		expected.setValue(1, 1, 2);
		expected.setValue(2, 0, 3);
		expected.setValue(2, 1, 4);
		assertEquals("Two TimeFrame not equal", expected, tf);
	}

	@Test
	public void testGetVector() {
		tf = new TimeFrame(3, 2);
		tf.setValue(0, 0, 1);
		tf.setValue(0, 1, 2);
		tf.setValue(1, 0, 3);
		tf.setValue(1, 1, 4);
		tf.setValue(2, 0, 5);
		tf.setValue(2, 1, 6);
		
		SimpleMatrix vector = tf.getVector();
		assertEquals(vector.getNumElements(), 6);
		assertEquals(tf.getValue(1, 0),vector.get(2),0.0);
	}

}
