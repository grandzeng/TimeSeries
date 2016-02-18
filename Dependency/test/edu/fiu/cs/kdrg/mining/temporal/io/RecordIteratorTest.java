package edu.fiu.cs.kdrg.mining.temporal.io;

import static org.junit.Assert.*;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

public class RecordIteratorTest {
	RecordIterator iterator = new RecordIterator("records.data");

	@Test
	public void testHasNext() {
		assertTrue(iterator.hasNext());
	}

	@Test
	public void testNext() {
		iterator.hasNext();
		SimpleMatrix matrix = iterator.next();
		assertNotNull(matrix);
		matrix.print();
	}

	@Test
	public void testClose() {
		try {
			iterator.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		assertFalse(iterator.isAvailable());
	}

}
