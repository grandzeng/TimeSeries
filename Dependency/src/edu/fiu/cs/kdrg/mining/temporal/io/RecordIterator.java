package edu.fiu.cs.kdrg.mining.temporal.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;

import org.ejml.simple.SimpleMatrix;

public class RecordIterator implements Iterator<SimpleMatrix> {

	private String dataFileName;

	private BufferedReader reader = null;

	private boolean isClosed = true;

	private SimpleMatrix current = null;

	public RecordIterator(String fileName) {
		this.dataFileName = fileName;
	}

	private void initIterator() throws FileNotFoundException {
		InputStream fileStream = new FileInputStream(dataFileName);
		reader = new BufferedReader(new InputStreamReader(fileStream));
		isClosed = false;
	}

	private SimpleMatrix readRecord() throws IOException {
		if (reader == null) {
			initIterator();
		}
		if (isClosed)
			return null;
		String line = reader.readLine();
		if (line == null) {
			close();
			return null;
		} else {
			return parseLine(line);
		}
	}

	private SimpleMatrix parseLine(String line) {
		line = line.trim();
		String[] array = line.split(",");
		SimpleMatrix vector = new SimpleMatrix(array.length, 1);
		for (int i = 0; i < array.length; i++) {
			vector.set(i, Double.parseDouble(array[i].trim()));
		}
		return vector;
	}

	public void close() throws IOException {
		if (reader != null && !isClosed) {
			reader.close();
			isClosed = true;
		}
	}

	@Override
	public boolean hasNext() {
		try {
			current = readRecord();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return current != null;
	}

	@Override
	public SimpleMatrix next() {
		return current;
	}

	public boolean isAvailable() {
		return !isClosed;
	}

}
