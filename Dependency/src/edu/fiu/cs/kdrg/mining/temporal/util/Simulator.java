package edu.fiu.cs.kdrg.mining.temporal.util;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Properties;
import java.util.Random;

import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;

/**
 * This class is used to simulate the multivariate time series data.
 * 
 * @author Chunqiu Zeng
 * @date 22nd Feb, 2016
 *
 */
public class Simulator {

	/**
	 * The number of time series
	 */
	private int dimension;

	/**
	 * The time lag
	 */
	private int lag;

	/**
	 * The number of instances
	 */
	private int length;

	/**
	 * An instance of TimeFrame
	 */
	private TimeFrame tf;

	/**
	 * The dependency matrix
	 */
	private SimpleMatrix dependency;

	/**
	 * The sparsity ratio.
	 */
	private double sparsity = 0.3;

	/**
	 * The mean of the coefficients
	 */
	private double mean = 10.0;

	/**
	 * The standard deviation
	 */
	private double std = 3.0;

	/**
	 * The random number generator.
	 */
	private Random rand = new Random(System.nanoTime());

	/**
	 * The writer for storing the simulation data.
	 */
	private Writer writer;

	/**
	 * The constructor
	 * 
	 * @param dim
	 * @param len
	 * @param lag
	 * @param writer
	 */
	public Simulator(int dim, int len, int lag, Writer writer) {
		dimension = dim;
		length = len;
		this.lag = lag;
		tf = new TimeFrame(lag, dimension);
		dependency = new SimpleMatrix(dimension, lag * dimension);
		this.writer = writer;

	}

	/**
	 * Meta data generation
	 * 
	 * @param timestamp
	 */
	public void metadata(int timestamp) {
		// constant dependency
		if (timestamp != 0)
			return;
		double p;
		for (int i = 0; i < dependency.numRows(); i++) {
			for (int j = 0; j < dependency.numCols(); j++) {
				p = rand.nextDouble();
				if (p < sparsity) {
					dependency.set(i, j, rand.nextGaussian() * std + mean);
				} else {
					dependency.set(i, j, 0.0);
				}
			}
		}
		writer.writeMetadata(timestamp, dependency);
	}

	/**
	 * The simulation function
	 */
	public void simulate() {
		for (int r = 0; r < length; r++) {
			metadata(r);
			SimpleMatrix past = tf.getVector();
			SimpleMatrix current = dependency.mult(past);
			SimpleMatrix noise = new SimpleMatrix(dimension, 1);
			for (int i = 0; i < dimension; i++) {
				noise.set(i, rand.nextGaussian());
			}
			current = current.plus(noise);
			tf.shift(current);
			writer.writeInstance(current);
		}
	}

	/**
	 * This inner class is used to store the simulation data.
	 *
	 */
	public static class Writer {
		private PrintWriter dataWriter;
		private PrintWriter metadataWriter;

		public Writer(String dataFile, String metaFile) {
			try {
				dataWriter = new PrintWriter(new FileOutputStream(dataFile));
				metadataWriter = new PrintWriter(new FileOutputStream(metaFile));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}

		public void close() {
			try {
				dataWriter.close();
				metadataWriter.close();
			} catch (Exception ignore) {
			}
		}

		public void writeInstance(SimpleMatrix instance) {
			outLine(instance, dataWriter);
		}

		public void writeMetadata(int timestamp, SimpleMatrix meta) {
			metadataWriter.print(timestamp + ",");
			outLine(meta, metadataWriter);

		}

		private void outLine(SimpleMatrix vec, PrintWriter pw) {
			for (int i = 0; i < vec.numRows(); i++) {
				for (int j = 0; j < vec.numCols(); j++) {
					if (i == 0 && j == 0)
						pw.print(vec.get(i, j));
					else
						pw.print("," + vec.get(i, j));
				}
			}
			pw.println();
		}
	}

	/**
	 * The main Entry
	 * @param args
	 * @throws IOException
	 */
	public static void main(String args[]) throws IOException {
		Properties prop = ConfUtil.loadConf();
		int dimension = Integer.parseInt(prop.getProperty("dimension"));
		int lag = Integer.parseInt(prop.getProperty("lag"));
		int length = Integer.parseInt(prop.getProperty("length"));
		String metaFile = prop.getProperty("metaFile");
		String dataFile = prop.getProperty("dataFile");
		Writer writer = new Writer(dataFile, metaFile);
		Simulator sim = new Simulator(dimension, length, lag, writer);
		sim.simulate();
		writer.close();
	}

}
