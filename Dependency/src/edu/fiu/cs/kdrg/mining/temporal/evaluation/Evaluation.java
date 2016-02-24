package edu.fiu.cs.kdrg.mining.temporal.evaluation;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Iterator;

import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.io.RecordIterator;

/**
 * This class is used to evaluate the performance of different algorithms with
 * respect to the ground truth.
 * 
 * @author Chunqiu Zeng
 * @date 23rd Feb, 2016
 *
 */
public class Evaluation implements Iterator<SimpleMatrix[]> {

	/**
	 * The file name of ground truth .
	 */
	private String metaFile;

	private RecordIterator metaIterator;

	private String[] evaluationFiles;

	private RecordIterator[] evaluationIterators;

	private String performanceFile;

	private PrintWriter performanceWriter;

	public Evaluation(String metaFile, String[] evaluationFiles, String performanceFile) {
		this.metaFile = metaFile;
		this.evaluationFiles = evaluationFiles;
		this.performanceFile = performanceFile;
		try {
			init();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public void init() throws FileNotFoundException {
		metaIterator = new RecordIterator(metaFile);
		performanceWriter = new PrintWriter(new FileOutputStream(performanceFile));
		evaluationIterators = new RecordIterator[evaluationFiles.length];
		for (int i = 0; i < evaluationFiles.length; i++) {
			evaluationIterators[i] = new RecordIterator(evaluationFiles[i]);
			if (i == 0)
				performanceWriter.print(evaluationFiles[i]);
			else
				performanceWriter.print("," + evaluationFiles[i]);
		}
		performanceWriter.println();
	}

	public void evaluate() {
		while (hasNext()) {
			SimpleMatrix[] matrices = next();
			for (int i = 0; i < matrices.length - 1; i++) {
				if (i == 0)
					performanceWriter.print(AUC.compute(matrices[0], matrices[i + 1]));
				else
					performanceWriter.print("," + AUC.compute(matrices[0], matrices[i + 1]));
			}
			
			performanceWriter.println();
		}

		performanceWriter.close();
	}

	@Override
	public boolean hasNext() {
		if (!metaIterator.hasNext())
			return false;
		for (RecordIterator iter : evaluationIterators) {
			if (!iter.hasNext())
				return false;
		}
		return true;
	}

	@Override
	public SimpleMatrix[] next() {
		SimpleMatrix[] matrices = new SimpleMatrix[evaluationFiles.length + 1];
		matrices[0] = metaIterator.next();
		for (int i = 0; i < evaluationFiles.length; i++) {
			matrices[i + 1] = evaluationIterators[i].next();
		}
		return matrices;
	}

}
