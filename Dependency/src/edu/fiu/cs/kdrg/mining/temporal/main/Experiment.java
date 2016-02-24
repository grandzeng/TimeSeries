package edu.fiu.cs.kdrg.mining.temporal.main;

import java.io.FileNotFoundException;

import edu.fiu.cs.kdrg.mining.temporal.alg.BayesianLinearRegression;
import edu.fiu.cs.kdrg.mining.temporal.evaluation.Evaluation;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordIterator;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordWriter;
import edu.fiu.cs.kdrg.mining.temporal.util.ConfUtil;

public class Experiment {

	public static void runBayesianLinearRegression() throws FileNotFoundException {
		BayesianLinearRegression blr = new BayesianLinearRegression(ConfUtil.getLag(), ConfUtil.getDimension());
		RecordWriter rw = new RecordWriter(ConfUtil.getEvaluationFiles()[0]);
		blr.setOutput(rw);
		RecordIterator recIter = new RecordIterator(ConfUtil.getDataFile());
		while (recIter.hasNext()) {
			blr.onlineTrain(recIter.next());
		}
		rw.close();
		
	}
	
	public static void evaluate(){
		Evaluation eval = new Evaluation(ConfUtil.getMetaFile(),ConfUtil.getEvaluationFiles(),ConfUtil.getPerformanceFile());
		eval.evaluate();
	}

	public static void main(String[] args) throws FileNotFoundException {
		runBayesianLinearRegression();
		evaluate();
	}

}
