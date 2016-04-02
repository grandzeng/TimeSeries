package edu.fiu.cs.kdrg.mining.temporal.main;

import java.io.FileNotFoundException;

import edu.fiu.cs.kdrg.mining.temporal.alg.BayesianLasso;
import edu.fiu.cs.kdrg.mining.temporal.alg.BayesianLinearRegression;
import edu.fiu.cs.kdrg.mining.temporal.alg.TimeVaryingLinearRegression;
import edu.fiu.cs.kdrg.mining.temporal.evaluation.Evaluation;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordIterator;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordWriter;
import edu.fiu.cs.kdrg.mining.temporal.util.ConfUtil;

public class Experiment {

	public static void runBayesianLinearRegression(int index) throws FileNotFoundException {
		BayesianLinearRegression blr = new BayesianLinearRegression(ConfUtil.getLag(), ConfUtil.getDimension());
		RecordWriter rw = new RecordWriter(ConfUtil.getEvaluationFiles()[index]);
		blr.setOutput(rw);
		RecordIterator recIter = new RecordIterator(ConfUtil.getDataFile());
		while (recIter.hasNext()) {
			blr.onlineTrain(recIter.next());
		}
		rw.close();
	}

	public static void runBayesianLasso(int index) throws FileNotFoundException {
		int numOfParticles = 20;
		double lambda = 1.0;
		BayesianLasso bl = new BayesianLasso(ConfUtil.getDimension(), numOfParticles, ConfUtil.getLag(), lambda);
		RecordWriter rw = new RecordWriter(ConfUtil.getEvaluationFiles()[index]);
		bl.setOutput(rw);
		RecordIterator recIter = new RecordIterator(ConfUtil.getDataFile());
		while (recIter.hasNext()) {
			bl.onlineTrain(recIter.next());
		}
		rw.close();
	}
	
	public static void runTimeVaryingBayesianLasso(int index) throws FileNotFoundException {
		int numOfParticles = 20;
		double lambda1 = 1.0;
		double lambda2 = 1.0;
		TimeVaryingLinearRegression bl = new TimeVaryingLinearRegression(ConfUtil.getDimension(), numOfParticles, ConfUtil.getLag(), lambda1, lambda2);
		RecordWriter rw = new RecordWriter(ConfUtil.getEvaluationFiles()[index]);
		bl.setOutput(rw);
		RecordIterator recIter = new RecordIterator(ConfUtil.getDataFile());
		while (recIter.hasNext()) {
			bl.onlineTrain(recIter.next());
		}
		rw.close();
	}

	public static void evaluate() {
		Evaluation eval = new Evaluation(ConfUtil.getMetaFile(), ConfUtil.getEvaluationFiles(),
				ConfUtil.getPerformanceFile());
		eval.evaluate();
	}

	public static void main(String[] args) throws FileNotFoundException {
		System.out.println("start");
		long start = System.currentTimeMillis();
		long previous = start;
		long current;
		runBayesianLinearRegression(0);
		current = System.currentTimeMillis();
		System.out.println("baysianLinearRegression=" +(current - previous));
//		previous = current;
//		runBayesianLasso(1);
//		current = System.currentTimeMillis();
//		System.out.println("bayesianLasso=" +(current - previous));
//		previous = current;
//		runTimeVaryingBayesianLasso(2);
//		current = System.currentTimeMillis();
//		System.out.println("timevaryingLasso=" +(current - previous));
//		previous = current;
//		evaluate();
//		current = System.currentTimeMillis();
//		System.out.println("evaluation=" +(current - previous));
//		previous = current;
	}

}
