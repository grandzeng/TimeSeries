package edu.fiu.cs.kdrg.mining.temporal.alg;

import static org.junit.Assert.*;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordIterator;

public class BayesianLinearRegressionApproximationTest {

BayesianLinearRegressionApproximation blr = new BayesianLinearRegressionApproximation(1,1);
	
	RecordIterator iterator = new RecordIterator("data/records.data");

	@Test
	public void test() {
		int lag = blr.getLag();
		int dim = blr.getDimension();
		TimeFrame tf = new TimeFrame(lag, dim);
		while(iterator.hasNext()){
			SimpleMatrix temp = iterator.next();
			
			if(lag>0){
				tf.shift(temp);
				lag --;
				continue;
			}else if(lag == 0){
				blr.setTimeFrame(tf);
				lag --;
			}
			
			blr.onlineTrain(temp);
			System.out.println("mean");
			blr.getPosteriorMean(0).print();
			System.out.println("variance");
			blr.getPosteriorVariance(0).print();
		}
	}

}
