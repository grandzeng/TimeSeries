package edu.fiu.cs.kdrg.mining.temporal.alg;

import static org.junit.Assert.*;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordIterator;

public class TimeVaryingLinearRegressionTest {
	int dim =1; 
	int numOfParticles =1;
	int lag = 1;
	double lambda1 = 1.0;
	double lambda2 = 1.0;
	
	TimeVaryingLinearRegression tvlr = new TimeVaryingLinearRegression(dim,numOfParticles,lag,lambda1,lambda2);
	
	RecordIterator iterator = new RecordIterator("data/records.data");

	@Test
	public void testTrain() {
		TimeFrame tf = new TimeFrame(dim,lag);
		while(iterator.hasNext()){
			SimpleMatrix temp = iterator.next();
			
			if(lag>0){
				tf.shift(temp);
				lag --;
				continue;
			}else if(lag == 0){
				tvlr.setTimeFrame(tf);
				lag --;
			}
			
			tvlr.onlineTrain(temp);
			System.out.println("mean");
			tvlr.getPosteriorMean(0).print();
			System.out.println("variance");
			tvlr.getPosteriorVariance(0).print();
		}
	}

	

}
