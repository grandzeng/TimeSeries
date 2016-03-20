package edu.fiu.cs.kdrg.mining.temporal.alg;

import static org.junit.Assert.*;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import edu.fiu.cs.kdrg.mining.temporal.core.TimeFrame;
import edu.fiu.cs.kdrg.mining.temporal.io.RecordIterator;

public class BayesianLassoTest {

	int dim =1; 
	int numOfParticles =1000;
	int lag = 1;
	double lambda = 1.0;
	double lambda2 = 1.0;
	
	BayesianLasso bayesianLasso = new BayesianLasso(dim,numOfParticles,lag,lambda);
	
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
				bayesianLasso.setTimeFrame(tf);
				lag --;
			}
			
			bayesianLasso.onlineTrain(temp);
			System.out.println("mean");
			bayesianLasso.getPosteriorMean(0).print();
			System.out.println("variance");
			bayesianLasso.getPosteriorVariance(0).print();
		}
	}

}
