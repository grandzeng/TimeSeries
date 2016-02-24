package edu.fiu.cs.kdrg.mining.temporal.evaluation;

import static org.junit.Assert.*;

import java.util.Properties;

import org.junit.Test;

import edu.fiu.cs.kdrg.mining.temporal.util.ConfUtil;

public class EvaluationTest {

	@Test
	public void testEvaluate() {
		Evaluation evaluation = new Evaluation(ConfUtil.getMetaFile(), ConfUtil.getEvaluationFiles(),
				ConfUtil.getPerformanceFile());
	}

}
