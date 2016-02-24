package edu.fiu.cs.kdrg.mining.temporal.util;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.Properties;

import org.junit.Test;

public class ConfUtilTest {

	@Test
	public void testLoadConf() {
		try {
			Properties prob = ConfUtil.loadConf();
			System.out.println(prob.keySet());
			System.out.println(prob.values());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
