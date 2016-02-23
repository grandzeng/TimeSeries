package edu.fiu.cs.kdrg.mining.temporal.util;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * @author Chunqiu Zeng
 * @date 22nd Feb, 2016
 *
 */
public class ConfUtil {
	
	public static String confFile = "conf/sim.conf";
	
	public static Properties loadConf() throws IOException{
		Properties prop = new Properties();
		InputStream pIn = new FileInputStream(confFile);
		prop.load(pIn);
		pIn.close();
		return prop;
	}
	
	
}
