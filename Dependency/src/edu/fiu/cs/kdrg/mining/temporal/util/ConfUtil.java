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
	
	private static Properties prop;
	
	static{
		try {
			prop = loadConf();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static Properties loadConf() throws IOException{
		Properties prop = new Properties();
		InputStream pIn = new FileInputStream(confFile);
		prop.load(pIn);
		pIn.close();
		return prop;
	}
	
	public static String getMetaFile(){
		return prop.getProperty("metaFile");
	}
	
	public static String getDataFile(){
		return prop.getProperty("dataFile");
	}
	
	public static String[] getEvaluationFiles(){
		String temp = prop.getProperty("evaluationFiles");
		return temp.split(",");
	}
	
	public static String getPerformanceFile(){
		return prop.getProperty("performanceFile");
	}
	
	public static int getDimension(){
		return Integer.parseInt(prop.getProperty("dimension"));
	}
	
	public static int getLag(){
		return Integer.parseInt(prop.getProperty("lag"));
	}
	
	public static int getLength(){
		return Integer.parseInt(prop.getProperty("length"));
	}
	
}
