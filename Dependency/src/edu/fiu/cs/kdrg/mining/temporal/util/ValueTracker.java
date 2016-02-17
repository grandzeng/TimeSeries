package edu.fiu.cs.kdrg.mining.temporal.util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * For debugging
 *
 * @author Liang Tang
 * @date May 26, 2014 3:47:47 PM
 */
public class ValueTracker {
	
	SortedMap<String, List<String>> attrValues = new TreeMap<String, List<String>>();
	
	public ValueTracker() {
		
	}
	
	public void add(String attrName, double val) {
		add(attrName, val+"");
	}
	
	public void add(String attrName, String val) {
		List<String> valueList = attrValues.get(attrName);
		if (valueList == null) {
			valueList = new ArrayList<String>();
			attrValues.put(attrName, valueList);
		}
		valueList.add(val);
	}
	
	public String toString(String attrName) {
		StringBuffer buf = new StringBuffer();
		buf.append(attrName + " : ");
		List<String> valueList = attrValues.get(attrName);
		if (valueList == null) {
			buf.append("null");
		}
		else {
			buf.append("["+ valueList.size()+"] "+valueList.toString());
		}
		return buf.toString();
	}
	
	public void print(String attrName) {
		System.out.println(toString(attrName));
	}
	
	public void save(String fileName) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			save(writer);
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	public void save(Writer writer) throws IOException {
		for (String attrName :  attrValues.keySet()) {
			writer.write(toString(attrName)+"\n");				
		}
	}
	

}
