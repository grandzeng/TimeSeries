package edu.fiu.cs.kdrg.mining.temporal.io;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

import org.ejml.simple.SimpleMatrix;

import edu.fiu.cs.kdrg.mining.temporal.core.SparseVector;

public class RecordWriter {

	private PrintWriter pw;

	public RecordWriter(String fileName) throws FileNotFoundException {
		pw = new PrintWriter(new FileOutputStream(fileName));
	}

	public void outLine(SimpleMatrix[] vecs) {

		for (int i = 0; i < vecs.length; i++) {
			for (int j = 0; j < vecs[i].getNumElements(); j++)
				if (i == 0 && j == 0)
					pw.print(vecs[i].get(j));
				else
					pw.print("," + vecs[i].get(j));
		}
		pw.println();
	}
	
	public void outLine(SparseVector[] vecs){
		for(int i=0;i<vecs.length;i++){
			for(int j=0;j<vecs[i].getDimension();j++){
				if (i == 0 && j == 0)
					pw.print(vecs[i].get(j));
				else
					pw.print("," + vecs[i].get(j));
			}
		}
		pw.println();
	}

	public void close() {
		try {
			pw.close();
		} catch (Exception ex) {
		}
	}

}
