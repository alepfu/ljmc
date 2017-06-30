package cs.cophy2.ljmc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;

public class AverageRDF {

	public static void main(String[] args) {
		
		int numBins = 100;
		
		try {
			File dir = new File("/home/alepfu/Desktop/LJMC");
			File[] files = dir.listFiles(new FilenameFilter() {
			    public boolean accept(File dir, String name) {
			        return name.startsWith("rdf");
			    }
			});

			double[] hist = new double[100];
			
			for (File f : files) {
				BufferedReader br = new BufferedReader(new FileReader(f));
				String line;
				int i = 0;
				while((line = br.readLine()) != null) {
					String[] split = line.split("\\s+");
					hist[i] += Double.parseDouble(split[1]) / files.length;
					++i;
				}
			}
			
			double binSize = 0.5 * 7.9370 / numBins;
			
			try {
				BufferedWriter w = new BufferedWriter(new FileWriter("/home/alepfu/Desktop/LJMC/avg_rdf.csv"));
				for (int i = 0; i < numBins; i++)
					w.write((i * binSize) + " " + hist[i] + "\n");
				w.flush();
				w.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

}
