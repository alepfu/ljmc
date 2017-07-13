package cs.cophy2.ljmc;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Verify results from Monte Carlo Lennard-Jones simulation
 * with restults from Johnson et al.
 * 
 * @author Alexander Pfundner
 * 
 */
public class Verification {

	public static String workDir = "/home/alepfu/Desktop/LJMC";
	
	/*public static void main(String[] args) {
		
		//Parameters set of Johnson et al.
		double[][] paramTable = {{0.5, 5.0, 4.0}, //density, temp, cutoff
				 {0.1, 2.0, 5.0},
				 {0.2, 2.0, 5.0},
				 {0.4, 2.0, 5.0},
				 {0.5, 2.0, 5.0},
				 {0.6, 2.0, 4.71},
				 {0.7, 2.0, 4.47},
				 {0.8, 2.0, 4.27},
				 {0.9, 2.0, 4.11},
				 {0.1, 1.2, 4.0},
				 {0.1, 1.15, 4.0},
				 {0.05, 1.0, 5.0},
				 {0.6, 1.0, 4.71},
				 {0.7, 1.0, 4.47},
				 {0.8, 1.0, 4.27},
				 {0.8, 1.0, 4.0},
				 {0.9, 1.0, 4.11}};
		
		double eps = 1.0;
		int n = 200;
		
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter(workDir + "/verification.txt"));
			w.write("density temp cutoff pressure energy\n");
			
			for (int i = 0; i < paramTable.length; i++) {
				
				LJMCSimulation sim = new LJMCSimulation(paramTable[i][0], paramTable[i][1], paramTable[i][2], eps, n);
				sim.run();

				w.write(sim.density + " " + sim.temp + " " + sim.cutoff + " "
						+ sim.pressure + " " + sim.avgEnergy + "\n");
			}
			
			w.flush();
			w.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}

	}*/
	
}
