package cs.cophy2.ljmc;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;


public class Simulation {

	public static int numDim;
	public static int numParticles;			
	public static double boxLength;			
	public static List<double[]> particles;			
	public static double temp;				
	public static int numEqSteps;
	public static int numSampSteps;
	public static int numTotalSteps;
	public static double epsilon;			
	public static double sigma;				
    public static double cutoff;				
    public static double radius;			
	public static double density;
	public static double beta;
	public static double boxVolume;
	public static double displ;
	public static double cutoffFactor;
	
	public static int printProgressFreq;
	public static int writeTrjFreq;
	
	public static long seed = 123;
	public static Random rand = new Random(seed);
	public static NumberFormat nf = DecimalFormat.getInstance();
	public static String outDir = "/home/alepfu/Desktop/LJMC";
	public static String outSuffix = "";
	public static BufferedWriter trjWriter;
	
	public static void main(String[] args) {
		
		//verifyWithJohsonEtAl();
		
		//RDF runs
		/*run(0.5, 0.5, 2.5);
		run(0.5, 1.0, 2.5);
		run(0.5, 2.0, 2.5);
		run(0.5, 5.0, 2.5);*/
		
		//Trj runs
		//run(0.2, 1.0, 2.5);
		
		//density temp cutofffactor
		
		//Email
		run(0.2, 1.0, 2.5);   //TODO work in progress
	}
	
	/**
	 * Verify with results from Johnson et al.
	 */
	public static void verifyWithJohsonEtAl() {
		double[][] paramTable = {{0.5, 5.0, 4.0}, //density, temp, cutoff/sigma
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
		try {
			BufferedWriter resultsWriter = new BufferedWriter(new FileWriter(outDir + "/results.txt"));
			
			for (int i = 0; i < paramTable.length; i++) {
				String result = run(paramTable[i][0], paramTable[i][1], paramTable[i][2]);
				log(result);
				resultsWriter.write(result + "\n");
			}
			
			resultsWriter.flush();
			resultsWriter.close();
			
			} catch (Exception e) {
				e.printStackTrace();
			}
	}
	
	/**
	 * Monte Carlo simulation of Lennard-Jones Particles.
	 */
	public static String run(double param1, double param2, double param3) {
		
		nf.setMaximumFractionDigits(4);
		nf.setMinimumFractionDigits(4);
		trjWriter = null;

		//Input parameters
		density = param1;
		temp = param2;
		cutoffFactor = param3;
		
		//Simulation parameters
		numDim = 3;
		numParticles = 100;					//TODO work in progress	
		epsilon = 5.0;						//TODO work in progress		
		sigma = 1.0;
		numEqSteps = 100000;
		numSampSteps = 100000;
		displ = 0.5;
		
		//Calculated parameters
		boxVolume = numParticles / density;
		boxLength = Math.pow(boxVolume, 1.0 / numDim);
		cutoff = cutoffFactor * sigma;
		beta = 1.0 / temp;
		radius = 0.5 * sigma;
		numTotalSteps = numEqSteps + numSampSteps;
		
		//Output parameters
		printProgressFreq = 100;
		writeTrjFreq = 100;
		
		//Log parameters
		log("Num. dimensions = " + numDim);
		log("Num. particles = " + numParticles);
		log("Epsilon = " + nf.format(epsilon));
		log("Sigma = " + nf.format(sigma));
		log("Cutoff = " + nf.format(cutoff));
		log("Temperature = " + nf.format(temp));
		log("Density = " + nf.format(density));
		log("Box volume = " + nf.format(boxVolume));
		log("Box length = " + nf.format(boxLength) + "\n");
		
		//Initialize system
		initSystemFrequentlyDistributed();
		
		writeTrajectory();
		
		//Calc initial energy and virial
		double energy = calcEnergyTotal();
		double virial = calcVirialTotal();
		
		//Init sampling quantities
		double energySum = 0.0;
		double virialSum = 0.0;
		int acceptCounter = 0;
		
		//Metropolis algorithm
		for (int step = 0; step < numTotalSteps; step++) {
			
			//Choose a particle at random and calculate its energy and virial
			int particleIndex = rand.nextInt(numParticles);
			double prevParticleEnergy = calcEnergyParticle(particleIndex);
			double prevParticleVirial = calcVirialParticle(particleIndex);
			
			//Move the choosen particle and calculate the new energy and virial
			double[] prevParticle = new double[numDim];
			for (int j = 0; j < prevParticle.length; j++)
				prevParticle[j] = particles.get(particleIndex)[j];
			moveParticle(particleIndex);
			double newParticleEnergy = calcEnergyParticle(particleIndex);
			double newParticleVirial = calcVirialParticle(particleIndex);
			
			//Acceptance rule
			double deltaParticleEnergy = newParticleEnergy - prevParticleEnergy;
			if ((deltaParticleEnergy < 0) || (rand.nextDouble() < Math.exp(-beta * deltaParticleEnergy))) {
				
				//Accept move
				energy += deltaParticleEnergy;
				virial += newParticleVirial - prevParticleVirial;
				++acceptCounter;
				
				if (step % writeTrjFreq == 0)
					writeTrajectory();
				
			} else {
				
				//Restore old configuration
				for (int j = 0; j < prevParticle.length; j++)
					particles.get(particleIndex)[j] = prevParticle[j];
			}
			
			energySum += energy;
			virialSum += virial;
			
			//Reset for sampling
			if (step == numEqSteps) {
				energySum = 0.0;
				virialSum = 0.0;
				acceptCounter = 0;
			}
			
			//Print progress
			if (step % printProgressFreq == 0)
				log("Running " + (step < numEqSteps ? "equilibration" : "sampling")
						+ " [ " + (int)((step * 1.0 / numTotalSteps) * 100) + "% ]");
		}
		
		//Calc and log results
	    double avgEnergy = energySum / numSampSteps / numParticles;
	    double pressureTailCorr = 0.0;
	    if (numDim == 3)
	    	pressureTailCorr = (16.0 / 3.0) * Math.PI * Math.pow(density, 2) * epsilon * Math.pow(sigma, 3) * (  (2.0/3.0) * Math.pow(sigma / cutoff, 9) - Math.pow(sigma / cutoff, 3)  );
	    double pressure = virialSum / 3.0 / numSampSteps / boxVolume + density * temp + pressureTailCorr;
	    double finalAcceptRate = acceptCounter * 1.0 / numSampSteps * 100.0;
		log("Avg. energy = " + nf.format(avgEnergy));
		log("Avg. pressure = " + nf.format(pressure));
		log("Acceptance rate = " + nf.format(finalAcceptRate));
	    
	    //Result line for verification with paper
	    String resultLine = density + "|" + temp + "|" + nf.format(pressure) + "|" + nf.format(avgEnergy) + "|" + cutoffFactor + "|" + nf.format(finalAcceptRate);
	    
	    writeRDFToFile();
	    
		//Plotting
		plotRDF();
//		plot2DSystemPeriodic();

		try {
			trjWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return resultLine; 
	}
	
	/**
	 * Initializes the system with randomly choosen positions.
	 */
	public static void initSystemRandomlyDistributed() {
		particles = new ArrayList<double[]>(numParticles);
		
		//Get particle positions randomly
		for (int i = 0; i < numParticles; i++) {
			double[] pos = new double[numDim];		
			for (int j = 0; j < pos.length; j++)
				pos[j] = rand.nextDouble() * boxLength;
			particles.add(pos);
		}
	}
	
	/**
	 * Initializes the system distributing the particles frequently.
	 */
	public static void initSystemFrequentlyDistributed() {
		particles = new ArrayList<double[]>(numParticles);
		
		if (numDim == 3) {
			//Get number of partitions for the lowest perfect cube
			int numPart = 2;
	        while (Math.pow(numPart, 3) < numParticles)
	            numPart = numPart + 1;
	
	        //Get particle positions by advancing an index
	        int[] index = {0, 0, 0};
			for (int i = 0; i < numParticles; i++) {
				
				double[] pos = new double[numDim];		
				for (int j = 0; j < pos.length; j++)
					pos[j] = (index[j] + radius) * (boxLength / numPart);
				particles.add(pos);
		
				index[0] = index[0] + 1;
				if (index[0] == numPart) {
					index[0] = 0;
					index[1] = index[1] + 1;
					if (index[1] == numPart) {
						index[1] = 0;
						index[2] = index[2] + 1;
					}
				}
			}
		}
		
		if (numDim == 2) {
			//Get number of partitions for the lowest perfect cube
			int numPart = 2;
	        while (Math.pow(numPart, 2) < numParticles)
	            numPart = numPart + 1;
	
	        //Get particle positions by advancing an index
	        int[] index = {0, 0};
			for (int i = 0; i < numParticles; i++) {
				
				double[] pos = new double[numDim];		
				for (int j = 0; j < pos.length; j++)
					pos[j] = (index[j] + radius) * (boxLength / numPart);
				particles.add(pos);
		
				index[0] = index[0] + 1;
				if (index[0] == numPart) {
					index[0] = 0;
					index[1] = index[1] + 1;
				}
			}
		}
	}
	
	/**
	 * Moves a particle at random with applying PBC.
	 */
	public static void moveParticle(int p) {
		
		//x
		particles.get(p)[0] += (rand.nextDouble() - 0.5) * displ;
		if (particles.get(p)[0] >= boxLength)
			particles.get(p)[0] -= boxLength;
		if (particles.get(p)[0] < 0.0)
			particles.get(p)[0] += boxLength;
		
		//y
		particles.get(p)[1] += (rand.nextDouble() - 0.5) * displ;
		if (particles.get(p)[1] >= boxLength)
			particles.get(p)[1] -= boxLength;
		if (particles.get(p)[1] < 0.0)
			particles.get(p)[1] += boxLength;
		
		//z
		if (numDim == 3) {
			particles.get(p)[2] += (rand.nextDouble() - 0.5) * displ;
			if (particles.get(p)[2] >= boxLength)
				particles.get(p)[2] -= boxLength;
			if (particles.get(p)[2] < 0.0)
				particles.get(p)[2] += boxLength;
		}
	}

	/**
	 * Calculates the energy of the system.
	 */
	public static double calcEnergyTotal() {
		double energyTotal = 0.0;
		
		for (int i = 0; i < numParticles; i++)
			energyTotal += calcEnergyParticle(i);
		
		energyTotal = 0.5 * energyTotal;
		
		//Tail correction
		if (numDim == 3) {
			double energyTailCorr = (8.0/3.0) * Math.PI * density * epsilon * Math.pow(sigma, 3) * (  (1.0/3.0) * Math.pow(sigma/cutoff, 9) - Math.pow(sigma/cutoff, 3)  );
			energyTotal += numParticles * energyTailCorr;
		}
		//TODO tail correction for 2D
		
		return energyTotal;
	}

	/**
	 * Calculates the energy of a single particle.
	 */
	public static double calcEnergyParticle(int p) {
		double energyParticle = 0.0;
		
		for (int j = 0; j < numParticles; j++) {
			if (j != p) {
				double dist = calcDist(p, j);
				if (dist <= cutoff)
					energyParticle += calcLJPotential(dist);
			}
		}
		
		return energyParticle;
	}
	
	/**
	 * Calculates the Lennard-Jones potential for a given distance.
	 */
	public static double calcLJPotential(double dist) {
		
		double pot = 4.0 * epsilon * (  Math.pow(sigma / dist, 12) - Math.pow(sigma / dist, 6)  ); 
		
		if (numDim == 3) {
			double energyShiftCorr = 4.0 * epsilon * (  Math.pow(sigma / cutoff, 12) - Math.pow(sigma / cutoff, 6)  ) ;
			pot -= energyShiftCorr;
		}
		//TODO shift correction for 2D
		
		return pot;
	}
	
	/**
	 * Calculates the virial of the system.
	 */
	public static double calcVirialTotal() {
		double virialTotal = 0.0;
		
		for (int i = 0; i < numParticles; i++)
			virialTotal += calcVirialParticle(i);

		return 0.5 * virialTotal;
	}
	
	/**
	 * Calculates the virial of a single particle.
	 */
	public static double calcVirialParticle(int p) {
		double virialParticle = 0.0;
		
		for (int j = 0; j < numParticles; j++) {
			if (j != p) {
				double dist = calcDist(p, j);
				if (dist <= cutoff)
					virialParticle += calcVirial(dist);
			}
		}
		
		return virialParticle;
	}
	
	/**
	 * Calculates the virial for a given distance
	 */
	public static double calcVirial(double dist) {
		
		return 48.0 * epsilon * (  Math.pow(sigma/dist, 12) - 0.5 * Math.pow(sigma/dist, 6)  );
	}
	
	/**
	 * Calculates the distance between two particles,
	 * restricting to the minimum image ("wrapping distances").
	 */
	public static double calcDist(int i, int j) {
		double dist = 0.0;
		
		double[] p1 = particles.get(i);
		double[] p2 = particles.get(j);
		
		if (numDim == 3) {
			double dx = p1[0] - p2[0];
			double dy = p1[1] - p2[1];
			double dz = p1[2] - p2[2];
			
			if (dx >  boxLength/2) dx -= boxLength;
			if (dx < -boxLength/2) dx += boxLength;
			if (dy >  boxLength/2) dy -= boxLength;
			if (dy < -boxLength/2) dy += boxLength;
			if (dz >  boxLength/2) dz -= boxLength;
			if (dz < -boxLength/2) dz += boxLength;
			
			dist = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2));
		}
		
		if (numDim == 2) {
			double dx = p1[0] - p2[0];
			double dy = p1[1] - p2[1];
			
			if (dx >  boxLength/2) dx -= boxLength;
			if (dx < -boxLength/2) dx += boxLength;
			if (dy >  boxLength/2) dy -= boxLength;
			if (dy < -boxLength/2) dy += boxLength;
			
			dist = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
		}
		
		return dist;
	}
	
	/**
	 * Writes the positions to a trajectory file readable by VMD.
	 */
	public static void writeTrajectory() {
		try {
			//Initialize writer and write header information
			if (trjWriter == null) {
				trjWriter =  new BufferedWriter(new FileWriter(outDir + "/trj_" + temp + ".vtf"));
				trjWriter.write("pbc " + boxLength + " " + boxLength + " " + boxLength +"\n" + 
						   "atom 0:" + (numParticles-1) + " radius " + radius + " name O type 0\n");
			}
			
			//Write system to file
			trjWriter.write("\ntimestep ordered\n");
			for (double[] p : particles) { 
				if (numDim == 3)
					trjWriter.write(p[0] + " " + p[1] + " " + p[2] + "\n");
				if (numDim == 2)
					trjWriter.write(p[0] + " " + p[1] + " 0.0\n");
			}
			
			trjWriter.flush();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates a histogram and exports a plot for the radial distribution function.
	 */
	public static void plotRDF() {
		
		int numBins = 100;
		double max = 0.5 * boxLength; 		//5.0 * sigma;		//TODO work in progress
		double[] hist = new double[numBins];
		double binSize = max / numBins;
		
		for (int i = 0; i < (numParticles - 1); i++)  {
        	for (int j = i + 1; j < numParticles; j++) {
        		double r = calcDist(i, j);
        		if (r < max) {
        			int bin = (int) (r / binSize);
        			if ((bin >= 0) && (bin < numBins))
        				hist[bin] += 2;
        		}
        	}
		}
		
		if (numDim == 3) {
			for (int i = 0; i < numBins; i++) {
				double vb = (Math.pow(i+1, 3) - Math.pow(i, 3)) * Math.pow(binSize, 3);
				double nid = (4.0/3.0) * Math.PI * vb * density;
				hist[i] = hist[i] / (numParticles * nid);
			}
		}
		
		if (numDim == 2) {
			for (int i = 0; i < numBins; i++) {
				double vb = (Math.pow(i+1, 2) - Math.pow(i, 2)) * Math.pow(binSize, 2);
				double nid = Math.PI * vb * density;
				hist[i] = hist[i] / (numParticles * nid);
			}
		}
		
		//Plotting
		@SuppressWarnings("unchecked")
		DataTable data = new DataTable(Double.class, Double.class);
		for (int i = 0; i < numBins; i++)
			data.add(i * binSize, hist[i]);
		
		XYPlot plot = new XYPlot(data);
		plot.setInsets(new Insets2D.Double(70.0, 70.0, 70.0, 70.0));
		plot.setBackground(Color.WHITE);
		 plot.getPointRenderers(data).get(0).setShape(null);
		LineRenderer lines = new DefaultLineRenderer2D();
        plot.setLineRenderers(data, lines);
        
        plot.getTitle().setText("Radial Distribution");
		plot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("r");
       
		DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
		try {
			writer.write(plot, new FileOutputStream(new File(outDir + "/rdf.png")), 800, 800);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Exports histogram from radial distribution to file.
	 */
	public static void writeRDFToFile() {
		
		int numBins = 100;
		double max = 5.0 * sigma;
		double[] hist = new double[numBins];
		double binSize = max / numBins;
		
		for (int i = 0; i < (numParticles - 1); i++)  {
        	for (int j = i + 1; j < numParticles; j++) {
        		double r = calcDist(i, j);
        		if (r < max) {
        			int bin = (int) (r / binSize);
        			if ((bin >= 0) && (bin < numBins))
        				hist[bin] += 2;
        		}
        	}
		}
		
		if (numDim == 3) {
			for (int i = 0; i < numBins; i++) {
				double vb = (Math.pow(i+1, 3) - Math.pow(i, 3)) * Math.pow(binSize, 3);
				double nid = (4.0/3.0) * Math.PI * vb * density;
				hist[i] = hist[i] / (numParticles * nid);
			}
		}
		
		if (numDim == 2) {
			for (int i = 0; i < numBins; i++) {
				double vb = (Math.pow(i+1, 2) - Math.pow(i, 2)) * Math.pow(binSize, 2);
				double nid = Math.PI * vb * density;
				hist[i] = hist[i] / (numParticles * nid);
			}
		}
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "/rdf_" + temp + ".txt"));
			for (int i = 0; i < numBins; i++)
				writer.write(i * binSize + " " + hist[i] + "\n");
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Plots a 2D system with its periodic images.
	 */
	public static void plot2DSystemPeriodic() {
		
		if (numDim == 2) {
		
			@SuppressWarnings("unchecked")
			DataTable dataTable = new DataTable(Double.class, Double.class);
			
			for (double[] p : particles) {
				dataTable.add(p[0], p[1]);
				dataTable.add(p[0] - boxLength, p[1]);
				dataTable.add(p[0] + boxLength, p[1]);
				dataTable.add(p[0], p[1] - boxLength);
				dataTable.add(p[0], p[1] + boxLength);
				dataTable.add(p[0] - boxLength, p[1] - boxLength);
				dataTable.add(p[0] + boxLength, p[1] + boxLength);
				dataTable.add(p[0] - boxLength, p[1] + boxLength);
				dataTable.add(p[0] + boxLength, p[1] - boxLength);
			}

			XYPlot plot = new XYPlot(dataTable);
			double margin = 0.25 * boxLength;			
			plot.getAxis(XYPlot.AXIS_X).setRange(-margin - boxLength, boxLength * 2 + margin);
			plot.getAxis(XYPlot.AXIS_Y).setRange(-margin - boxLength, boxLength * 2 + margin);
			
			DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
			try {
				writer.write(plot, new FileOutputStream(new File(outDir + "/system_2D.png")), 800, 800);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	public static void log(String s) {
		System.out.println(s);
	}
	
}
