package cs.cophy2.ljmc;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
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


public class LJMC {

	public static int numDim;
	public static int numParticles;			
	public static double boxLength;			
	public static List<double[]> particles;			
	public static double temp;				
	public static int numSteps;				
	public static double epsilon;			
	public static double sigma;				
    public static double trunc;				
    public static double truncSq;
    public static double radius;			
	public static double density;
	public static double beta;
	
	public static List<Double> acceptRates;
	public static List<Double> energies;
	
	public static int printFreq;
	public static int trjFreq;
	
	public static long seed = 123;
	public static Random rand = new Random(seed);
	public static NumberFormat nf = DecimalFormat.getInstance();
	public static String outDir = "/home/alepfu/Desktop/LJMC";
	public static String outSuffix = "";
	public static BufferedWriter trjWriter;
	
	/**
	 * Monte Carlo simulation of Lennard-Jones Particles.
	 */
	public static void main(String[] args) {
		
		nf.setMaximumFractionDigits(4);
		nf.setMinimumFractionDigits(4);
		acceptRates = new ArrayList<Double>();
		energies = new ArrayList<Double>();
		
		try {
			
			//Parameters
			numDim = 2;
			numParticles = (int) Math.pow(10, numDim);
			epsilon = 5.0;							
			sigma = 1.0;
			radius = 0.5 * sigma;			
			numSteps = 1000000;
			printFreq = (int) (numSteps * 0.01);
			trjFreq = (int) (numSteps * 0.00001);
			temp = 1.0;
			beta = 1.0 / temp;
			density = 0.1;
			boxLength = Math.pow(numParticles / density, 1.0 / numDim);
			trunc = 2.5 * sigma;
			truncSq = Math.pow(trunc, 2);
			
			//Log parameters
			System.out.println("Num. dimensions = " + numDim);
			System.out.println("Num. particles = " + numParticles);
			System.out.println("Epsilon = " + nf.format(epsilon));
			System.out.println("Temperature = " + nf.format(temp));
			System.out.println("Density = " + nf.format(density));
			System.out.println("Box length = " + nf.format(boxLength));
			System.out.println("Dist. truncation = " + nf.format(trunc));
			
			//Initialize system
			initSystemFrequentlyDistributed();
			writeTrajectory();
			
			//Calc initial energy
			double energy = calcEnergySystem();
			System.out.println("Initial energy = " + energy);
			
			//MC loop
			int numAcceptSteps = 0;
			for (int i = 0; i < numSteps; i++) {
				
				//Choose a particle at random and calculate its energy
				int particleIndex = rand.nextInt(numParticles);
				double prevEnergy = calcEnergyParticle(particleIndex);
				
				//Move the choosen particle and calculate the new energy of the particle
				double[] prevParticle = new double[numDim];
				for (int j = 0; j < prevParticle.length; j++)
					prevParticle[j] = particles.get(particleIndex)[j];
				moveParticle(particleIndex);
				double newEnergy = calcEnergyParticle(particleIndex);
				
				//Acceptance
				double deltaEnergy = newEnergy - prevEnergy;
				if ((deltaEnergy < 0) || (rand.nextDouble() < Math.exp(-beta * deltaEnergy))) {
					energy += deltaEnergy;
					++numAcceptSteps;
					
					if (i % trjFreq == 0)
						writeTrajectory();
					
				} else {
					particles.set(particleIndex, prevParticle);
				}
				
				//Collect the total energies
				energies.add(energy);
				
				//Calc and collect acceptance rates
				acceptRates.add(numAcceptSteps * 1.0 / (i+1));
				
				//Print progress
				if (i % printFreq == 0)
					System.out.print("\r" + (int)((i * 1.0 / numSteps) * 100) + "% ");
			}
			
			//Results
			System.out.println("\nFinal total energy = " + nf.format(energy));
			System.out.println("Avg. total energy = " + nf.format(calcAvgEnergy()));
			plotEnergies();
			System.out.println("Final acceptance rate = " + nf.format(numAcceptSteps * 1.0 / numSteps));
			System.out.println("Avg. acceptance rate = " + nf.format(calcAvgAcceptRate()));
			plotAcceptanceRate();
			plotRadialDistribution();

			
			
			trjWriter.close();
			System.out.println("\nFinished.");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
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
		
		if (numDim == 3) {
			particles.get(p)[0] += rand.nextDouble() * 2 - 1;
			particles.get(p)[1] += rand.nextDouble() * 2 - 1;
			particles.get(p)[2] += rand.nextDouble() * 2 - 1;
			
			if (particles.get(p)[0] > boxLength) particles.get(p)[0] -= boxLength;
			if (particles.get(p)[0] < 0) particles.get(p)[0] += boxLength;
			if (particles.get(p)[1] > boxLength) particles.get(p)[1] -= boxLength;
			if (particles.get(p)[1] < 0) particles.get(p)[1] += boxLength;
			if (particles.get(p)[2] > boxLength) particles.get(p)[2] -= boxLength;
			if (particles.get(p)[2] < 0) particles.get(p)[2] += boxLength;
		}
		
		if (numDim == 2) {
			particles.get(p)[0] += rand.nextDouble() * 2 - 1;
			particles.get(p)[1] += rand.nextDouble() * 2 - 1;
			
			if (particles.get(p)[0] > boxLength) particles.get(p)[0] -= boxLength;
			if (particles.get(p)[0] < 0) particles.get(p)[0] += boxLength;
			if (particles.get(p)[1] > boxLength) particles.get(p)[1] -= boxLength;
			if (particles.get(p)[1] < 0) particles.get(p)[1] += boxLength;
		}
	}

	/**
	 * Calculates the energy of the system.
	 */
	public static double calcEnergySystem() {
		double energySystem = 0.0;
		
		for (int i = 0; i < numParticles; i++)
			energySystem += calcEnergyParticle(i);
		
		return 0.5 * energySystem;
	}

	/**
	 * Calculates the energy of a single particle.
	 */
	public static double calcEnergyParticle(int p) {
		double energyParticle = 0.0;
		
		for (int j = 0; j < numParticles; j++) {
			if (j != p) {
				double dist = calcDist(p, j);
				if (dist <= trunc)
					energyParticle += calcLJPotential(dist);
			}
		}
		
		return energyParticle;
	}
	
	/**
	 * Calculates the Lennard-Jones potential for a given distance.
	 * U(r) = 4 * epsilon * [ (sigma/r)^12 - (sigma/r)^6 ]
	 */
	public static double calcLJPotential(double dist) {
		
		return 4.0 * epsilon * (  Math.pow(sigma / dist, 12) - Math.pow(sigma / dist, 6)  );
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
			
			dist = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2)  + Math.pow(dz, 2));
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
				trjWriter =  new BufferedWriter(new FileWriter(outDir + "/ljmc.vtf"));
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
	 * Plots the evolution of the total energy.
	 */
	public static void plotEnergies() {
		@SuppressWarnings("unchecked")
		DataTable data = new DataTable(Integer.class, Double.class);
		for (int i = 0; i < energies.size(); i += 100 )
			data.add(i, energies.get(i));
		
		XYPlot plot = new XYPlot(data);
		plot.setInsets(new Insets2D.Double(70.0, 70.0, 70.0, 70.0));
		plot.setBackground(Color.WHITE);
		 plot.getPointRenderers(data).get(0).setShape(null);
		LineRenderer lines = new DefaultLineRenderer2D();
        plot.setLineRenderers(data, lines);
        
        plot.getTitle().setText("Total energy");
		plot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("Iteration");
       
		DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
		try {
			writer.write(plot, new FileOutputStream(new File(outDir + "/energy.png")), 800, 800);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Plots the evolution of the acceptance rate.
	 */
	public static void plotAcceptanceRate() {
		@SuppressWarnings("unchecked")
		DataTable data = new DataTable(Integer.class, Double.class);
		for (int i = 0; i < acceptRates.size(); i += 100 )
			data.add(i, acceptRates.get(i));
		
		XYPlot plot = new XYPlot(data);
		plot.setInsets(new Insets2D.Double(70.0, 70.0, 70.0, 70.0));
		plot.setBackground(Color.WHITE);
		 plot.getPointRenderers(data).get(0).setShape(null);
		LineRenderer lines = new DefaultLineRenderer2D();
        plot.setLineRenderers(data, lines);
        
        plot.getTitle().setText("Acceptance rate");
		plot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("Iteration");
       
		DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
		try {
			writer.write(plot, new FileOutputStream(new File(outDir + "/ar.png")), 800, 800);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates a histogram and exports a plot for the radial distribution function.
	 */
	public static void plotRadialDistribution() {
		
		int numBins = 100;
		double max = 5.0 * sigma;  //0.5 * boxLength;
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
	
	public static double calcAvgEnergy() {
		double sum = 0.0;
		for (double e : energies)
			sum += e;
		return sum / numSteps;
	}
	
	public static double calcAvgAcceptRate() {
		double sum = 0.0;
		for (double ar : acceptRates)
			sum += ar;
		return sum / numSteps;
	}
}
