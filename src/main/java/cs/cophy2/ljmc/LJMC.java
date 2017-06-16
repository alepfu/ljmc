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


public class LJMC {

	public static int numParticles;			
	public static double boxLength;			
	public static List<double[]> particles;	
	public static double energy;			
	public static double temp;				
	public static int numSteps;				
	public static double epsilon;			
	public static double sigma;				
    public static double trunc;				
    public static double truncSq;
    public static double radius;			
	public static double density;
	public static double beta;
	
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
		
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);
		
		try {
			
			//Parameters
			numParticles = 50;
			epsilon = 5.0;			//1.0 to 10.0 (no to strong interaction)
			sigma = 1.0;
			radius = 0.5 * sigma;			
			numSteps = 1000000;
			printFreq = 10000;
			trjFreq = 10000;
			temp = 2.0;
			beta = 1.0 / temp;
			density = 0.15;
			boxLength = Math.cbrt(numParticles / density);
			System.out.println("Box length = " + nf.format(boxLength));
			trunc = 0.5 * boxLength; //3.0 * sigma;
			truncSq = Math.pow(trunc, 2);
			
			//Initialize system
			initSystemFrequentlyDistributed();
			
			//Calc initial energy
			double energy = calcEnergySystem();
			System.out.println("Initial energy = " + energy);
			System.out.println("Initial avg. energy = " + (energy / numParticles));
			
			calcRDF("start");
			
			//MC loop
			int acceptedSteps = 0;
			for (int i = 0; i < numSteps; i++) {
				
				//Write sytem to trajectory file
				if (i % trjFreq == 0)
					writeTrajectory();
				
				//Choose a particle at random and calculate its energy
				int particle = rand.nextInt(numParticles);
				double prevEnergy = calcEnergyParticle(particle);
				
				//Move the choosen particle and calculate the new energy of the particle
				double[] prevParticle = particles.get(particle).clone();
				moveParticle(particle);
				double newEnergy = calcEnergyParticle(particle);
				
				//Accept the new configuration with a probability corresponding to Boltzmann statistics
				double deltaEnergy = newEnergy - prevEnergy;
				if (deltaEnergy < 0) {
					energy += deltaEnergy;
					++acceptedSteps;
				}
				else {
					double prob = rand.nextDouble();
					if (prob < Math.exp(-beta * deltaEnergy)) {
						energy += deltaEnergy;
						++acceptedSteps;
					} else {
						particles.set(particle, prevParticle);  //Don't move, reset configuration
					}
				}
				
				//Print progress
				if (i % printFreq == 0)
					System.out.print("\r" + nf.format((i * 1.0 / numSteps)*100) + "%");
			}
			
			System.out.println("\nFinal energy = " + energy);
			System.out.println("Final avg. energy = " + (energy / numParticles));
			System.out.println("Acceptance rate = " + nf.format(acceptedSteps * 1.0 / numSteps));
			
			
			calcRDF("end");
			
			
			System.out.println("Finished.");
			
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
			double[] pos = new double[3];		
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
		
		//Get number of partitions for the lowest perfect cube
		int numPart = 2;
        while (Math.pow(numPart, 3) < numParticles)
            numPart = numPart + 1;

        //Get particle positions by advancing an index
        int[] index = {0, 0, 0};
		for (int i = 0; i < numParticles; i++) {
			
			double[] pos = new double[3];		
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
	
	/**
	 * Moves a particle at random with applying PBC.
	 */
	public static void moveParticle(int p) {
		
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

	/**
	 * Calculates the energy of the system.
	 */
	public static double calcEnergySystem() {
		double energy = 0.0;
		for (int i = 0; i < (numParticles-1); i++) {
			for (int j = i + 1; j < numParticles; j++) {
				double distSq = calcSquaredDistance(i, j);
				if (distSq <= truncSq)
					energy += 4 * epsilon * (  Math.pow(sigma / distSq, 6)) - (Math.pow(sigma / distSq, 3)  );
			}
		}
		return energy;
	}

	/**
	 * Calculates the energy of a single particle.
	 */
	public static double calcEnergyParticle(int p) {
		double energyParticle = 0.0;
		for (int j = 0; j < numParticles; j++) {
			if (j != p) {
				double distSq = calcSquaredDistance(p, j);
				if (distSq <= truncSq)
					energyParticle += 4 * epsilon * (  Math.pow(sigma / distSq, 6)) - (Math.pow(sigma / distSq, 3)  );
			}
		}
		return energyParticle;
	}
	
	/**
	 * Calculates the squared distance between two particles,
	 * restricting to the minimum image ("wrapping distances").
	 */
	public static double calcSquaredDistance(int i, int j) {
		double[] p1 = particles.get(i);
		double[] p2 = particles.get(j);
		double dx = p1[0] - p2[0];
		double dy = p1[1] - p2[1];
		double dz = p1[2] - p2[2];
		if (dx >  boxLength/2) dx -= boxLength;
		if (dx < -boxLength/2) dx += boxLength;
		if (dy >  boxLength/2) dy -= boxLength;
		if (dy < -boxLength/2) dy += boxLength;
		if (dz >  boxLength/2) dz -= boxLength;
		if (dz < -boxLength/2) dz += boxLength;
		return Math.pow(dx, 2) + Math.pow(dy, 2)  + Math.pow(dz, 2);
	}
	
	/**
	 * Calculates the distance between two particles,
	 * restricting to the minimum image ("wrapping distances").
	 */
	public static double calcDistance(int i, int j) {
		return Math.sqrt(calcSquaredDistance(i, j));
	}
	
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
			for (double[] x : particles) 
				trjWriter.write(x[0] + " " + x[1] + " " + x[2] + "\n");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates a histogram and exports a plot for the radial distribution function.
	 */
	public static void calcRDF(String plotLabel) {
		
		int numBins = 200;
		double max = 0.5 * boxLength;
		double[] hist = new double[numBins];
		double binSize = max / numBins;
		for (int i = 0; i < (numParticles - 1); i++)  {
        	for (int j = i + 1; j < numParticles; j++) {
        		double r = calcDistance(i, j);
        		if (r < max) {
        			int bin = (int) (r / binSize);
        			if ((bin >= 0) && (bin < numBins))
        				hist[bin] += 2;
        		}
        	}
		}
		
		for (int i = 0; i < numBins; i++) {
			double r = binSize * (i + 0.5);
			double vb = (Math.pow(i+1, 3) - Math.pow(i, 3)) * Math.pow(binSize, 3);
			double nid = (4.0/3.0) * Math.PI * vb * r;   //Fehler im Buch??
			hist[i] = hist[i] / (numParticles * nid);
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
        
        plot.getTitle().setText("Radial Distribution Function");
		plot.getAxisRenderer(XYPlot.AXIS_X).getLabel().setText("r");
		plot.getAxisRenderer(XYPlot.AXIS_Y).getLabel().setText("g(r)");
       
		DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
		try {
			writer.write(plot, new FileOutputStream(new File(outDir + "/rdf_" + plotLabel +".png")), 800, 800);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
