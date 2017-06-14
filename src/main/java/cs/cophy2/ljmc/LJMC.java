package cs.cophy2.ljmc;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;


public class LJMC {

	public static int dim = 2;
	public static int numParticles;			//Number of particles
	public static double boxLength;			//Box length, has to be even
	public static List<double[]> particles;	//Positions of the particles
	public static double energy;			//Energy of system
	public static double temp;				//Temperature
	public static int numSteps;				//Number of MC steps to do
	public static double epsilon;			//Epsilon
	public static double sigma;				//Sigma
	public static double trunc;				//Truncation for interactions
	public static double radius;			//Radius of the spheres
	public static double volSphere;			//Volume of each sphere
	public static double volBox;			//Volume of the box
	public static double volFraction;		//Volume fraction
	public static int numBins;				//Number of bins for RDF calculation
	
	public static long seed;
	public static Random rand;
	public static NumberFormat nf;
	public static String outDir;
	public static BufferedWriter trajectoryWriter;
	
	/**
	 * Monte Carlo simulation of Lennard-Jones Particles.
	 */
	public static void main(String[] args) {
		
		setUp();
		
		//Parameters
		boxLength = 100;
		radius = 1.0;		
		volSphere = Math.PI	* Math.pow(radius, 2);;		
		volFraction = 0.2;
		temp = 2.0;
		numSteps = 1000000;
		epsilon = 5.0;			//1.0 to 10.0 (no to strong interaction)
		sigma = 1.0;
		trunc = 3 * sigma;
		numBins = 200;
		volBox  = Math.pow(boxLength, 2);	
		numParticles = new Double((volFraction * volBox) / volSphere).intValue();
		
		//Initialize system
		System.out.println("Initializing the system ...");
		particles = new ArrayList<double[]>(numParticles);
		for (int i = 0; i < numParticles; i++) {
			double[] pos = new double[dim];
			for (int j = 0; j < pos.length; j++)
				pos[j] = rand.nextDouble() * boxLength;
			particles.add(pos);
		}
		
		writeToTrajectoryFile();
	
		//Calc initial energy
		energy = calcEnergy();
		System.out.println("Initial energy = " + energy);
		
		//MC loop
		for (int i = 0; i < numSteps; i++) {
			
			//Choose a particle at random and calcuate its energy
			int particle = rand.nextInt(numParticles);
			double prevEnergy = calcEnergy(particle);
			
			//Move the choosen particle and calculate the new energy of the particle
			double[] prevParticle = particles.get(particle).clone();
			moveParticle(particle);
			double newEnergy = calcEnergy(particle);
			
			//Accept the new configuration with a probability corresponding to Boltzmann statistics
			double deltaEnergy = newEnergy - prevEnergy;
			boolean doMove = false;
			if (deltaEnergy < 0) 
				doMove = true;
			else {
				double prob = rand.nextDouble();
				if (prob < Math.exp(-(1.0/temp) * deltaEnergy))
					doMove = true;
				else
					particles.set(particle, prevParticle);  //Don't move, reset configuration
			}
			
			//Update energy
			if (doMove) 
				energy += deltaEnergy;
		}
		
		System.out.println("Final energy = " + energy);
		
		writeToTrajectoryFile();
		
		calcRDF();
		
		tearDown();
		
		System.out.println("Finished.");
	}
	
	/**
	 * Initialization routine.
	 */
	public static void setUp() {
		seed = 123;
		rand = new Random(seed);
		outDir = "/home/alepfu/Desktop/LJMC";
		nf = DecimalFormat.getInstance();
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);
		try {
			trajectoryWriter = new BufferedWriter(new FileWriter(outDir + "/ljmc.vtf"));
			trajectoryWriter.write("pbc " + boxLength + " " + boxLength + " " + boxLength +"\n" + 
								   "atom 0:" + (numParticles-1) +  radius + " name O type 0\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Cleanup routine.
	 */
	public static void tearDown() {
		try {
			trajectoryWriter.flush();
			trajectoryWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Appends the acutal configuration to the trajectory file.
	 */
	public static void writeToTrajectoryFile() {
		try {
			trajectoryWriter.write("\ntimestep ordered\n");
			for (double[] x : particles) 
				trajectoryWriter.write(x[0] + " " + x[1] + " 1.0\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Moves a particle at random with applying PBC
	 */
	public static void moveParticle(int p) {
		particles.get(p)[0] += rand.nextDouble() * 2 - 1;
		particles.get(p)[1] += rand.nextDouble() * 2 - 1;
		if (particles.get(p)[0] > boxLength) particles.get(p)[0] -= boxLength;
		if (particles.get(p)[0] < 0) particles.get(p)[0] += boxLength;
		if (particles.get(p)[1] > boxLength) particles.get(p)[1] -= boxLength;
		if (particles.get(p)[1] < 0) particles.get(p)[1] += boxLength;
	}

	/**
	 * Calculates the energy of the system.
	 */
	public static double calcEnergy() {
		double energy = 0.0;
		for (int i = 0; i < (numParticles-1); i++) {
			for (int j = i + 1; j < numParticles; j++) {
				double distSq = calcSquaredDistance(i, j);
				if (distSq <= Math.pow(trunc, 2))
					energy += 4 * epsilon * (1 / Math.pow(distSq, 6)) - (1 / Math.pow(distSq, 3));
			}
		}
		return energy;
	}
	
	/**
	 * Calculates the energy of a single particle.
	 */
	public static double calcEnergy(int p) {
		double energy = 0.0;
		int i = 0;
		for (int j = 0; j < numParticles; j++) {
			if (i != p) {
				double distSq = calcSquaredDistance(p, j);
				if (distSq <= Math.pow(trunc, 2))
					energy += 4 * epsilon * (1 / Math.pow(distSq, 6)) - (1 / Math.pow(distSq, 3));
			}
			i++;
		}
		return energy;
	}
	
	/**
	 * Calculates the squared distances between two particles,
	 * restricting to the minimum image ("wrapping distances").
	 */
	public static double calcSquaredDistance(int i, int j) {
		double[] p1 = particles.get(i);
		double[] p2 = particles.get(j);
		double dx = p1[0] - p2[0];
		double dy = p1[1] - p2[1];
		if (dx >  boxLength/2) dx -= boxLength;
		if (dx < -boxLength/2) dx += boxLength;
		if (dy >  boxLength/2) dy -= boxLength;
		if (dy < -boxLength/2) dy += boxLength;
		return Math.pow(dx, 2) + Math.pow(dy, 2) /* + Math.pow(dz, 2)*/;
	}
	
	/**
	 * Calculates the distances between two particles,
	 * restricting to the minimum image ("wrapping distances").
	 */
	public static double calcDistance(int i, int j) {
		return Math.sqrt(calcSquaredDistance(i, j));
	}
	
	/**
	 * Calculates a histogram and exports a plot for the radial distribution function.
	 */
	public static void calcRDF() {
		
		System.out.println("Calculating RDF ...");
		
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
			writer.write(plot, new FileOutputStream(new File(outDir + "/rdf.png")), 800, 800);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
