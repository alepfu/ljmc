package cs.cophy2.ljmc;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Monte Carlo simulation of Lennard-Jones particles.
 * 
 * @author Alexander Pfundner
 * 
 */
public class LJMCSimulation {

	public int numParticles;			
	public double boxLength;			
	public List<double[]> particles;			
	public double temp;				
	public int numEqSteps;
	public int numSampSteps;
	public int numTotalSteps;
	public double epsilon;			
	public double sigma;				
    public double cutoff;				
    public double radius;			
	public double density;
	public double beta;
	public double boxVolume;
	public double displ;
	
	public int printProgressFreq;
	public int writeTrjFreq;
	
	public double avgEnergy;
	public double pressure;
	public double finalAcceptRate;
	
	public NumberFormat nf;
	public static BufferedWriter trjWriter;
	public Random rand;
	
	public String workDir = "/home/alepfu/Desktop/LJMC";
	
	/**
	 * Main method, use for testing purposes only.
	 */
	public static void main(String[] args) {
		
		double d = 0.2;
		double t = 1.0;
		double co = 2.5;
		double eps = 5.0;
		int n = 100;
		
		LJMCSimulation sim = new LJMCSimulation(d, t, co, eps, n);
		sim.run();
		
	}
	
	/**
	 * An instance of a Monte Carlo Lennard-Jones simulation.
	 * Call run() method for running the simulation, evaluate results
	 * via attributes avgEnergy, pressure and finalAcceptRate.
	 * 
	 * @param d The density.
	 * @param t The temperature.
	 * @param co The cutoff distance.
	 */
	public LJMCSimulation(double d, double t, double co, double eps, int n) {
		
		density = d;
		temp = t;
		cutoff = co;
		numParticles = n;						
		epsilon = eps;
		
		sigma = 1.0;
		numEqSteps = 100000;
		numSampSteps = 100000;
		displ = 0.5;
		
		boxVolume = numParticles / density;
		boxLength = Math.pow(boxVolume, 1.0 / 3);
		beta = 1.0 / temp;
		radius = 0.5 * sigma;
		numTotalSteps = numEqSteps + numSampSteps;
		
		printProgressFreq = 1000;
		writeTrjFreq = 100;
		
		trjWriter = null;
		rand = new Random();
		nf = DecimalFormat.getInstance();
		nf.setMaximumFractionDigits(4);
		nf.setMinimumFractionDigits(4);
		
	}
	
	/**
	 * Run the simulation with the given parameters from the constructor.
	 */
	public void run() {
		
		//Log parameters
		log("Num. particles = " + numParticles);
		log("Epsilon = " + nf.format(epsilon));
		log("Sigma = " + nf.format(sigma));
		log("Cutoff = " + nf.format(cutoff));
		log("Temp. = " + nf.format(temp));
		log("Density = " + nf.format(density));
		log("Box length = " + nf.format(boxLength) + "\n");
		
		//Initialize system
		initSystemFrequentlyDistributed();
		
		//Setup trajectory file and add initial configuration
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
			double[] prevParticle = new double[3];
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
			
			//Reset sums/counter for sampling
			if (step == numEqSteps) {
				energySum = 0.0;
				virialSum = 0.0;
				acceptCounter = 0;
			}
			
			//Print progress
			if (step % printProgressFreq == 0)
				log((step < numEqSteps ? "Equilibration" : "Sampling")
						+ " [ " + (int)((step * 1.0 / numTotalSteps) * 100) + "% ]");
		}
		
		//Calc and log results
	    avgEnergy = energySum / numSampSteps / numParticles;
	    double pressureTailCorr = (16.0 / 3.0) * Math.PI * Math.pow(density, 2) * epsilon * Math.pow(sigma, 3) * (  (2.0/3.0) * Math.pow(sigma / cutoff, 9) - Math.pow(sigma / cutoff, 3)  );
	    pressure = virialSum / 3.0 / numSampSteps / boxVolume + density * temp + pressureTailCorr;
	    finalAcceptRate = acceptCounter * 1.0 / numSampSteps * 100.0;
		log("Avg. energy = " + nf.format(avgEnergy));
		log("Avg. pressure = " + nf.format(pressure));
		log("Acceptance rate = " + nf.format(finalAcceptRate));
	    
		//Write the last configuration and then close trajectory file
		try {
			writeTrajectory();
			trjWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	/**
	 * Initializes the system with randomly choosen positions.
	 */
	public void initSystemRandomlyDistributed() {
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
	public void initSystemFrequentlyDistributed() {
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
	public void moveParticle(int p) {
		
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
		particles.get(p)[2] += (rand.nextDouble() - 0.5) * displ;
		if (particles.get(p)[2] >= boxLength)
			particles.get(p)[2] -= boxLength;
		if (particles.get(p)[2] < 0.0)
			particles.get(p)[2] += boxLength;
		
	}

	/**
	 * Calculates the energy of the system.
	 */
	public double calcEnergyTotal() {
		double energyTotal = 0.0;
		
		for (int i = 0; i < numParticles; i++)
			energyTotal += calcEnergyParticle(i);
		
		energyTotal = 0.5 * energyTotal;
		
		//Tail correction
		double energyTailCorr = (8.0/3.0) * Math.PI * density * epsilon * Math.pow(sigma, 3) * (  (1.0/3.0) * Math.pow(sigma/cutoff, 9) - Math.pow(sigma/cutoff, 3)  );
		energyTotal += numParticles * energyTailCorr;
		
		return energyTotal;
	}

	/**
	 * Calculates the energy of a single particle.
	 */
	public double calcEnergyParticle(int p) {
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
	public double calcLJPotential(double dist) {
		
		double pot = 4.0 * epsilon * (  Math.pow(sigma / dist, 12) - Math.pow(sigma / dist, 6)  ); 

		//Shift correction
		double energyShiftCorr = 4.0 * epsilon * (  Math.pow(sigma / cutoff, 12) - Math.pow(sigma / cutoff, 6)  ) ;
		pot -= energyShiftCorr;
		
		return pot;
	}
	
	/**
	 * Calculates the virial of the system.
	 */
	public double calcVirialTotal() {
		double virialTotal = 0.0;
		
		for (int i = 0; i < numParticles; i++)
			virialTotal += calcVirialParticle(i);

		return 0.5 * virialTotal;
	}
	
	/**
	 * Calculates the virial of a single particle.
	 */
	public double calcVirialParticle(int p) {
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
	public double calcVirial(double dist) {
		
		return 48.0 * epsilon * (  Math.pow(sigma/dist, 12) - 0.5 * Math.pow(sigma/dist, 6)  );
	}
	
	/**
	 * Calculates the distance between two particles,
	 * restricting to the minimum image ("wrapping distances").
	 */
	public double calcDist(int i, int j) {

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
		
		return Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2));
	}
	
	/**
	 * Writes the positions to a trajectory file.
	 */
	public void writeTrajectory() {
		try {
			//Initialize writer and write header information
			if (trjWriter == null) {
				trjWriter =  new BufferedWriter(new FileWriter(workDir + "/trj.vtf"));
				trjWriter.write("pbc " + boxLength + " " + boxLength + " " + boxLength +"\n" + 
						   "atom 0:" + (numParticles-1) + " radius " + radius + " name O type 0\n");
			}
			
			//Write system to file
			trjWriter.write("\ntimestep ordered\n");
			for (double[] p : particles)  
					trjWriter.write(p[0] + " " + p[1] + " " + p[2] + "\n");
			
			trjWriter.flush();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void log(String s) {
		System.out.println(s);
	}
	
}
