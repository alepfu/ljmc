package cs.cophy2.ljmc;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;

/**
 * Calculates and plots the RDF.
 * 
 * @author Alexander Pfundner
 *
 */
public class RadialDistribution {

	public int numParticles;
	public List<double[]> particles;
	public double density;
	public double boxLength;
	public int numBins;
	double maxr;
	double[] hist;
	double binSize;
	
	public static String workDir = "/home/alepfu/Desktop/LJMC";
	
	public RadialDistribution(int np, double d, double bl) {
		
		numParticles = np;
		density = d;
		boxLength = bl;
		
		particles = new ArrayList<double[]>(np);
		numBins = 100;
		maxr = 0.5 * boxLength;
		hist = new double[numBins];
		binSize = maxr / numBins;
		
		//Get final configuration from trajectory file
        try {
			BufferedReader br = new BufferedReader(new FileReader(workDir + "/trj.vtf"));
			String l;
			List<String> lines = new ArrayList<String>();
			while((l = br.readLine()) != null)
				lines.add(l);
			for(int i = lines.size()-1; i >= 0; i--) 
				if (lines.get(i).equals("timestep ordered"))
					break;
				else {
					double[] p = new double[3];
					String[] split = lines.get(i).split("\\s+");
					p[0] = Double.parseDouble(split[0]);
					p[1] = Double.parseDouble(split[1]);
					p[2] = Double.parseDouble(split[2]);
					particles.add(p);
				}
			
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates RDF and exports data to file.
	 */
	public void run() {

		//Calc histogram
		for (int i = 0; i < (numParticles - 1); i++)  {
        	for (int j = i + 1; j < numParticles; j++) {
        		double r = calcDist(i, j);
        		if (r < maxr) {
        			int bin = (int) (r / binSize);
        			if ((bin >= 0) && (bin < numBins))
        				hist[bin] += 2;
        		}
        	}
		}
		for (int i = 0; i < numBins; i++) {
			double vb = (Math.pow(i+1, 3) - Math.pow(i, 3)) * Math.pow(binSize, 3);
			double nid = (4.0/3.0) * Math.PI * vb * density;
			hist[i] = hist[i] / (numParticles * nid);
		}
		
		//Export to file
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter(workDir + "/rdf.txt"));
			for (int i = 0; i < numBins; i++)
				w.write(i * binSize + " " + hist[i] + "\n");
			w.flush();
			w.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Plots the RDF.
	 */
	public void plot() {
		
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
			writer.write(plot, new FileOutputStream(new File(workDir + "/rdf.png")), 800, 800);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
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

}