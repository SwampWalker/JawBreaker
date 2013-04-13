package ca.tonita.jawbreaker.shenzerotemperature.drivers;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import ca.tonita.jawbreaker.shenzerotemperature.drivers.interpolators.Polynomials;

/**
 * Writes a simple cpp class to file. Also changes the units
 * to something more suitable for general relativity.
 */
public class CppClassMaker implements InputOutput {
	
	private static final int NPOINTS = 110;
	private double[] numberDensity;
	private double[] pressure;
	private double[] internalEnergy;
	private double[] lameR;
	private double[] lameQ;

	private double[] dpressure;
	private double[] dinternalEnergy;
	private double[] dlameR;
	private double[] dlameQ;
	
	private boolean polytropeOverride = false;
	private double gamma = 2;
	private double kappa = 100;
	
	private double particleMass;
	private double dn;
	
	private UnitSystem inputUnits = null;
	private UnitSystem outputUnits = null;
	
	public CppClassMaker() {
		numberDensity = new double[NPOINTS];
		pressure = new double[NPOINTS];
		internalEnergy = new double[NPOINTS];
		lameR = new double[NPOINTS];
		lameQ = new double[NPOINTS];
		dpressure = new double[NPOINTS];
		dinternalEnergy = new double[NPOINTS];
		dlameR = new double[NPOINTS];
		dlameQ = new double[NPOINTS];
		
		inputUnits = CommonUnits.MEV;
		outputUnits = CommonUnits.GEOMETRICASTRO;
		particleMass = UnitSystem.convert(931.494, Dimension.MASS, inputUnits, outputUnits); // Value is from Shen.
	}
	
	public void setPolytrope(double gamma, double kappa) {
		this.kappa = kappa;
		this.gamma = gamma;
		this.polytropeOverride = true;
	}

	@Override
	public int readAndWrite(File input, File output) {
		try {
			readData(input);
			if (polytropeOverride) {
				overwriteEOS();
			}
			computeDerived();
			outputData(output);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return InputOutput.SUCCESS;
	}

	/**
	 * Overwrite the pressure and energy terms with a polytropic pressure.
	 */
	private void overwriteEOS() {
		for (int i = 0; i < numberDensity.length; i++) {
			double rho = Math.pow(10, numberDensity[i])*particleMass;
			pressure[i] = kappa*Math.pow(rho, gamma);
			internalEnergy[i] = pressure[i]/((gamma - 1)*rho);
		}
	}

	/**
	 * Write the class to file.
	 * @param output The file to output to.
	 * @throws IOException 
	 */
	private void outputData(File output) throws IOException {
		String name = output.getName().replace(".cpp", "");
		FileWriter fstream = new FileWriter(output);
		BufferedWriter out = new BufferedWriter(fstream);
		
		out.write("double " + name + "::particleMass = " + String.valueOf(particleMass) + ";\n");
		out.write("double " + name + "::dchidlogn = " + String.valueOf(1./dn) + ";\n");
		arrayOut(name, "log10numberDensity", numberDensity, out);
		arrayOut(name, "log10pressure", pressure, out);
		arrayOut(name, "dlog10pressure", dpressure, out);
		arrayOut(name, "internalEnergy", internalEnergy, out);
		arrayOut(name, "dinternalEnergy", dinternalEnergy, out);
		arrayOut(name, "lameR", lameR, out);
		arrayOut(name, "dlameR", dlameR, out);
		arrayOut(name, "lameQ", lameQ, out);
		arrayOut(name, "dlameQ", dlameQ, out);
		
		out.flush();
		fstream.flush();
		out.close();
		fstream.close();
	}

	private void arrayOut(String name, String string, double[] array, BufferedWriter out) throws IOException {
		out.write("double " + name + "::" + string + "[] = {");
		for (int i = 0; i < array.length; i++) {
			out.write(String.valueOf(array[i]));
			if (i != array.length - 1) out.write(", ");
		}
		out.write("};\n");
	}

	/**
	 * Computes the derived quantities: derivatives, bulk modulus.
	 */
	private void computeDerived() {
		double[][] fields = {pressure, internalEnergy, lameR};
		double[][] derivatives = {dpressure, dinternalEnergy, dlameR};
		
		dn = (numberDensity[numberDensity.length-1] - numberDensity[0])/(numberDensity.length - 1.);
		
		for (int iField = 0; iField < 3; iField++) {
			double[] field = fields[iField];
			double[] derivative = derivatives[iField];
			for (int i = 0; i < field.length; i++) {
				int shift = getShift(i, field.length);
				double[] neighbours = {field[i + shift], field[i + shift + 1], field[i + shift + 2],
						field[i + shift + 3], field[i + shift + 4]};
				double[] coeffs = Polynomials.interpolatingCoefficients(neighbours);
				derivative[i] = Polynomials.differentiate(coeffs, -shift, dn);
			}
		}
		
		// We can compute the bulk modulus and lameQ now.
		for (int i = 0; i < lameQ.length; i++) {
			// K = n dp/dn, and lambda = K - 2/3mu, then we divide by m*n to get q and r.
			lameQ[i] = Math.pow(10, pressure[i] - numberDensity[i])*dpressure[i]/particleMass - 2./3.*lameR[i];
		}
		for (int i = 0; i < lameQ.length; i++) {
			int shift = getShift(i, lameQ.length);
			double[] neighbours = {lameQ[i + shift], lameQ[i + shift + 1], lameQ[i + shift + 2],
					lameQ[i + shift + 3], lameQ[i + shift + 4]};
			double[] coeffs = Polynomials.interpolatingCoefficients(neighbours);
			dlameQ[i] = Polynomials.differentiate(coeffs, -shift, dn);
		}
	}

	/**
	 * Given an 1-dimensional array that we want to take a spline of,
	 * and given an index i, what is the beginning index such that for indices
	 * i + shift, i + shift + 1, ... i + shift + 4, we both have data in the array and the points
	 * are i's nearest neighbours. This method returns the shift.
	 * @param i The index of the point that we want the nearest neighbours of.
	 * @param length The length of the array.
	 * @return The shift I spoke of above.
	 */
	private int getShift(int i, int length) {
		if (i < 0 || i >= length) throw new IllegalArgumentException("Index out of bounds.");
		if (i == 0) return 0;
		else if (i == 1) return -1;
		else if (i == length - 1) return -4;
		else if (i == length - 2) return -3;
		return -2; 
	}

	/**
	 * Reads the data from file, converts the units.
	 * @param input The beta-equilibrium equation of state with shear added file.
	 * @throws IOException 
	 */
	private void readData(File input) throws IOException {
		// A lot of work just to readLine() >:( That's boiler plate for you.
		FileInputStream rawInput = new FileInputStream(input);
		DataInputStream inData = new DataInputStream(rawInput);
		BufferedReader in = new BufferedReader(new InputStreamReader(inData));

		in.readLine(); // Header line describes order of interpolation.
		// with 110 density lines
		for (int i = 0; i < 110; i++) {
			try {
				String[] tokens = in.readLine().split("\\s+");
				numberDensity[i] = Math.log10(UnitSystem.convert(Math.pow(10, Double.valueOf(tokens[InputOutput.numberdensity])), Dimension.NUMBERDENSITY, inputUnits, outputUnits));
				pressure[i] = Math.log10(UnitSystem.convert(Math.pow(10, Double.valueOf(tokens[InputOutput.pressure])), Dimension.PRESSURE, inputUnits, outputUnits));
				internalEnergy[i] = UnitSystem.convert(Double.valueOf(tokens[InputOutput.internalEnergy]), Dimension.ENERGY, inputUnits, outputUnits);
				// What we read in is actually mu/n
				lameR[i] = UnitSystem.convert(Double.valueOf(tokens[InputOutput.shear]), Dimension.ENERGY, inputUnits, outputUnits);
				lameR[i] /= particleMass;
			} catch (ArrayIndexOutOfBoundsException e) {
				System.err.println("i = " + i);
				throw e;
			}
		}
	}

}
