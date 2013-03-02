package jaw.breaker.shenzerotemperature.drivers;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import jaw.breaker.shenzerotemperature.drivers.interpolators.Interpolator;
import jaw.breaker.shenzerotemperature.drivers.interpolators.QuadraticInterpolator;

/**
 * Arises from equilibration of the chemical equation:<br>
 * <br>
 * n + \nu <-> p + e-<br>
 * <br>
 * which leads to the physical equation of equilibrium<br>
 * <br>
 * m_n + mu_n + m_nu + mu_nu = m_p + mu_p + m_e + mu_e<br>
 * <br>
 * Beta equilibrium here is the condition that mu_nu = 0.
 * 
 * @author atonita
 *
 */
public class BetaEquilibrium implements InputOutput {
	
	/**
	 * The array for the data. The first index is the
	 * charge fraction index, the second is the density
	 * index and the third is the field index. The indexing
	 * of the fields is according to Shen's hyperon table
	 * with the final field being the added electron chemical
	 * potential.
	 */
	private double data[][][];
	
	/**
	 * The array to contain the interpolated values of the fields
	 * for beta equilibrium.
	 */
	private double betaEquilibrium[][];
	
	/**
	 * The mass difference: m_p + m_e - m_n.
	 */
	private double massGap;
	
	/**
	 * The interpolator to use.
	 */
	private Interpolator interpolator;
	
	public BetaEquilibrium() {
		data = new double[110][19][66];
		betaEquilibrium = new double[110][19];
		massGap = SIConstants.mproton.getValue() 
				+ SIConstants.melectron.getValue() 
				- SIConstants.mneutron.getValue();
		massGap = UnitSystem.convert(massGap, Dimension.MASS, CommonUnits.SI, CommonUnits.MEV);
		interpolator = new QuadraticInterpolator();
	}

	/**
	 * Entry point of the routine.
	 * @param input The file to read input from, must have had leptons added.
	 * @param output The file to store output in.
	 * @return Any return value other than SUCCESS is clearly not one.
	 */
	@Override
	public int readAndWrite(File input, File output) {
		try {
			readInput(input);
			findEquilibrium();
			output(output);
		} catch (IOException e) {
			// This should never happen.
			e.printStackTrace();
		} catch (IllegalArgumentException e1) {
			return InputOutput.INCORRECTINPUT;
		}
		return InputOutput.SUCCESS;
	}

	private void output(File output) throws IOException {
		FileWriter fstream = new FileWriter(output);
		BufferedWriter out = new BufferedWriter(fstream);
		
		out.write("# " + interpolator.getDescriptor() + "\n");
		for (int iRho = 0; iRho < 110; iRho++) {
			for (int iField = 0; iField < nFields-1; iField++) {
				out.write(betaEquilibrium[iRho][iField] + "  ");
			}
			out.write(betaEquilibrium[iRho][nFields-1] + "\n");
		}
		
		out.flush();
		fstream.flush();
		out.close();
		fstream.close();
	}

	/**
	 * Finds the beta equilibrium sequence.
	 */
	private void findEquilibrium() {
		for (int iRho = 0; iRho < 110; iRho++) {
			double[] muNu = new double[66];
			for (int iYp = 0; iYp < 66; iYp++) {
				double muN = data[iRho][neutronPotential][iYp];
				double muP = data[iRho][protonPotential][iYp];
				double muE = data[iRho][electronPotential][iYp];
				muNu[iYp] = muP + muE - muN + massGap;
			}
			double zero = interpolator.zeroPoint(muNu);
			// Make sure we found it.
			if (zero < 0) {
				System.err.println("Could not find beta equilibrium for density " + data[iRho][density][0]);
				for (int i = 0; i < 66; i++) {
					System.out.println(i + " " + muNu[i]);
				}
				System.exit(-1);
			}
			
			// interpolate the fields
			for (int iField = 0; iField < nFields; iField++) {
				betaEquilibrium[iRho][iField] = interpolator.interpolate(zero, data[iRho][iField]);
			}
		}
	}

	/**
	 * Main routine. I made it so that I could have the try-catch
	 * in the readAndWrite() function.
	 * @param input The file to read input from, must have had leptons added.
	 * @param output The file to store output in.
	 * @throws IOException 
	 */
	private void readInput(File input) throws IOException {
		FileInputStream rawInput = new FileInputStream(input);
		DataInputStream inData = new DataInputStream(rawInput);
		BufferedReader in = new BufferedReader(new InputStreamReader(inData));

		// First 3 lines are comments;
		for (int i = 0; i < 4; i++) in.readLine();
		
		// 66 Yp blocks
		for (int j = 0; j < 66; j++) {
			// with 110 density lines
			for (int i = 0; i < 110; i++) {
				processLine(in.readLine(), j, i);
			}
			// Read the block separating line.
			in.readLine();
		}
	}

	/**
	 * Reads all the data from a line into the data array.
	 * @param readLine The line read from the file.
	 * @param iYp The proton number index (the block number)
	 * @param iRho The density index (the line number inside a block)
	 */
	private void processLine(String readLine, int iYp, int iRho) {
		String[] tokens = readLine.split("\\s+");
		for (int i = 0; i < 19; i++) {
			data[iRho][i][iYp] = Double.valueOf(tokens[i]);
		}
	}

}
