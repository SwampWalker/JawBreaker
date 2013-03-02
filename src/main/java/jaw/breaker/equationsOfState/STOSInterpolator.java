package jaw.breaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.PhysicalQuantity;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class STOSInterpolator {
	private double[] numberDensity;
	private double[] logNumberDensity;
	private double[] protonFraction;
	private double[] logProtonFraction;
	private double[][] freeEnergy;
	private double[][] neutronChemicalPotential;
	private double[][] protonChemicalPotential;
	private double[][] pressure;
	private double[] betaEquilibriumPressure;
	private int[] betaEquilibriumProtonFractionLowerIndex;
	
	// I will work in cgs units, Shen et al. work in both of the following for some quantities.
	private UnitSystem meVUnits = null;
	private UnitSystem cgs = null;
	
	private ElectronPositronPlasmaEOS epPlasma = null;
	private BlackBodyEOS blackBody = null;
	
	// This exists in case at some point I want to output multiple index files for different temperatures.
	private String betaEquilibriumIndexFilename = "betaIndexes.ind";
	private double inputTemperature = 0.1;
	
	private static final int NDENSITY = 104;	// The number of density values in the Shen et al. equation of state table
	private static final int NFRACTION = 71;	// The number of proton fraction values in the Shen et al. equation of state table.
	private static final int NUMBERLENGTH = 13;	// The length in characters of each quantity in a line of the Shen et al. equation of state table.
	private static final int NUMBERFIELDS = 17;	// The number of quantities printed on each line of the Shen et al. equation of state table.
	
	// The following give the index of the various physical quantities on a given line of the Shen et al. equation of state table.
	public static final int DENSITYINDEX = 1;
	public static final int NUMBERDENSITYINDEX = 2;
	public static final int PROTONFRACTIONINDEX = 3;
	public static final int FREEENERGYINDEX = 5;
	public static final int PRESSUREINDEX = 15;
	public static final int NEUTRONCHEMICALPOTENTIALINDEX = 16;
	public static final int PROTONCHEMICALPOTENTIALINDEX = 17;
	
	public static void main(String[] args) {
		new STOSInterpolator();
	}
	
	/**
	 * Creates the EOS: opens a hard coded table file and reads what I think is important, storing that information.
	 * Also creates the electron-positron and photon contributions.
	 */
	public STOSInterpolator() {
		try {
			// Initialise data structures.
	        numberDensity = new double[STOSInterpolator.NDENSITY];
	        logNumberDensity = new double[STOSInterpolator.NDENSITY];
	        protonFraction = new double[STOSInterpolator.NFRACTION];
	        logProtonFraction = new double[STOSInterpolator.NFRACTION];
	        betaEquilibriumProtonFractionLowerIndex = new int[STOSInterpolator.NDENSITY];
	        betaEquilibriumPressure = new double[STOSInterpolator.NDENSITY];
	        freeEnergy = new double[STOSInterpolator.NDENSITY][STOSInterpolator.NFRACTION];
	        pressure = new double[STOSInterpolator.NDENSITY][STOSInterpolator.NFRACTION];
	        protonChemicalPotential = new double[STOSInterpolator.NDENSITY][STOSInterpolator.NFRACTION];
	        neutronChemicalPotential = new double[STOSInterpolator.NDENSITY][STOSInterpolator.NFRACTION];
	        meVUnits = new UnitSystem(SIConstants.METRE.times(1E-15),SIConstants.eV.times(1E6),SIConstants.SECOND,SIConstants.COULOMB,SIConstants.k);
	        // Read the file.
	        BufferedReader in = new BufferedReader(new FileReader("C:\\Users\\atonita\\Workspace\\Spectral Methods\\eosSTOS\\lowTeos.tab"));
	        cgs = CommonUnits.CGS;
	        String str;
	        int index = 0;
	        while ((str = in.readLine()) != null) {
	            index = processLine(str, index);
	        }
	        in.close();
	        epPlasma = new ElectronPositronPlasmaEOS();
	        epPlasma.setUnits(cgs);
	        blackBody = new BlackBodyEOS();
	        blackBody.setInputUnits(cgs);
	        getBetaEquilibriumBracket();
	    } catch (IOException e) {
	    	System.out.println("Could not read file.");
	    }
	}
	
	/**
	 * Computes the beta equilibrium series, by first bracketing the beta equilibrium.
	 */
	private void getBetaEquilibriumBracket() {
		File betaEquilibriumIndexFile = new File(betaEquilibriumIndexFilename);
		if (betaEquilibriumIndexFile.exists()) {
			FileInputStream fis;
			try {
				fis = new FileInputStream(betaEquilibriumIndexFilename);
				ObjectInputStream ois = new ObjectInputStream(fis);
				betaEquilibriumProtonFractionLowerIndex = (int[]) ois.readObject();
				ois.close();
			} catch (FileNotFoundException e) {
				System.err.println("How is that possible?");
			} catch (IOException e) {
				System.err.println("Couldn't read index array.");
			} catch (ClassNotFoundException e) {
				System.err.println("Couldn't recognize that object.");
			}
		} else {
			computeBetaEquilibrium();
			outputIndexFile();
		}
		// Actually interpolate to the beta equilibrium EOS.
		interpolateBetaEquilibrium();
	}
	
	/**
	 * This subroutine interpolates the beta equilibrium values.
	 */
	private void interpolateBetaEquilibrium() {
		for (int densityIndex = 0; densityIndex < STOSInterpolator.NDENSITY; densityIndex++) {
			int fractionIndex = this.betaEquilibriumProtonFractionLowerIndex[densityIndex];
			PhysicalQuantity T = new PhysicalQuantity(inputTemperature,Dimension.TEMPERATURE,meVUnits);
			T = cgs.convert(T);
			double temperature = T.getValue();
			double[] electronChemicalPotential = new double[2];
			double[] equilibriumChemicalPotential = new double[2];
			double[] chargeNumberDensity = new double[2];
			for (int i = 0; i < 2; i++) {
				chargeNumberDensity[i] = protonFraction[fractionIndex+i]*numberDensity[densityIndex];
				electronChemicalPotential[i] = epPlasma.chemicalPotential(chargeNumberDensity[i],temperature);
				equilibriumChemicalPotential[i] = betaEquilibriumElectronChemicalPotential(densityIndex,fractionIndex+i);
			}
			double [] slope = new double[2];
			slope[0] = (equilibriumChemicalPotential[1]-equilibriumChemicalPotential[0])/(chargeNumberDensity[1]-chargeNumberDensity[0]);
			slope[1] = (electronChemicalPotential[1] - electronChemicalPotential[0])/(chargeNumberDensity[1]-chargeNumberDensity[0]);
			double numberDensityRoot = chargeNumberDensity[1]-(equilibriumChemicalPotential[1]-electronChemicalPotential[1])/(slope[0]-slope[1]);
			double exactMu = epPlasma.chemicalPotential(numberDensityRoot,temperature);
			double pressureSlope = (pressure[densityIndex][fractionIndex+1]-pressure[densityIndex][fractionIndex])
									/(chargeNumberDensity[1]-chargeNumberDensity[0]);
			betaEquilibriumPressure[densityIndex] = pressureSlope*(numberDensityRoot-chargeNumberDensity[1]) + pressure[densityIndex][fractionIndex+1];
			betaEquilibriumPressure[densityIndex] += epPlasma.pressure(exactMu, temperature);
			betaEquilibriumPressure[densityIndex] += blackBody.pressure(temperature);
			//double amu = UnitSystem.convert(SIConstants.amu.getValue(), Dimension.MASS, CommonUnits.SI, cgs);
			//System.out.println(""+numberDensity[densityIndex]*amu+"\t"+betaEquilibriumPressure[densityIndex]);
		}
	}
	
	/**
	 * If there was no index file, then we find the beta equilibrium curve through computation.
	 */
	private void computeBetaEquilibrium() {
		for (int densityIndex = 0; densityIndex < STOSInterpolator.NDENSITY; densityIndex++) {
			for (int fractionIndex = 0; fractionIndex < STOSInterpolator.NFRACTION-1; fractionIndex++) {
				double n_e = protonFraction[fractionIndex]*numberDensity[densityIndex];
				PhysicalQuantity T = new PhysicalQuantity(0.1,Dimension.TEMPERATURE,meVUnits);
				T = cgs.convert(T);
				double temperature = T.getValue();
				double mu_eOfn = epPlasma.chemicalPotential(n_e,temperature);
				double mu_e = betaEquilibriumElectronChemicalPotential(densityIndex,fractionIndex);
				double n_ep1 = protonFraction[fractionIndex+1]*numberDensity[densityIndex];
				double mu_eOfnp1 = epPlasma.chemicalPotential(n_ep1,temperature);
				double mu_ep1 = betaEquilibriumElectronChemicalPotential(densityIndex,fractionIndex+1);
				if (mu_e > mu_eOfn && mu_ep1 < mu_eOfnp1) {
					betaEquilibriumProtonFractionLowerIndex[densityIndex] = fractionIndex;
				}
			}
		}
	}
	
	/**
	 * Writes the index file which gives the lower index of the beta equilibrium proton fraction.
	 */
	private void outputIndexFile() {
		try {
			FileOutputStream fos = new FileOutputStream(betaEquilibriumIndexFilename);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(betaEquilibriumProtonFractionLowerIndex);
			oos.close();
		} catch (IOException e) {
			System.err.println("IOException");
		}
	}
	
	/**
	 * Reads a single line from the Shen et al. equation of state table and retrieves specific variables.
	 * @param str the input line
	 * @param lineIndex the index of the line
	 * @return the index of the next line unless there was a problem
	 */
	private int processLine(String str, int lineIndex) {
		if (str.length() == (STOSInterpolator.NUMBERLENGTH+1)*NUMBERFIELDS-1) {
			// Read off the required quantities.
			double baryonNumberDensity = getDouble(str,STOSInterpolator.NUMBERDENSITYINDEX);
			double logProtonFraction = getDouble(str,STOSInterpolator.PROTONFRACTIONINDEX);
			double freeEnergy = getDouble(str,STOSInterpolator.FREEENERGYINDEX);
			double neutronChemicalPotential = getDouble(str,STOSInterpolator.NEUTRONCHEMICALPOTENTIALINDEX);
			double protonChemicalPotential = getDouble(str,STOSInterpolator.PROTONCHEMICALPOTENTIALINDEX);
			double pressure = getDouble(str,STOSInterpolator.PRESSUREINDEX);
			// We needed the line index so that we could compute the numberDensity index and proton fraction index.
			int densityIndex = lineIndex%(STOSInterpolator.NDENSITY);
			int fractionIndex = lineIndex/(STOSInterpolator.NDENSITY);
			// This will constantly overwrite the old values, but I don't care.
			// Convert everything to centimetre-gram-second units.
			this.numberDensity[densityIndex] = 
				UnitSystem.convert(baryonNumberDensity, Dimension.NUMBERDENSITY, meVUnits, cgs);
			this.logNumberDensity[densityIndex] = 
				Math.log10(UnitSystem.convert(baryonNumberDensity, Dimension.NUMBERDENSITY, meVUnits, cgs));
			this.protonFraction[fractionIndex] = Math.pow(10.0,logProtonFraction);
			this.logProtonFraction[fractionIndex] = logProtonFraction;
			this.freeEnergy[densityIndex][fractionIndex] =
				UnitSystem.convert(freeEnergy, Dimension.ENERGY, meVUnits, cgs);
			this.neutronChemicalPotential[densityIndex][fractionIndex] = 
				UnitSystem.convert(neutronChemicalPotential, Dimension.ENERGY, meVUnits, cgs);
			this.protonChemicalPotential[densityIndex][fractionIndex] = 
				UnitSystem.convert(protonChemicalPotential, Dimension.ENERGY, meVUnits, cgs);
			this.pressure[densityIndex][fractionIndex] = 
				UnitSystem.convert(pressure, Dimension.PRESSURE, meVUnits, cgs);
			return lineIndex+1;
		} else {
			return lineIndex;
		}
	}
	
	/**
	 * Computes the require chemical potential of electrons for beta equilibrium.
	 * @param densityIndex the index of the density
	 * @param fractionIndex the index of the proton fraction
	 * @return the chemical potential of electrons required for beta equilibrium
	 */
	private double betaEquilibriumElectronChemicalPotential(int densityIndex, int fractionIndex) {
		return (neutronChemicalPotential[densityIndex][fractionIndex] - protonChemicalPotential[densityIndex][fractionIndex]);
	}
	
	/**
	 * Reads a line from the Shen et al. equation of state table and pulls out the quantity at a given index.
	 * @param str the line from the equation of state table
	 * @param i the index of the quantity to read
	 * @return the double value of that quantity
	 */
	private double getDouble(String str, int i) {
		int beginIndex = (i-1)*(STOSInterpolator.NUMBERLENGTH+1);
		int endIndex = beginIndex + STOSInterpolator.NUMBERLENGTH;
		String subString = str.substring(beginIndex, endIndex);
		return Double.parseDouble(subString);
	}
}
