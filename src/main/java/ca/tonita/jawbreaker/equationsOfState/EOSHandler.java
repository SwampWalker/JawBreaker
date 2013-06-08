/*
 * Copyright (C) 2012 atonita
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package ca.tonita.jawbreaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * A Singleton which creates TabulatedHermite equations of state.
 *
 * @author atonita
 */
public class EOSHandler {

    /**
     * The singleton.
     */
    private EOSHandler eosHandler;

    /**
     * Force singleton.
     */
    private EOSHandler() {
    }

    /**
     * Creates a TabulatedHermite equation of state. The equation of state is
     * logarithmic in the particle number density spacing.
     *
     * @param k The polytropic constant k or kappa to use.
     * @param gamma The adiabatic index to use.
     * @param particleMass The average mass of a particle to use.
     * @param nPoints The number of points to use in the table.
     * @param logn0 The initial (log10) number density of the table.
     * @param lognStep The step size to use between different table points.
     * @return A polytropic equation of state to use.
     */
    public static TabulatedHermite polytrope(double k, double gamma, double particleMass, int nPoints, double logn0, double lognStep) {
        double[] logn = new double[nPoints];
        double[] logp = new double[nPoints];
        double[] spE = new double[nPoints];
        double[] A = new double[nPoints];
        double[] Z = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            logn[i] = logn0 + i * lognStep;
            logp[i] = Math.log10(polyP(Math.pow(10, logn[i]), k, gamma));
            spE[i] = polySpE(Math.pow(10, logn[i]), k, gamma);
            A[i] = 56;
            Z[i] = 56;
        }
        return new TabulatedHermite(logn, logp, spE, particleMass, A, Z);
    }

    /**
     * Returns the polytropic pressure. Assumes a particle mass of 1.
     *
     * @param n The particle number.
     * @param k The polytropic constant k or kappa.
     * @param gamma The adiabatic index gamma.
     * @return The pressure.
     */
    private static double polyP(double n, double k, double gamma) {
        return k * Math.pow(n, gamma);
    }

    /**
     * Returns the polytropic specific energy. Assumes a particle mass of 1.
     *
     * @param n The particle number.
     * @param k The polytropic constant k or kappa.
     * @param gamma The adiabatic index gamma.
     * @return The specific energy.
     */
    private static double polySpE(double n, double k, double gamma) {
        return k * Math.pow(n, gamma - 1) / (gamma - 1);
    }
    
    
    /**
     * Reads an equation of state from file, returns a <code>TabulateHermite</code> object.
     * @param eosFile The file to open.
     * @return the <code>TabulatedHermite</code> from the file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static TabulatedHermite readTableFromFile(File eosFile) throws FileNotFoundException, IOException {
        // Read the file.
        FileInputStream fis = new FileInputStream(eosFile);
        BufferedReader br = new BufferedReader(new InputStreamReader(fis));
        String line = br.readLine();
        ArrayList<String[]> tokenSet = new ArrayList<String[]>();
        int nLineOne = -1;
        while (line != null) {
            if (!line.startsWith("#")) {
                String[] tokens = line.split("\\s+");
                if (nLineOne == -1) {
                    nLineOne = tokens.length;
                }
                if (tokens.length == nLineOne) {
                    tokenSet.add(tokens);
                }
            }
            line = br.readLine();
        }
        fis.close();

        double[] logn = new double[tokenSet.size()];
        double[] logp = new double[tokenSet.size()];
        double[] energyPerParticle = new double[tokenSet.size()];
        double[] A = new double[tokenSet.size()];
        double[] Z = new double[tokenSet.size()];
        
        // Unit conversion
        double nConversion = Math.log10(UnitSystem.convert(1, Dimension.NUMBERDENSITY, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO));
        double pConversion = Math.log10(UnitSystem.convert(1, Dimension.PRESSURE, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO));
        double eConversion = UnitSystem.convert(1, Dimension.ENERGY, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO);
        double particleMass = UnitSystem.convert(931.494, Dimension.MASS, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO); // value from shen guide
        
        for (int i = 0; i < tokenSet.size(); i++) {
            String[] tokens = tokenSet.get(i);
            logn[i] = Double.valueOf(tokens[1]) + nConversion;
            logp[i] = Double.valueOf(tokens[13]) + pConversion;
            A[i] = Double.valueOf(tokens[6]); // Tolerance added to avoid NaN in the A=0 core.
            Z[i] = Double.valueOf(tokens[7]);
            if (A[i] == 0.0) {
                A[i] = A[i-1];
            }
            if (Z[i] == 0.0) {
                Z[i] = Z[i-1];
            }
            energyPerParticle[i] = Double.valueOf(tokens[4])*eConversion;
            if (i != 0) {
                if (logp[i] < logp[i-1]) {
                    System.out.println("Problem reading table at entry " + (i-1) + "-" + i + ", non-monotonic pressure.");
                    System.out.println(logp[i-1] - pConversion + " " + tokens[13]);
                }
            }
        }
        
        TabulatedHermite eos = new TabulatedHermite(logn, logp, energyPerParticle, particleMass, A, Z);
        return eos;
    }
}
