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
     * Returns the single instance of the eosHandler.
     *
     * @return the instance of the EOSHandler
     */
    public EOSHandler getInstance() {
        if (eosHandler == null) {
            eosHandler = new EOSHandler();
        }
        return eosHandler;
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
}
