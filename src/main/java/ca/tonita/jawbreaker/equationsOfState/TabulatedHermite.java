package ca.tonita.jawbreaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;
import ca.tonita.jawbreaker.shenzerotemperature.drivers.interpolators.Polynomials;

/**
 *
 * @author atonita
 */
public class TabulatedHermite {

    /**
     * An identifier to differentiate this table from others.
     */
    private String identifier;
    private final double particleMass;
    private double edgePressure;
    private double edgeDensity;

    /**
     * Returns the identifier of this equation of state.
     *
     * @return the identifier
     */
    public String getIdentifier() {
        return identifier;
    }

    /**
     * Sets the identifier for this equation of state.
     *
     * @param identifier the identifier to set
     */
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    @Override
    public String toString() {
        return identifier;
    }
    /**
     * The pressure array.
     */
    private double[] pressure;
    /**
     * The number density array.
     */
    private double[] numberDensity;
    /**
     * Total energy density.
     */
    private double[] rho;
    /**
     * The coefficient of the energyDensity polynomial in pressure. The array
     * energyDensityA[i] is the coefficients of {p^0,p^1,p^2,p^3} that lead to
     * continuity of value and derivative through piecewise Hermite
     * interpolation.
     */
    private double[][] energyDensityA;
    /**
     * The coefficient of the number density polynomial in pressure. The array
     * energyDensityA[i] is the coefficients of {p^0,p^1,p^2,p^3} that lead to
     * continuity of value and derivative through piecewise Hermite
     * interpolation.
     */
    private double[][] numberDensityA;
    /**
     * The coefficients of the energy per particle in pressure.
     */
    private final double[][] energyA;
    /**
     * Coefficients for extrapolating the energy per particle. A polytrope has
     * energy per particle given by e = k^(1/g) * m / (g - 1) * p ^ (1 - 1/g).
     * This array holds the coefficients {k^(1/g)*m/(g-1), 1 - 1/g}.
     */
    private double[][] energyExtrap;
    /**
     * Coefficients for extrapolating the energy density. The array
     * energyDensityExtrap[0] gives the 2 coefficients {g, c, 1/g, 1/(g-1)} that
     * continuously fit the value and derivative of the input energy density to
     * a polytrope having energy density: rho = c*p^(1/g) + p/(g-1)
     */
    private double[][] energyDensityExtrap;
    /**
     * Coefficients for extrapolating the number density. The array
     * numberDensityExtrap[0] gives the 2 coefficients {g, c, 1/g} that
     * continuously fit the value and derivative of the input energy density to
     * a polytrope having number density: n = c*p^(1/g)
     */
    private double[][] numberDensityExtrap;
    /**
     * The number of points.
     */
    private int N;
    /**
     * The index to the element of the pressure array which is a lower bound for
     * the memo pressure.
     */
    private int memoIndex;
    /**
     * The last pressure seen by the tabulated equation of state. To avoid
     * recomputing the bounding pressure.
     */
    private double memoPressure;

    /**
     * Construct an equation of state using tabulated data. Assumes units where
     * c=1.
     *
     * @param logn The log_10 of the particle number density. Dimensions = 1/L^3
     * @param logp The log_10 of the pressure. Dimensions of pressure.
     * @param energyPerParticle The internal energy per particle (excludes rest
     * mass). Dimensions of energy.
     * @param particleMass The mass of a particle.
     */
    public TabulatedHermite(double[] logn, double[] logp, double[] energyPerParticle, double particleMass) {
        N = logp.length;
        pressure = new double[N];
        this.particleMass = particleMass;
        for (int i = 0; i < N; i++) {
            pressure[i] = Math.pow(10, logp[i]);
        }

        // Energy density.
        numberDensity = new double[N];
        rho = new double[N];
        for (int i = 0; i < N; i++) {
            numberDensity[i] = Math.pow(10, logn[i]);
            rho[i] = numberDensity[i] * (particleMass + energyPerParticle[i]);
        }

        // The derivatives.
        double[] drho = dvariable(rho);
        double[] dn = dvariable(numberDensity);
        double[] denergyPerParticle = dvariable(energyPerParticle);
        // Now we compute the interpolating coefficients.
        energyDensityA = computeHermiteCoefficients(rho, drho);
        numberDensityA = computeHermiteCoefficients(numberDensity, dn);
        energyA = computeHermiteCoefficients(energyPerParticle, denergyPerParticle);

        // Compute the extrapolating coefficients.
        energyDensityExtrap = new double[2][4];
        energyDensityExtrap[0][0] = (rho[0] + pressure[0]) / (pressure[0] * drho[0]);
        energyDensityExtrap[1][0] = (rho[N - 1] + pressure[N - 1]) / (pressure[N - 1] * drho[N - 1]);
        energyDensityExtrap[0][2] = 1 / energyDensityExtrap[0][0];
        energyDensityExtrap[1][2] = 1 / energyDensityExtrap[1][0];
        energyDensityExtrap[0][3] = 1 / (energyDensityExtrap[0][0] - 1);
        energyDensityExtrap[1][3] = 1 / (energyDensityExtrap[1][0] - 1);
        energyDensityExtrap[0][1] = -energyDensityExtrap[0][0] * energyDensityExtrap[0][3] * (pressure[0] * drho[0] - rho[0]) * Math.pow(pressure[0], -energyDensityExtrap[0][2]);
        energyDensityExtrap[1][1] = -energyDensityExtrap[1][0] * energyDensityExtrap[1][3] * (pressure[N - 1] * drho[N - 1] - rho[N - 1]) * Math.pow(pressure[N - 1], -energyDensityExtrap[1][2]);
        numberDensityExtrap = new double[2][3];
        numberDensityExtrap[0][0] = energyDensityExtrap[0][0];
        numberDensityExtrap[0][2] = energyDensityExtrap[0][2];
        numberDensityExtrap[0][1] = numberDensity[0] / Math.pow(pressure[0], numberDensityExtrap[0][2]);
        numberDensityExtrap[1][0] = energyDensityExtrap[1][0];
        numberDensityExtrap[1][2] = energyDensityExtrap[1][2];
        numberDensityExtrap[1][1] = numberDensity[N - 1] / Math.pow(pressure[N - 1], numberDensityExtrap[1][2]);
        // Goal as always is to extrapolate a polytrope through the last two points...
        energyExtrap = new double[2][2];
        energyExtrap[0][1] = Math.log(energyPerParticle[0] / energyPerParticle[1]) / Math.log(pressure[0] / pressure[1]);
        energyExtrap[0][0] = energyPerParticle[0] * Math.pow(pressure[0], -energyExtrap[0][1]);
        energyExtrap[1][1] = Math.log(energyPerParticle[N - 1] / energyPerParticle[N - 2]) / Math.log(pressure[N - 1] / pressure[N - 2]);
        energyExtrap[1][0] = energyPerParticle[N - 1] * Math.pow(pressure[N - 1], -energyExtrap[1][1]);
        setEdgePressure();
    }

    /**
     * Returns the energy density, dimensions = energy/volume [ML^-1T^-2]
     *
     * @param p The pressure p.
     * @return The energy density.
     */
    public double energyDensity(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return p * energyDensityExtrap[0][3] + energyDensityExtrap[0][1] * Math.pow(p, energyDensityExtrap[0][2]);
        } else if (memoIndex == N - 1) {
            return p * energyDensityExtrap[1][3] + energyDensityExtrap[1][1] * Math.pow(p, energyDensityExtrap[1][2]);
        }
        return Polynomials.interpolate(energyDensityA[memoIndex], p);
    }

    /**
     * Returns the derivative of energy density with respect to pressure,
     * dimensions = acceleration/volume [ML^-2T^-2]
     *
     * @param p The pressure p.
     * @return The energy density.
     */
    public double denergyDensity(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return energyDensityExtrap[0][3] + energyDensityExtrap[0][1] * energyDensityExtrap[0][2] * Math.pow(p, energyDensityExtrap[0][2] - 1);
        } else if (memoIndex == N - 1) {
            return energyDensityExtrap[1][3] + energyDensityExtrap[1][1] * energyDensityExtrap[1][2] * Math.pow(p, energyDensityExtrap[1][2] - 1);
        }
        return Polynomials.differentiate(energyDensityA[memoIndex], p);
    }

    /**
     * Clones the data from this table into the argument. The fields copied are
     * in order: <ol> <li>particle number density</li> <li>pressure</li>
     * <li>total energy density (rest mass + specific internal)</li>
     * <li>derivative of total energy density w.r.t pressure</li></ol>
     *
     * @param table A 2-d array, at least <code>double[4][]</code>
     */
    public void cloneTable(double[][] table) {
        for (int i = 0; i < 4; i++) {
            table[i] = new double[N];
        }
        System.arraycopy(numberDensity, 0, table[0], 0, N);
        System.arraycopy(pressure, 0, table[1], 0, N);
        System.arraycopy(rho, 0, table[2], 0, N);
        for (int i = 0; i < N; i++) {
            table[3][i] = denergyDensity(pressure[i]);
        }
    }

    private void memoize(double p) {
        memoPressure = p;
        if (p < pressure[0]) {
            memoIndex = -1;
        } else if (p > pressure[N - 1]) {
            memoIndex = N - 1;
        } else {
            // Binary search for bounding indices.
            memoIndex = 0;
            int upperI = N - 1;
            while (upperI - memoIndex != 1) {
                int midI = (upperI + memoIndex) / 2;
                if (pressure[midI] <= p) {
                    memoIndex = midI;
                } else {
                    upperI = midI;
                }
            }
        }
    }

    /**
     * Returns the number density as a function of pressure.
     *
     * @param p the pressure
     * @return the number density
     */
    public double numberDensity(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return numberDensityExtrap[0][1] * Math.pow(p, numberDensityExtrap[0][2]);
        } else if (memoIndex == N - 1) {
            return numberDensityExtrap[1][1] * Math.pow(p, numberDensityExtrap[1][2]);
        }
        return Polynomials.interpolate(numberDensityA[memoIndex], p);
    }

    /**
     * Returns the derivative of number density as a function of pressure.
     *
     * @param p the pressure
     * @return the number density
     */
    public double dnumberDensity(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return numberDensityExtrap[0][1] * numberDensityExtrap[0][2] * Math.pow(p, numberDensityExtrap[0][2] - 1);
        } else if (memoIndex == N - 1) {
            return numberDensityExtrap[1][1] * numberDensityExtrap[1][2] * Math.pow(p, numberDensityExtrap[1][2] - 1);
        }
        return Polynomials.differentiate(numberDensityA[memoIndex], p);
    }

    /**
     * Computes the derivative of the input variable using splines.
     *
     * @param variable the variable to differentiate
     * @return the derivative of the variable at the same points.
     */
    private double[] dvariable(double[] variable) {
        // drho/dpressure
        double[] dvariable = new double[N];
        dvariable[0] = (variable[1] - variable[0]) / (pressure[1] - pressure[0]);
        dvariable[1] = Polynomials.differentiate(
                Polynomials.interpolatingCoefficients(
                new double[]{variable[0], variable[1], variable[2]},
                new double[]{pressure[0], pressure[1], pressure[2]}), pressure[1]);
        for (int i = 2; i < N - 2; i++) {
            int nSpline = 5;
            double[] localRho = new double[nSpline];
            double[] localPressure = new double[nSpline];
            for (int j = 0; j < nSpline; j++) {
                int iSpline = i - nSpline / 2 + j;
                localRho[j] = variable[iSpline];
                localPressure[j] = pressure[iSpline];
            }
            dvariable[i] = Polynomials.differentiate(
                    Polynomials.interpolatingCoefficients(
                    localRho, localPressure), pressure[i]);
        }
        dvariable[N - 2] = Polynomials.differentiate(
                Polynomials.interpolatingCoefficients(
                new double[]{variable[N - 3], variable[N - 2], variable[N - 1]},
                new double[]{pressure[N - 3], pressure[N - 2], pressure[N - 1]}), pressure[N - 2]);
        dvariable[N - 1] = (variable[N - 1] - variable[N - 2]) / (pressure[N - 1] - pressure[N - 2]);
        return dvariable;
    }

    /**
     * Computes the hermite interpolating coefficients of the variable.
     *
     * @param variable the variable to get interpolating coefficients of
     * @param dvariable the derivative of the variable
     * @return the interpolating coefficients in each interval
     */
    private double[][] computeHermiteCoefficients(double[] variable, double[] dvariable) {
        double[][] coefficients = new double[N - 1][];
        for (int i = 0; i < N - 1; i++) {
            double[] localVariable = {variable[i], variable[i + 1]};
            double[] localDVariable = {dvariable[i], dvariable[i + 1]};
            double[] localPressure = {pressure[i], pressure[i + 1]};
            coefficients[i] = Polynomials.interpolatingCoefficients(localVariable, localDVariable, localPressure);
        }
        return coefficients;
    }

    public double getParticleMass() {
        return particleMass;
    }

    /**
     * Computes the value of the edge density.
     */
    private void setEdgePressure() {
        double nEdge = UnitSystem.convert(0.08, Dimension.NUMBERDENSITY, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO) * UnitSystem.convert(SIConstants.mneutron.getValue(), Dimension.MASS, CommonUnits.MKS, CommonUnits.GEOMETRICASTRO) / particleMass;
        int iGuess = 0;
        while (numberDensity[iGuess] < nEdge) {
            iGuess++;
        }
        double guess = pressure[iGuess];
        double nGuess = numberDensity[iGuess];
        while (Math.abs(nGuess - nEdge) > 1.0E-10) {
            guess += (nEdge - nGuess) / dnumberDensity(guess);
            nGuess = numberDensity(guess);
        }
        edgePressure = guess;
        edgeDensity = nEdge;
    }

    /**
     * Returns the energy per particle at the pressure p.
     *
     * @param p the pressure p
     * @return the energy per particle
     */
    public double energyPerParticle(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return energyExtrap[0][0] * Math.pow(p, energyExtrap[0][1]);
        } else if (memoIndex == N - 1) {
            return energyExtrap[1][0] * Math.pow(p, energyExtrap[1][1]);
        }
        return Polynomials.interpolate(energyA[memoIndex], p);
    }

    /**
     * Returns the derivative of the energy per particle at the pressure p.
     *
     * @param p the pressure p
     * @return the derivative of energy per particle
     */
    public double dEnergyPerParticle(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return energyExtrap[0][0] * energyExtrap[0][1] * Math.pow(p, energyExtrap[0][1] - 1);
        } else if (memoIndex == N - 1) {
            return energyExtrap[1][0] * energyExtrap[1][1] * Math.pow(p, energyExtrap[1][1] - 1);
        }
        return Polynomials.differentiate(energyA[memoIndex], p);
    }

    public double getEdgePressure() {
        return edgePressure;
    }

    public double getEdgeDensity() {
        return edgeDensity;
    }
}
