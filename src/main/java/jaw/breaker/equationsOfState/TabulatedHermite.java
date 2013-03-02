package jaw.breaker.equationsOfState;

import jaw.breaker.shenzerotemperature.drivers.interpolators.Polynomials;

/**
 *
 * @author atonita
 */
public class TabulatedHermite {
    /**
     * An identifier to differentiate this table from others.
     */
    private String identifier;

    /**
     * Returns the identifier of this equation of state.
     * @return the identifier
     */
    public String getIdentifier() {
        return identifier;
    }

    /**
     * Sets the identifier for this equation of state.
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
     * The coefficient of the energyDensity polynomial in pressure.
     * The array energyDensityA[i] is the coefficients of {p^0,p^1,p^2,p^3}
     * that lead to continuity of value and derivative through piecewise
     * Hermite interpolation.
     */
    private double[][] energyDensityA;
    
    /**
     * Coefficients for extrapolating the energy density. The array
     * energyDensityExtrap[0] gives the 2 coefficients {g, c, 1/g, 1/(g-1)} that
     * continuously fit the value and derivative of the input energy density to
     * a polytrope having energy density: rho = c*p^(1/g) + p/(g-1)
     */
    private double[][] energyDensityExtrap;
    
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
     * Construct an equation of state using tabulated data. Assumes units where c=1.
     * @param logn The log_10 of the particle number density. Dimensions = 1/L^3
     * @param logp The log_10 of the pressure. Dimensions of pressure.
     * @param energyPerParticle The internal energy per particle (excludes rest mass). Dimensions of energy.
     * @param particleMass The mass of a particle.
     */
    public TabulatedHermite(double[] logn, double[] logp, double[] energyPerParticle, double particleMass) {
        N = logp.length;
        pressure = new double[N];
        for (int i = 0; i < N; i++) {
            pressure[i] = Math.pow(10, logp[i]);
        }
        
        // Energy density.
        numberDensity = new double[N];
        rho = new double[N];
        for (int i = 0; i < N; i++) {
            numberDensity[i] = Math.pow(10, logn[i]);
            rho[i] = numberDensity[i]*(particleMass + energyPerParticle[i]);
        }
        // drho/dpressure
        double[] drho = new double[N];
        drho[0] = (rho[1] - rho[0])/(pressure[1] - pressure[0]);
        drho[1] = Polynomials.differentiate(
                Polynomials.interpolatingCoefficients(
                new double[] {rho[0], rho[1], rho[2]}, 
                new double[] {pressure[0], pressure[1], pressure[2]}), pressure[1]);
        for (int i = 2; i < N-2; i++) {
            int nSpline = 5;
            double[] localRho = new double[nSpline];
            double[] localPressure = new double[nSpline];
            for (int j = 0; j < nSpline; j++) {
                int iSpline = i - nSpline/2 + j;
                localRho[j] = rho[iSpline];
                localPressure[j] = pressure[iSpline];
            }
            drho[i] = Polynomials.differentiate(
                    Polynomials.interpolatingCoefficients(
                    localRho, localPressure), pressure[i]);
        }
        drho[N-2] = Polynomials.differentiate(
                Polynomials.interpolatingCoefficients(
                new double[] {rho[N-3], rho[N-2], rho[N-1]}, 
                new double[] {pressure[N-3], pressure[N-2], pressure[N-1]}), pressure[N-2]);
        drho[N-1] = (rho[N-1] - rho[N-2])/(pressure[N-1] - pressure[N-2]);
        
        // Now we compute the interpolating coefficients.
        energyDensityA = new double[N-1][];
        for (int i = 0; i < N-1; i++) {
            double[] localRho = {rho[i], rho[i+1]};
            double[] localDrho = {drho[i], drho[i+1]};
            double[] localPressure = {pressure[i], pressure[i+1]};
            energyDensityA[i] = Polynomials.interpolatingCoefficients(localRho, localDrho, localPressure);
        }
        
        // Compute the extrapolating coefficients.
        energyDensityExtrap = new double[2][4];
        energyDensityExtrap[0][0] = (rho[0] + pressure[0])/(pressure[0]*drho[0]);
        energyDensityExtrap[1][0] = (rho[N-1] + pressure[N-1])/(pressure[N-1]*drho[N-1]);
        energyDensityExtrap[0][2] = 1/energyDensityExtrap[0][0];
        energyDensityExtrap[1][2] = 1/energyDensityExtrap[1][0];
        energyDensityExtrap[0][3] = 1/(energyDensityExtrap[0][0] - 1);
        energyDensityExtrap[1][3] = 1/(energyDensityExtrap[1][0] - 1);
        energyDensityExtrap[0][1] = -energyDensityExtrap[0][0]*energyDensityExtrap[0][3]*(pressure[0]*drho[0] - rho[0])*Math.pow(pressure[0], -energyDensityExtrap[0][2]);
        energyDensityExtrap[1][1] = -energyDensityExtrap[1][0]*energyDensityExtrap[1][3]*(pressure[N-1]*drho[N-1] - rho[N-1])*Math.pow(pressure[N-1], -energyDensityExtrap[1][2]);
    }

    /**
     * Returns the energy density, dimensions = energy/volume [ML^-1T^-2]
     * @param p The pressure p.
     * @return The energy density.
     */
    public double energyDensity(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return p*energyDensityExtrap[0][3] + energyDensityExtrap[0][1]*Math.pow(p, energyDensityExtrap[0][2]);
        } else if (memoIndex == N-1) {
            return p*energyDensityExtrap[1][3] + energyDensityExtrap[1][1]*Math.pow(p, energyDensityExtrap[1][2]);
        }
        return Polynomials.interpolate(energyDensityA[memoIndex], p);
    }
    
    /**
     * Returns the derivative of energy density, dimensions = acceleration/volume [ML^-2T^-2]
     * @param p The pressure p.
     * @return The energy density.
     */
    public double denergyDensity(double p) {
        memoize(p);
        if (memoIndex == -1) {
            return energyDensityExtrap[0][3] + energyDensityExtrap[0][1]*energyDensityExtrap[0][2]*Math.pow(p, energyDensityExtrap[0][2]-1);
        } else if (memoIndex == N-1) {
            return energyDensityExtrap[1][3] + energyDensityExtrap[1][1]*energyDensityExtrap[1][2]*Math.pow(p, energyDensityExtrap[1][2]-1);
        }
        return Polynomials.differentiate(energyDensityA[memoIndex], p);
    }
    
    /**
     * Clones the data from this table into the argument. The three fields copied are
     * in order:
     * <ol>
     * <li>particle number density</li>
     * <li>pressure</li>
     * <li>total energy density (rest mass + specific internal)</li>
     * </ol>
     * @param table A 2-d array, at least <code>double[3][]</code>
     */
    public void cloneTable(double[][] table) {
        for (int i = 0; i < 3; i++) {
            table[i] = new double[N];
        }
        System.arraycopy(numberDensity, 0, table[0], 0, N);
        System.arraycopy(pressure, 0, table[1], 0, N);
        System.arraycopy(rho, 0, table[2], 0, N);
    }
    
    private void memoize(double p) {
        memoPressure = p;
        if (p < pressure[0]) {
            memoIndex = -1;
        } else if (p > pressure[N-1]) {
            memoIndex = N-1;
        } else {
            // Binary search for bounding indices.
            memoIndex = 0;
            int upperI = N-1;
            while (upperI - memoIndex != 1) {
                int midI = (upperI + memoIndex)/2;
                if (pressure[midI] <= p) {
                    memoIndex = midI;
                } else {
                    upperI = midI;
                }
            }
        }
    }
}

