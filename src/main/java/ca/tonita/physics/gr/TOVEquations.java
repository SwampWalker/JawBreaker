package ca.tonita.physics.gr;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.math.numerical.NewtonRaphsonPostProcessor;
import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;
import ca.tonita.math.numerical.QuasiLinearFirstOrderODESystem;
import ca.tonita.math.numerical.spectral.SpectralVector1D;
import ca.tonita.math.polynomials.LinearlyMappedBasis;

/**
 *
 * @author atonita
 */
public class TOVEquations implements QuasiLinearFirstOrderODESystem, NonLinearFirstOrderODESystem, NewtonRaphsonPostProcessor {

    private TabulatedHermite eos;
    private NonLinearFirstOrderODEBean bean;
    private double centralPressure;
    private int iPressure = 0;
    private int iMass = 1;
    private int iLambda = 2;
    private int nEquations = 3;

    public TOVEquations(TabulatedHermite eos) {
        this.eos = eos;
        bean = new NonLinearFirstOrderODEBean();
        bean.setResidue(new double[3]);
        bean.setJacobian(new double[3][6]);
    }

    public double getCentralPressure() {
        return centralPressure;
    }

    public void setCentralPressure(double centralPressure) {
        this.centralPressure = centralPressure;
    }

    public TabulatedHermite getEos() {
        return eos;
    }

    public void setEos(TabulatedHermite eos) {
        this.eos = eos;
    }

    public double[] rightHandSide(double r, double[] y) {
        return new double[]{dpdr(r, y), dmdr(r, y), dlambdadr(r, y), dndr(r, y)};
    }

    /**
     * Derivative of conserved number of particles.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of particle number wrt r
     */
    public double dndr(double r, double[] y) {
        double n = eos.numberDensity(y[iPressure]);
        if (r == 0) {
            return 0;
        }
        return 4 * Math.PI * n * r * r / Math.sqrt(1 - 2 * y[iMass] / r);
    }

    /**
     * Returns the derivative of mass with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of mass wrt r
     */
    public double dmdr(double r, double[] y) {
        double rho = eos.energyDensity(y[iPressure]);
        return 4 * Math.PI * rho * r * r;
    }

    /**
     * Returns the derivative of lambda with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of lambda wrt r
     */
    public double dlambdadr(double r, double[] y) {
        if (r == 0.0) {
            return 0;
        }
        double mass = y[1];
        double pressure = y[0];
        return (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                / (r * (r - 2 * mass));
    }

    /**
     * Returns the derivative of pressure with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of pressure wrt r
     */
    public double dpdr(double r, double[] y) {
        if (r == 0.0) {
            return 0;
        }
        double mass = y[1];
        double pressure = y[0];
        double dlambdadr = (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                / (r * (r - 2 * mass));
        double rho = eos.energyDensity(pressure);
        return -(pressure + rho) * dlambdadr;
    }

    /**
     * Returns the residue and Jacobian for the TOV equations.
     *
     * @param r the areal radial coordinate
     * @param y an array containing {mass, pressure, lambda}
     * @param dy an array of the derivatives of the variables
     * @param parameters parameters of the problem
     * @return the residue and Jacobian
     */
    @Override
    public NonLinearFirstOrderODEBean equations(int iR, double r, double[] y, double[] dy, double[] parameters, int type) {
        double[] residue = new double[3];
        double[][] jacobian = new double[3][7];
        if (type == NonLinearFirstOrderODESystem.LEFTBOUNDARY) {
            for (int i = 0; i < residue.length; i++) {
                residue[i] = dy[i];
                for (int j = 0; j < jacobian[0].length; j++) {
                    jacobian[i][j] = 0;
                }
            }
            residue[iPressure] = y[iPressure] - centralPressure;
            residue[iMass] = y[iMass];
            residue[iLambda] = y[iLambda];
            jacobian[iPressure][iPressure] = 1;
            jacobian[iMass][iMass] = 1;
            jacobian[iLambda][iLambda] = 1;
        } else if (type == NonLinearFirstOrderODESystem.RIGHTBOUNDARY) {
            for (int i = 0; i < residue.length; i++) {
                residue[i] = dy[i];
                for (int j = 0; j < jacobian[0].length; j++) {
                    jacobian[i][j] = 0;
                }
                jacobian[i][nEquations + i] = 1;
            }
            double mass = y[1];
            double dlambdadr = mass / (r * (r - 2 * mass));
            double d2lambdadrdm = (1. + 2 * r * dlambdadr) / (r * (r - 2 * mass));
            residue[iLambda] = dy[iLambda] - dlambdadr;
            jacobian[iLambda][iMass] = -d2lambdadrdm;
            double d2lambdadr2 = -(2 * r - 2 * mass) * mass / Math.pow(r * (r - 2 * mass), 2);
            jacobian[iPressure][2 * nEquations] = -dy[iPressure] / parameters[0];
            jacobian[iMass][2 * nEquations] = -dy[iMass] / parameters[0];
            jacobian[iLambda][2 * nEquations] = (-dy[iLambda] - d2lambdadr2 * r) / parameters[0];

        } else {
            // Residue.
            double mass = y[1];
            double pressure = y[0];
            double dlambdadr = (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                    / (r * (r - 2 * mass));
            double rho = eos.energyDensity(pressure);
            residue[iPressure] = dy[iPressure] + (pressure + rho) * dlambdadr;
            residue[iMass] = dy[iMass] - 4 * Math.PI * rho * r * r;
            residue[iLambda] = dy[iLambda] - dlambdadr;

            // Jacobian.
            double drho = eos.denergyDensity(pressure);
            double d2lambdadrdp = 4 * Math.PI * r * r
                    / (r - 2 * mass);
            double d2lambdadrdm = (1. + 2 * r * dlambdadr) / (r * (r - 2 * mass));
            double d2lambdadr2 = 12 * Math.PI * r * pressure / (r - 2 * mass)
                    - (2 * r - 2 * mass) * (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                    / Math.pow(r * (r - 2 * mass), 2);
            // Pressure equation.
            jacobian[iPressure][iPressure] = (1 + drho) * dlambdadr + (pressure + rho) * d2lambdadrdp;
            jacobian[iPressure][iMass] = (pressure + rho) * d2lambdadrdm;
            jacobian[iPressure][iLambda] = 0;
            jacobian[iPressure][nEquations + iPressure] = 1;
            jacobian[iPressure][nEquations + iMass] = 0;
            jacobian[iPressure][nEquations + iLambda] = 0;
            jacobian[iPressure][2 * nEquations] = (-dy[iPressure] + (pressure + rho) * d2lambdadr2 * r) / parameters[0];
            // Mass equation.
            jacobian[iMass][iPressure] = - 4 * Math.PI * drho * r * r;
            jacobian[iMass][iMass] = 0;
            jacobian[iMass][iLambda] = 0;
            jacobian[iMass][nEquations + iPressure] = 0;
            jacobian[iMass][nEquations + iMass] = 1;
            jacobian[iMass][nEquations + iLambda] = 0;
            jacobian[iMass][2 * nEquations] = (-dy[iMass] - 8 * Math.PI * rho * r * r) / parameters[0];
            // Lambda equation.
            jacobian[iLambda][iPressure] = -d2lambdadrdp;
            jacobian[iLambda][iMass] = -d2lambdadrdm;
            jacobian[iLambda][iLambda] = 0;
            jacobian[iLambda][nEquations + iPressure] = 0;
            jacobian[iLambda][nEquations + iMass] = 0;
            jacobian[iLambda][nEquations + iLambda] = 1;
            jacobian[iLambda][2 * nEquations] = (-dy[iLambda] - d2lambdadr2 * r) / parameters[0];
        }
        bean.setResidue(residue);
        bean.setJacobian(jacobian);
        return bean;
    }

    /**
     * Sets the pressure at the surface to be zero.
     *
     * @param iConstraint The index of the constraint.
     * @param indexer The indexer.
     * @param vector The data vector.
     * @param dconstraint The vector that dotted with a perturbation gives the
     * derivative of the constraint in the direction of that perturbation. I.E.
     * dconstraint.
     * @return the residual of the constraint.
     */
    public double zeroPressureConstraint(int iConstraint, ODEIndexer1D indexer, SpectralVector1D vector, double[] dconstraint) {
        for (int i = 0; i < dconstraint.length; i++) {
            dconstraint[i] = 0;
        }
        int rank = vector.getRank(0);
        dconstraint[indexer.index(0, iPressure, rank - 1)] = 1;
        return vector.getVariable(0, iPressure)[rank - 1];
    }

    @Override
    public double computeConstraint(int iConstraint, ODEIndexer1D indexer, SpectralVector1D vector, double[] dconstraint) {
        if (iConstraint == 0) {
            return zeroPressureConstraint(iConstraint, indexer, vector, dconstraint);
        }
        throw new IllegalArgumentException("Unknown constraint index: " + iConstraint);
    }

    public int getNDomains() {
        return 1;
    }

    public int[] getNVariables() {
        return new int[]{3};
    }

    @Override
    public int getNVariables(int iDomain) {
        return 3;
    }

    @Override
    public int getNParameters() {
        return 1;
    }

    /**
     * Fixes pressure after each NR step -> pressure cannot be negative.
     *
     * @param variables The variables from the SpectralVector.
     * @param parameters The parameters.
     */
    @Override
    public void postStepProcessing(SpectralVector1D vector, LinearlyMappedBasis[] bases) {
        double[][][] variables = vector.getVariables();
        double[] parameters = vector.getParameters();
        for (int iAbscissa = 0; iAbscissa < variables[0][iPressure].length; iAbscissa++) {
            if (variables[0][iPressure][iAbscissa] < 0) {
                variables[0][iPressure][iAbscissa] = -variables[0][iPressure][iAbscissa];
            }
        }
        bases[0].setRight(parameters[0]);
        vector.differentiate();
    }
}
