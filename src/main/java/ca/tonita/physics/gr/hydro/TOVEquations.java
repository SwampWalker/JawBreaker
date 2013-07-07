package ca.tonita.physics.gr.hydro;

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
        return new double[]{dpdr(r, y), dmdr(r, y), dlambdadr(r, y), dm0dr(r, y)};
    }

    /**
     * Derivative of rest mass.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of rest mass wrt r
     */
    public double dm0dr(double r, double[] y) {
        double n = eos.numberDensity(y[TOVIndex.PRESSURE]);
        if (r == 0) {
            return 0;
        }
        return eos.getParticleMass()*4 * Math.PI * n * r * r / Math.sqrt(1 - 2 * y[TOVIndex.MASS] / r);
    }

    /**
     * Returns the derivative of mass with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of mass wrt r
     */
    public double dmdr(double r, double[] y) {
        double rho = eos.energyDensity(y[TOVIndex.PRESSURE]);
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
        double mass = y[TOVIndex.MASS];
        double pressure = y[TOVIndex.PRESSURE];
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
        double mass = y[TOVIndex.MASS];
        double pressure = y[TOVIndex.PRESSURE];
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
            residue[TOVIndex.PRESSURE] = y[TOVIndex.PRESSURE] - centralPressure;
            residue[TOVIndex.MASS] = y[TOVIndex.MASS];
            residue[TOVIndex.LAMBDA] = y[TOVIndex.LAMBDA];
            jacobian[TOVIndex.PRESSURE][TOVIndex.PRESSURE] = 1;
            jacobian[TOVIndex.MASS][TOVIndex.MASS] = 1;
            jacobian[TOVIndex.LAMBDA][TOVIndex.LAMBDA] = 1;
        } else if (type == NonLinearFirstOrderODESystem.RIGHTBOUNDARY) {
            for (int i = 0; i < residue.length; i++) {
                residue[i] = dy[i];
                for (int j = 0; j < jacobian[0].length; j++) {
                    jacobian[i][j] = 0;
                }
                jacobian[i][nEquations + i] = 1;
            }
            double mass = y[TOVIndex.MASS];
            double dlambdadr = mass / (r * (r - 2 * mass));
            double d2lambdadrdm = (1. + 2 * r * dlambdadr) / (r * (r - 2 * mass));
            residue[TOVIndex.LAMBDA] = dy[TOVIndex.LAMBDA] - dlambdadr;
            jacobian[TOVIndex.LAMBDA][TOVIndex.MASS] = -d2lambdadrdm;
            double d2lambdadr2 = -(2 * r - 2 * mass) * mass / Math.pow(r * (r - 2 * mass), 2);
            jacobian[TOVIndex.PRESSURE][2 * nEquations] = -dy[TOVIndex.PRESSURE] / parameters[0];
            jacobian[TOVIndex.MASS][2 * nEquations] = -dy[TOVIndex.MASS] / parameters[0];
            jacobian[TOVIndex.LAMBDA][2 * nEquations] = (-dy[TOVIndex.LAMBDA] - d2lambdadr2 * r) / parameters[0];

        } else {
            // Residue.
            double mass = y[TOVIndex.MASS];
            double pressure = y[TOVIndex.PRESSURE];
            double dlambdadr = (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                    / (r * (r - 2 * mass));
            double rho = eos.energyDensity(pressure);
            residue[TOVIndex.PRESSURE] = dy[TOVIndex.PRESSURE] + (pressure + rho) * dlambdadr;
            residue[TOVIndex.MASS] = dy[TOVIndex.MASS] - 4 * Math.PI * rho * r * r;
            residue[TOVIndex.LAMBDA] = dy[TOVIndex.LAMBDA] - dlambdadr;

            // Jacobian.
            double drho = eos.denergyDensity(pressure);
            double d2lambdadrdp = 4 * Math.PI * r * r
                    / (r - 2 * mass);
            double d2lambdadrdm = (1. + 2 * r * dlambdadr) / (r * (r - 2 * mass));
            double d2lambdadr2 = 12 * Math.PI * r * pressure / (r - 2 * mass)
                    - (2 * r - 2 * mass) * (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                    / Math.pow(r * (r - 2 * mass), 2);
            // Pressure equation.
            jacobian[TOVIndex.PRESSURE][TOVIndex.PRESSURE] = (1 + drho) * dlambdadr + (pressure + rho) * d2lambdadrdp;
            jacobian[TOVIndex.PRESSURE][TOVIndex.MASS] = (pressure + rho) * d2lambdadrdm;
            jacobian[TOVIndex.PRESSURE][TOVIndex.LAMBDA] = 0;
            jacobian[TOVIndex.PRESSURE][nEquations + TOVIndex.PRESSURE] = 1;
            jacobian[TOVIndex.PRESSURE][nEquations + TOVIndex.MASS] = 0;
            jacobian[TOVIndex.PRESSURE][nEquations + TOVIndex.LAMBDA] = 0;
            jacobian[TOVIndex.PRESSURE][2 * nEquations] = (-dy[TOVIndex.PRESSURE] + (pressure + rho) * d2lambdadr2 * r) / parameters[0];
            // Mass equation.
            jacobian[TOVIndex.MASS][TOVIndex.PRESSURE] = - 4 * Math.PI * drho * r * r;
            jacobian[TOVIndex.MASS][TOVIndex.MASS] = 0;
            jacobian[TOVIndex.MASS][TOVIndex.LAMBDA] = 0;
            jacobian[TOVIndex.MASS][nEquations + TOVIndex.PRESSURE] = 0;
            jacobian[TOVIndex.MASS][nEquations + TOVIndex.MASS] = 1;
            jacobian[TOVIndex.MASS][nEquations + TOVIndex.LAMBDA] = 0;
            jacobian[TOVIndex.MASS][2 * nEquations] = (-dy[TOVIndex.MASS] - 8 * Math.PI * rho * r * r) / parameters[0];
            // Lambda equation.
            jacobian[TOVIndex.LAMBDA][TOVIndex.PRESSURE] = -d2lambdadrdp;
            jacobian[TOVIndex.LAMBDA][TOVIndex.MASS] = -d2lambdadrdm;
            jacobian[TOVIndex.LAMBDA][TOVIndex.LAMBDA] = 0;
            jacobian[TOVIndex.LAMBDA][nEquations + TOVIndex.PRESSURE] = 0;
            jacobian[TOVIndex.LAMBDA][nEquations + TOVIndex.MASS] = 0;
            jacobian[TOVIndex.LAMBDA][nEquations + TOVIndex.LAMBDA] = 1;
            jacobian[TOVIndex.LAMBDA][2 * nEquations] = (-dy[TOVIndex.LAMBDA] - d2lambdadr2 * r) / parameters[0];
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
        dconstraint[indexer.index(0, TOVIndex.PRESSURE, rank - 1)] = 1;
        return vector.getVariable(0, TOVIndex.PRESSURE)[rank - 1];
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
        for (int iAbscissa = 0; iAbscissa < variables[0][TOVIndex.PRESSURE].length; iAbscissa++) {
            if (variables[0][TOVIndex.PRESSURE][iAbscissa] < 0) {
                variables[0][TOVIndex.PRESSURE][iAbscissa] = -variables[0][TOVIndex.PRESSURE][iAbscissa];
            }
        }
        bases[0].setRight(parameters[0]);
        vector.differentiate();
    }
}
