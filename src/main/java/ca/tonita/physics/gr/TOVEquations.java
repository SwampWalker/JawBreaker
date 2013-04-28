package ca.tonita.physics.gr;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;
import ca.tonita.math.numerical.QuasiLinearFirstOrderODESystem;
import ca.tonita.math.numerical.spectral.SpectralVector1D;

/**
 *
 * @author atonita
 */
public class TOVEquations implements QuasiLinearFirstOrderODESystem, NonLinearFirstOrderODESystem {

    private TabulatedHermite eos;
    private NonLinearFirstOrderODEBean bean;
    private int iPressure = 0;
    private int iMass = 1;
    private int iLambda = 2;
    private int nEquations = 3;

    TOVEquations(TabulatedHermite eos) {
        this.eos = eos;
        bean = new NonLinearFirstOrderODEBean();
        bean.setResidue(new double[3]);
        bean.setJacobian(new double[3][6]);
    }

    public TabulatedHermite getEos() {
        return eos;
    }

    public void setEos(TabulatedHermite eos) {
        this.eos = eos;
    }

    public double[] rightHandSide(double r, double[] y) {
        return new double[]{dpdr(r, y), dmdr(r, y), dlambdadr(r, y)};
    }

    /**
     * Returns the derivative of mass with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of mass wrt r
     */
    public double dmdr(double r, double[] y) {
        double rho = eos.energyDensity(y[0]);
        return 4 * Math.PI * rho * r * r;
    }

    /**
     * Returns the derivative of lambda with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of lambda wrt r
     */
    private double dlambdadr(double r, double[] y) {
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
    private double dpdr(double r, double[] y) {
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
                jacobian[i][residue.length + i] = 1;
            }
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
            double d2lambdadrdp = (mass + 4 * Math.PI * Math.pow(r, 3))
                    / (r * (r - 2 * mass));
            double d2lambdadrdm = 2 * (r - mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                    / (r * (r - 2 * mass) * (r - 2 * mass));
            // Pressure equation.
            jacobian[iPressure][iPressure] = (1 + drho) * dlambdadr + (pressure + rho) * d2lambdadrdp;
            jacobian[iPressure][iMass] = (pressure + rho) * d2lambdadrdm;
            jacobian[iPressure][iLambda] = 0;
            jacobian[iPressure][nEquations + iPressure] = 1;
            jacobian[iPressure][nEquations + iMass] = 0;
            jacobian[iPressure][nEquations + iLambda] = 0;
            jacobian[iPressure][2 * nEquations] = 0;
            // Mass equation.
            jacobian[iMass][iPressure] = - 4 * Math.PI * drho * r * r;
            jacobian[iMass][iMass] = 0;
            jacobian[iMass][iLambda] = 0;
            jacobian[iMass][nEquations + iPressure] = 0;
            jacobian[iMass][nEquations + iMass] = 1;
            jacobian[iMass][nEquations + iLambda] = 0;
            jacobian[iMass][2 * nEquations] = 0;
            // Lambda equation.
            jacobian[iLambda][iPressure] = -d2lambdadrdp;
            jacobian[iLambda][iMass] = -d2lambdadrdm;
            jacobian[iLambda][iLambda] = 0;
            jacobian[iLambda][nEquations + iPressure] = 0;
            jacobian[iLambda][nEquations + iMass] = 0;
            jacobian[iLambda][nEquations + iLambda] = 1;
            jacobian[iLambda][2 * nEquations] = 0;
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
     * @param dconstraint The vector that dotted with a perturbation gives the derivative of the constraint in the direction of that perturbation. I.E. dconstraint.
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
}
