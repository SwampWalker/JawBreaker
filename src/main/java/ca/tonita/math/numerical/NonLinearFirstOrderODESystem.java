package ca.tonita.math.numerical;

import ca.tonita.math.numerical.spectral.SpectralVector1D;

/**
 * Interface for a Nonlinear ODE system.
 *
 * @author atonita
 */
public interface NonLinearFirstOrderODESystem {

    /**
     * Should return a filled
     * <code>NonLinearFirstOrderODEBean</code> with the residual and jacobian.
     *
     * @param iX The coordinate index.
     * @param x The coordinate value.
     * @param y The variables.
     * @param dy The derivative of the variables.
     * @param parameters The parameters.
     * @return the residual and jacobian
     */
    public NonLinearFirstOrderODEBean equations(int iX, double x, double[] y, double[] dy, double[] parameters);

    /**
     * Computes and returns the constraint, as well as fills the jacobian of the
     * constraint with respect to the variables.
     *
     * @return
     */
    public double computeConstraint(int iConstraint, ODEIndexer1D indexer, SpectralVector1D vector, double[] dconstraint);
}
