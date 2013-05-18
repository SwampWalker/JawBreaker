package ca.tonita.math.numerical;

import ca.tonita.math.numerical.spectral.SpectralVector1D;

/**
 * Interface for a Nonlinear ODE system.
 *
 * @author atonita
 */
public interface NonLinearFirstOrderODESystem {
    
    public static int BULK = 0;
    public static int LEFTBOUNDARY = 1;
    public static int RIGHTBOUNDARY = 2;

    /**
     * Should return a filled
     * <code>NonLinearFirstOrderODEBean</code> with the residual and jacobian.
     *
     * @param iX The coordinate index.
     * @param x The coordinate value.
     * @param y The variables.
     * @param dy The derivative of the variables.
     * @param parameters The parameters.
     * @param type The type of the equation, one of bulk or boundary. Use the static members of this class.
     * @return the residual and jacobian
     */
    public NonLinearFirstOrderODEBean equations(int iX, double x, double[] y, double[] dy, double[] parameters, int type);

    /**
     * Computes and returns the constraint, as well as fills the jacobian of the
     * constraint with respect to the variables.
     *
     * @return
     */
    public double computeConstraint(int iConstraint, ODEIndexer1D indexer, SpectralVector1D vector, double[] dconstraint);
    
    /**
     * Returns the number of domains of the problem.
     * @return the number of domains.
     */
    public int getNDomains();
    
    /**
     * Returns the number of variables per domain.
     * @return the number of variables, indexed per domain
     */
    public int[] getNVariables();
    
    /**
     * Returns the number of variables for the specified domain.
     * @param iDomain The domain to get the number of variables in.
     * @return the number of variables.
     */
    public int getNVariables(int iDomain);
    
    /**
     * Returns the number of parameters of the problem.
     * @return the number of parameters.
     */
    public int getNParameters();
}
