package ca.tonita.math.numerical;

/**
 * Interface for a Nonlinear ODE system.
 * @author atonita
 */
public interface NonLinearFirstOrderODESystem {
    
    static int BULK = 0;
    static int LEFTBOUNDARY = 1;
    static int RIGHTBOUNDARY = 2;

    public NonLinearFirstOrderODEBean equations(double r, double[] y, double[] dy, double[] parameters, int type);
    
    public void computeConstraints(double[] x, double[][] variables, double[][] dvariables, double[] parameters, double[] dconstraint, int iConstraint);
}
