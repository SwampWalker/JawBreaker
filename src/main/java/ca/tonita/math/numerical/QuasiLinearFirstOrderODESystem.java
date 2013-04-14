package ca.tonita.math.numerical;

/**
 * A quasi-linear ODE system is any differential equation of the form:<br>
 * dy/dt = f(t,y)<br>
 * Where y can be multi-dimensional.
 * 
 * @author atonita
 */
public interface QuasiLinearFirstOrderODESystem {
    /**
     * Computes the right hand side of the system.
     * @param t the independent variable
     * @param y the dependent variable(s)
     * @return the right hand side, f
     */
    public double[] rightHandSide(double t, double[] y);
}
