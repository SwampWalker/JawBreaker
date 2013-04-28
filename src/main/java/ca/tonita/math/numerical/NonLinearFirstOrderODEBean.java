package ca.tonita.math.numerical;

/**
 * Stores the data for a non-linear ODE system. <br> The residue is given by:
 * f(y,dy) = 0<br> The jacobian is given by J(y,dy) = {Df/Dy, Df/Ddy}.
 *
 * @author atonita
 */
public class NonLinearFirstOrderODEBean {

    private double[] residue;
    private double[][] jacobian;

    /**
     * The residue R_i would be the residual of the i'th equation.
     *
     * @return the residue.
     */
    public double[] getResidue() {
        return residue;
    }

    /**
     * Sets the value of the residue.
     *
     * @param residue
     */
    public void setResidue(double[] residue) {
        this.residue = residue;
    }

    /**
     * The Jacobian. J_{i,j} is the derivative of the i'th equation with respect
     * to the j'th variable. This should roll over and become derivative of the
     * i'th equation with respect to the derivative of the (j-nVar)'th variable.
     * After that, it should be derivatives with respect to parameters.
     *
     * @return The jacobian.
     */
    public double[][] getJacobian() {
        return jacobian;
    }

    /**
     * Sets the jacobian.
     * @param jacobian The jacobian to set.
     */
    public void setJacobian(double[][] jacobian) {
        this.jacobian = jacobian;
    }
}
