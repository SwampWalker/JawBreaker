package ca.tonita.math.numerical;

/**
 * Stores the data for a non-linear ODE system. <br>
 * The residue is given by: f(y,dy) = 0<br>
 * The jacobian is given by J(y,dy) = {Df/Dy, Df/Ddy}.
 * @author atonita
 */
public class NonLinearFirstOrderODEBean {
    private double[] residue;
    private double[][] jacobian;

    public double[] getResidue() {
        return residue;
    }

    public void setResidue(double[] residue) {
        this.residue = residue;
    }

    public double[][] getJacobian() {
        return jacobian;
    }

    public void setJacobian(double[][] jacobian) {
        this.jacobian = jacobian;
    }
}
