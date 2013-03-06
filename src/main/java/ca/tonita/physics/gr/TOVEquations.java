package ca.tonita.physics.gr;

import ca.tonita.math.numerical.QuasiLinearODESystem;
import jaw.breaker.equationsOfState.TabulatedHermite;

/**
 *
 * @author atonita
 */
public class TOVEquations implements QuasiLinearODESystem {

    private TabulatedHermite eos;

    TOVEquations(TabulatedHermite eos) {
        this.eos = eos;
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
        return  4 * Math.PI * rho * r * r;
    }

    /**
     * Returns the derivative of lambda with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param y The TOV variables {p,m,lambda}
     * @return the derivative of lambda wrt r
     */
    private double dlambdadr(double r, double[] y) {
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
        double mass = y[1];
        double pressure = y[0];
        double dlambdadr = (mass + 4 * Math.PI * Math.pow(r, 3) * pressure)
                / (r * (r - 2 * mass));
        double rho = eos.energyDensity(pressure);
        return -(pressure + rho) * dlambdadr;
    }
}
