package jaw.breaker.tov;

import jaw.breaker.equationsOfState.TabulatedHermite;
import static java.lang.Math.PI;

/**
 *
 * @author atonita
 */
public class TOVRK4 {

    /**
     * Computes the TOV solution using standard Runge-Kutta 4th order.
     * @param eos The equation of state.
     * @param centralPressure The central pressure of the solution.
     * @param terminationPressure The pressure to terminate integration at.
     * @param h The step size.
     * @return An array containing the mass and radius: {m, r}
     */
    public double[] getSolution(TabulatedHermite eos, double centralPressure, double terminationPressure, double h) {
        double w[] = {centralPressure, 0}; // Solution at the center.
        double wOld[] = {w[0], w[1]};
        double r = 0;
        double rNext = 0.1;
        int i = 0;
        while (w[0] > terminationPressure) {
            double K1[] = {h * (dpdr(r, w[0], w[1], eos)),
                h * (dmdr(r, w[0], w[1], eos))};
            double K2[] = {h * dpdr(r + h * 0.5, w[0] + K1[0] * 0.5, w[1] + K1[1] * 0.5, eos),
                h * dmdr(r + h * 0.5, w[0] + K1[0] * 0.5, w[1] + K1[1] * 0.5, eos)};
            double K3[] = {h * dpdr(r + h * 0.5, w[0] + K2[0] * 0.5, w[1] + K2[1] * 0.5, eos),
                h * dmdr(r + h * 0.5, w[0] + K2[0] * 0.5, w[1] + K2[1] * 0.5, eos)};
            double K4[] = {h * dpdr(r + h, w[0] + K3[0], w[1] + K3[1], eos),
                h * dmdr(r + h, w[0] + K3[0], w[1] + K3[1], eos)};
            wOld[0] = w[0];
            wOld[1] = w[1];
            w[0] += (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]) / 6.;
            w[1] += (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]) / 6.;
            r = r + h;
            i++;
        }
        return new double[]{w[1],r};
    }

    /**
     * The TOV pressure equation for dp/dr.
     * @param r The radius.
     * @param p The pressure.
     * @param m The mass potential.
     * @param eos The equation of state.
     * @return The value of dp/dr that satisfies the TOV equation for those values.
     */
    private double dpdr(double r, double p, double m, TabulatedHermite eos) {
        if (r == 0.) {
            return 0;
        }
        double dlambdadr = (m + 4 * PI * r * r * r * p)
                / (r * (r - 2 * m));
        double rho = eos.energyDensity(p);
        return -(p + rho) * dlambdadr;
    }

    /**
     * The TOV pressure equation for dm/dr.
     * @param r The radius.
     * @param p The pressure.
     * @param m The mass potential.
     * @param eos The equation of state.
     * @return The value of dm/dr that satisfies the TOV equation for those values.
     */
    private double dmdr(double r, double p, double m, TabulatedHermite eos) {
        if (r == 0.) {
            return 0;
        }
        double rho = eos.energyDensity(p);
        return 4 * PI * rho * r * r;
    }
}
