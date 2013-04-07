/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.physics.gr;

import ca.tonita.math.numerical.RK4;
import java.util.ArrayList;
import jaw.breaker.equationsOfState.TabulatedHermite;
import jaw.breaker.models.TOVData;

/**
 *
 * @author atonita
 */
public class TOVBuilder {

    /**
     * No constructing this guy.
     */
    private TOVBuilder() {
    }

    /**
     * Evolves the TOV equations using RK4 method.
     *
     * @param tov a <code>TOVData</code> object to store the result
     * @param eos the equation of state to use
     * @param centralPressure the central pressure of the star
     * @param stepSize the step size to use for the RK4
     * @param outputEvery how often to output data
     * @param minPressure the value of pressure to terminate evolution at (the
     * effective surface pressure)
     */
    public static void evolve(TOVData tov, TabulatedHermite eos, double centralPressure, double stepSize, int outputEvery, double minPressure) {
        TOVEquations eqns = new TOVEquations(eos);
        ArrayList<double[]> points = tov.getVariables();
        ArrayList<Double> radii = tov.getRadii();
        points.clear();
        points.add(new double[]{centralPressure, 0, 0});
        radii.clear();
        radii.add(0.);
        int maxSteps = Integer.MAX_VALUE;
        TOVTerminator terminator = new TOVTerminator(maxSteps*stepSize, maxSteps, minPressure);
        RK4.evolve(points, radii, eqns, stepSize, outputEvery, terminator);
    }
}
