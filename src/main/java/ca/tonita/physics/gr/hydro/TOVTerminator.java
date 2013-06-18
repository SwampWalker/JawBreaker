package ca.tonita.physics.gr.hydro;

import ca.tonita.math.numerical.EvolutionTerminator;

/**
 * Terminates TOV integration when a maximum radius, maximum number of steps or
 * minimum pressure is reached.
 *
 * @author atonita
 */
public class TOVTerminator extends EvolutionTerminator {

    protected double minPressure;

    private TOVTerminator(double tMax, int maxSteps) {
        super(tMax, maxSteps);
    }

    public TOVTerminator(double maxRadius, int maxSteps, double minPressure) {
        super(maxRadius, maxSteps);
        this.minPressure = minPressure;
    }

    @Override
    public boolean terminate(double[] y, double t, int steps) {
        boolean terminate = super.terminate(y, t, steps);
        if (y[0] < minPressure) {
            terminate = true;
        }
        return terminate;
    }
}
