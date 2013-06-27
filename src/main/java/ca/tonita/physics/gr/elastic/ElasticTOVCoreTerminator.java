package ca.tonita.physics.gr.elastic;

import ca.tonita.math.numerical.EvolutionTerminator;

/**
 * Terminates TOV integration when a maximum radius, maximum number of steps or
 * minimum pressure is reached.
 *
 * @author atonita
 */
public class ElasticTOVCoreTerminator extends EvolutionTerminator {

    protected double coreMass;
    private final double minPressure;

    private ElasticTOVCoreTerminator(double tMax, int maxSteps, double minPressure) {
        super(tMax, maxSteps);
        this.minPressure = minPressure;
    }

    public ElasticTOVCoreTerminator(double maxRadius, int maxSteps, double coreMass, double minPressure) {
        super(maxRadius, maxSteps);
        this.coreMass = coreMass;
        this.minPressure = minPressure;
    }

    @Override
    public boolean terminate(double[] y, double t, int steps) {
        boolean terminate = super.terminate(y, t, steps);
        if (y[ElasticTOVIndex.RESTMASS] > coreMass) {
            terminate = true;
        } else if (y[ElasticTOVIndex.PRESSURE] < minPressure) {
            terminate = true;
        }
        return terminate;
    }
}
