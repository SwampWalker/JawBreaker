package ca.tonita.physics.gr.elastic;

import ca.tonita.math.numerical.EvolutionTerminator;

/**
 * Terminates TOV integration when a maximum radius, maximum number of steps or
 * minimum pressure is reached.
 *
 * @author atonita
 */
public class ElasticTOVCrustTerminator extends EvolutionTerminator {

    private double rMax;
    
    public ElasticTOVCrustTerminator(double rMax) {
        super(Double.MAX_VALUE, Integer.MAX_VALUE);
        this.rMax = rMax;
    }

    @Override
    public boolean terminate(double[] y, double t, int steps) {
        boolean terminate = super.terminate(y, t, steps);
        if (y[ElasticTOVIndex.XI] > rMax) {
            terminate = true;
        }
        return terminate;
    }
}
