package ca.tonita.math.numerical;

import ca.tonita.math.numerical.spectral.SpectralVector1D;
import ca.tonita.math.polynomials.LinearlyMappedBasis;

/**
 *
 * @author atonita
 */
public interface NewtonRaphsonPostProcessor {

    /**
     * A function called after a Newton-Raphson step to do post processing.
     * @param vector The data vector used in the NR step.
     * @param bases The polynomical bases used to solve the problem.
     */
    void postStepProcessing(SpectralVector1D vector, LinearlyMappedBasis[] bases);
    
}
