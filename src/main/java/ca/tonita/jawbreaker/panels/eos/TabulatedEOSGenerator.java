/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.panels.eos;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;

/**
 *
 * @author tonita
 */
public interface TabulatedEOSGenerator {
    /**
     * Gets the configured equation of state.
     * @return a <code>TabulatedHermite</code> equation of state.
     */
    public TabulatedHermite getTabulatedEOS();
};
