package ca.tonita.jawbreaker.models;

import ca.tonita.physics.gr.hydro.TOVData;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class TOVFamily extends ArrayList<TOVData> {
    private TabulatedHermite eos;
    private int iMaximum;
    
    /**
     * Constructs the TOVFamily without any data.
     * @param eos the equation of state that will identify this family.
     */
    public TOVFamily(TabulatedHermite eos) {
        this.eos = eos;
    }
    
    @Override
    public String toString() {
        return eos.getIdentifier();
    }

    public TabulatedHermite getEos() {
        return eos;
    }

    /**
     * Adds the TOV solution to this family as the maximum.
     * @param iMaximum The index of the maximum in the family.
     * @param max The solution to use as the maximum.
     */
    public void setMaximum(int iMaximum, TOVData max) {
        this.add(iMaximum, max);
        this.iMaximum = iMaximum;
    }
    
    /**
     * Returns the TOV with maximum rest mass.
     * @return the Maximum rest mass TOV for the equation of state.
     */
    public TOVData getMaximum() {
        return get(iMaximum);
    }
    
    /**
     * Returns the index of the maximum.
     * @return the index of the maximum mass TOV
     */
    public int getIMaximum() {
        return iMaximum;
    }
}
