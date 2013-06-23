package ca.tonita.jawbreaker.models;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class TOVFamily extends ArrayList<TOVData> {
    private TabulatedHermite eos;
    
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
}
