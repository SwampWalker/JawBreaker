package jaw.breaker.eoswindows;

import jaw.breaker.equationsOfState.TabulatedHermite;

/**
 * Interface between the NewEOSDialog and the EOSPanel.
 * @author tonita
 */
public interface EOSStorage {
    /**
     * Add an equation of state to the storage.
     * @param eos the equation of state to add.
     */
    public void add(TabulatedHermite eos);
}
