package ca.tonita.jawbreaker.eoswindows;

import javax.swing.event.ChangeListener;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;

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
    
    /**
     * Returns the number of equations of state stored.
     * @return the number of equations of state stored.
     */
    public int nEOS();
    
    /**
     * Returns the names of the equations of state.
     * @return 
     */
    public String[] getEOSNames();
    
    /**
     * Returns the indexed equation of state.
     * @param i the index of the equation of state to get
     * @return the indexed equation of state
     */
    public TabulatedHermite getEos(int i);
    
    /**
     * Receive notifications of additional equations of state.
     * @param cl the change listener.
     */
    public void addEOSChangeListener(ChangeListener cl);
    
    /**
     * Removes the indexed equation of state.
     * @param i the index of the equation of state to remove.
     */
    public void removeEOS(int i);
}
