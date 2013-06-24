/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.models;

import java.util.ArrayList;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import ca.tonita.jawbreaker.datasets.EOSDataset;
import ca.tonita.jawbreaker.datasets.TOVFamilyDataset;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;

/**
 *
 * @author atonita
 */
public class JawBreakerModel implements EOSStorage {

    TOVFamilyDataset tovFamilies;
    ArrayList<TabulatedHermite> eosStorage;
    ArrayList<ChangeListener> eosChangeListeners;
    EOSDataset eosDataset;

    /**
     * Constructs the data model.
     */
    public JawBreakerModel() {
        eosStorage = new ArrayList<TabulatedHermite>();
        eosChangeListeners = new ArrayList<ChangeListener>();
        tovFamilies = new TOVFamilyDataset(this);
    }

    /**
     * Adds an equation of state. The equations of state are sorted by name.
     * Duplicate names are not added.
     *
     * @param eos the equation of state to add.
     */
    public void add(TabulatedHermite eos) {
        // Order eos by identifier, ignoring case.
        String id = eos.getIdentifier();
        String comparator = id.toLowerCase();
        int index = 0;
        boolean duplicate = false;
        for (int i = 0; i < eosStorage.size() && !duplicate; i++) {
            int compareTo = eosStorage.get(i).getIdentifier().toLowerCase().compareTo(comparator);
            if (compareTo == 0 && eosStorage.get(i).getIdentifier().compareTo(id) == 0) {
                duplicate = true;
            } else if (compareTo <= 0) {
                index = i + 1;
            } // Could terminate after we get > 0...
        }

        // Add eos.
        if (!duplicate) {
            eosStorage.add(index, eos);
            if (eosDataset != null) {
                eosDataset.add(index, eos);
            }
            if (tovFamilies != null) {
                System.out.println("Adding TOV Family.");
                tovFamilies.add(index, new TOVFamily(eos));
            }
        }

        fireEOSChange();
    }

    /**
     * Returns the number of equations of state.
     *
     * @return the number of equations of state.
     */
    public int nEOS() {
        return eosStorage.size();
    }

    /**
     * Returns the names of the equations of state.
     *
     * @return the names of the equstions of state.
     */
    public String[] getEOSNames() {
        String[] names = new String[eosStorage.size()];
        for (int i = 0; i < eosStorage.size(); i++) {
            names[i] = eosStorage.get(i).getIdentifier();
        }
        return names;
    }

    /**
     * Returns the equation of state by the index.
     *
     * @param i the index of the equation of state to get.
     * @return the indexed equation of state.
     */
    public TabulatedHermite getEos(int i) {
        return eosStorage.get(i);
    }

    /**
     * Adds a changelistener that will be informed of additions or deletions to
     * the equation of state storage.
     *
     * @param cl the change listener to add
     */
    public void addEOSChangeListener(ChangeListener cl) {
        eosChangeListeners.add(cl);
    }

    public void removeEOS(int i) {
        TabulatedHermite toRemove = eosStorage.remove(i);
        if (eosDataset != null) {
            eosDataset.remove(i);
        }
        if (tovFamilies != null) {
            tovFamilies.remove(i);
        }
        fireEOSChange();
    }

    /**
     * Notifies all eosChangeListeners of a change.
     */
    private void fireEOSChange() {
        for (ChangeListener cl : eosChangeListeners) {
            cl.stateChanged(new ChangeEvent(this));
        }
    }

    public EOSDataset getEOSDataset() {
        if (eosDataset == null) {
            eosDataset = new EOSDataset();
            for (int i = 0; i < eosStorage.size(); i++) {
                eosDataset.add(i, eosStorage.get(i));
            }
        }
        return eosDataset;
    }

    public TOVFamilyDataset getTOVFamilies() {
        return tovFamilies;
    }

    public TOVFamily getTOVFamily(int iFamily) {
        return tovFamilies.get(iFamily);
    }

    public void setEOSDataset(EOSDataset eosDataset) {
        this.eosDataset = eosDataset;
    }
}
