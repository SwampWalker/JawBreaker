/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.datasets;

import ca.tonita.jawbreaker.models.JawBreakerModel;
import ca.tonita.jawbreaker.models.TOVFamily;

/**
 *
 * @author atonita
 */
public class TOVFamilyDataset extends MultiDataset<TOVFamily> {
    
    private JawBreakerModel model;

    public TOVFamilyDataset(JawBreakerModel model) {
        variableNames = tovFamilyNames;
    }
    
    private static final String[] tovFamilyNames = {"Radius", "Central pressure", "Gravitational mass", "Rest mass", "Core radius", "Core rest mass"};

    /**
     * Returns the specified data element.
     * @param iItem the index of the point to get the variable of
     * @param iVariable the variable to get, as specified by the variableNames field
     * @return the variable at the specified point
     */
    protected double getItem(int iSeries, int iVariable, int iItem) {
        double variable;
        if (iVariable == 0) {
            variable = getDataset(iSeries).get(iItem).getRadius();
        } else if (iVariable == 1) {
            variable = getDataset(iSeries).get(iItem).getPressure(0);
        } else if (iVariable == 2) {
            variable = getDataset(iSeries).get(iItem).getGravitationalMass();
        } else if (iVariable == 3) {
            variable = getDataset(iSeries).get(iItem).getConservedMass();
        } else if (iVariable == 4) {
            variable = getDataset(iSeries).get(iItem).getCoreRadius();
        } else if (iVariable == 5) {
            variable = getDataset(iSeries).get(iItem).getCoreRestMass();
        } else {
            throw new IllegalArgumentException("No variable with index " + iVariable + " in TOVFamilyDataset");
        }
        return variable;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        return datasets.get(getDatasetIndex(series)).toString();
    }

    public int getItemCount(int series) {
        return datasets.get(getDatasetIndex(series)).size();
    }
    
    public TOVFamily get(int i) {
        return datasets.get(i);
    }
}
