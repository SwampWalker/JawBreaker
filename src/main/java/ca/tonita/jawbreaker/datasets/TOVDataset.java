/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.datasets;

import ca.tonita.physics.gr.hydro.TOVData;

/**
 *
 * @author atonita
 */
public class TOVDataset extends MultiDataset<TOVData> {

    public TOVData getTov(int i) {
        return datasets.get(i);
    }

    public TOVDataset() {
        variableNames = tovVariableNames;
    }
    
    private static final String[] tovVariableNames = {"Radius", "Pressure", "Gravitational Mass", "Lambda", "Interior mass", "Number Density", "Energy Density", "Mass Density"};

    /**
     * Returns the specified data element.
     * @param iItem the index of the point to get the variable of
     * @param iVariable the variable to get, as specified by the variableNames field
     * @return the variable at the specified point
     */
    protected double getItem(int iSeries, int iVariable, int iItem) {
        double variable;
        if (iVariable == 0) {
            variable = getDataset(iSeries).getRadius(iItem);
        } else if (iVariable < 5) {
            variable = getDataset(iSeries).getVariables(iItem)[iVariable - 1];
        } else {
            variable = getDataset(iSeries).getSecondaries(iItem)[iVariable - 4];
        }
        return variable;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        return datasets.get(getDatasetIndex(series)).getIdentifier();
    }

    public int getItemCount(int series) {
        return datasets.get(getDatasetIndex(series)).getNPoints();
    }
}
