/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.datasets;

import java.util.ArrayList;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import org.jfree.data.DomainInfo;
import org.jfree.data.Range;
import org.jfree.data.RangeInfo;
import org.jfree.data.xy.AbstractXYDataset;

/**
 *
 * @author atonita
 */
public class EOSDataset extends AbstractXYDataset implements DomainInfo, RangeInfo {

    private ArrayList<double[][]> dataSets = null;
    private ArrayList<Boolean> activated = null;
    private ArrayList<String> names = null;
    private int iX = 1;
    private int iY = 2;

    /**
     * Sets the variable being used for the range.
     *
     * @param iY the index of the variable to use for the range
     */
    public void setRangeVariable(int iY) {
        if (iY < variableNames.length) {
            this.iY = iY;
        } else {
            throw new IllegalArgumentException("Range variable index (" + iY + ") larger than number of variables supported");
        }
    }

    /**
     * Sets the variable being used for the domain.
     *
     * @param iX the index of the variable to use for the domain
     */
    public void setDomainVariable(int iX) {
        if (iX < variableNames.length) {
            this.iX = iX;
        } else {
            throw new IllegalArgumentException("Domain variable index (" + iX + ") larger than number of variables supported");
        }
    }

    public EOSDataset() {
        dataSets = new ArrayList<double[][]>();
        activated = new ArrayList<Boolean>();
        names = new ArrayList<String>();
    }
    
    private String[] variableNames = {"Number Density", "Pressure", "Total Energy Density"};

    public String[] getVariableNames() {
        return variableNames;
    }

    /**
     * Returns the name of the Range variable.
     *
     * @return The name of the data being used for the range.
     */
    public String getRangeName() {
        return variableNames[iY];
    }

    /**
     * Returns the name of the domain variable.
     *
     * @return The name of the data being used for the domain.
     */
    public String getDomainName() {
        return variableNames[iX];
    }

    /**
     * Adds an equation of state to this data set.
     *
     * @param index the index into the datasets to add the table
     * @param eos the equation of state to add.
     */
    public void add(int index, TabulatedHermite eos) {
        double[][] table = new double[3][];
        eos.cloneTable(table);
        dataSets.add(index, table);
        activated.add(index, Boolean.TRUE);
        names.add(eos.getIdentifier());
    }
    
    /**
     * Gets the index into dataSets given the series number. Necessary since
     * not all dataSets will be activated.
     * 
     * @param series the number of the series
     * @return the index into the dataSets.
     */
    private int getDatasetIndex(int series) {
        int index = -1;
        int count = -1;
        while (count < series) {
            index++;
            if (activated.get(index)) {
                count++;
            }
        }
        return index;
    }

    @Override
    public int getSeriesCount() {
        int count = 0;
        for (Boolean active: activated) {
            if (active) {
                count++;
            }
        }
        return count;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        return names.get(getDatasetIndex(series));
    }

    public int getItemCount(int series) {
        return this.dataSets.get(getDatasetIndex(series))[0].length;
    }

    public Number getX(int series, int item) {
        return dataSets.get(getDatasetIndex(series))[iX][item];
    }

    public Number getY(int series, int item) {
        return dataSets.get(getDatasetIndex(series))[iY][item];
    }

    public double getDomainLowerBound(boolean includeInterval) {
        double lb = 0;
        if (!dataSets.isEmpty()) {
            lb = Double.MAX_VALUE;
            for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
                if (activated.get(iDataSet) && dataSets.get(iDataSet)[iX][0] < lb) {
                    lb = dataSets.get(iDataSet)[iX][0];
                }
            }
        }
        return lb;
    }

    public double getDomainUpperBound(boolean includeInterval) {
        double ub = 1;
        if (!dataSets.isEmpty()) {
            ub = Double.MIN_VALUE;
            for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
                int last = dataSets.get(iDataSet)[iX].length - 1;
                if (activated.get(iDataSet) && dataSets.get(iDataSet)[iX][last] > ub) {
                    ub = dataSets.get(iDataSet)[iX][last];
                }
            }
        }
        return ub;
    }

    public Range getDomainBounds(boolean includeInterval) {
        return new Range(getDomainLowerBound(includeInterval), getDomainUpperBound(includeInterval));
    }

    public double getRangeLowerBound(boolean includeInterval) {
        if (!dataSets.isEmpty()) {
            double lb = Double.MAX_VALUE;
            for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
                if (activated.get(iDataSet) && dataSets.get(iDataSet)[iY][0] < lb) {
                    lb = dataSets.get(iDataSet)[iY][0];
                }
            }
            return lb;
        }
        return 0;
    }

    public double getRangeUpperBound(boolean includeInterval) {
        if (!dataSets.isEmpty()) {
            double ub = Double.MIN_VALUE;
            for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
                int last = dataSets.get(iDataSet)[iY].length - 1;
                if (activated.get(iDataSet) && dataSets.get(iDataSet)[iY][last] > ub) {
                    ub = dataSets.get(iDataSet)[iY][last];
                }
            }
            return ub;
        }
        return 1;
    }

    public Range getRangeBounds(boolean includeInterval) {
        return new Range(getRangeLowerBound(includeInterval), getRangeUpperBound(includeInterval));
    }
    
    /**
     * Sets the series as activated or not.
     * 
     * @param index The index of the eos to activate/disactivate.
     * @param activated the state to set.
     */
    public void setActivated(int index, boolean activated) {
        this.activated.set(index, activated);
    }

    /**
     * Removes the indexed dataset.
     * @param index The index of the dataset to remove. 
     */
    public void remove(int index) {
        dataSets.remove(index);
        activated.remove(index);
        names.remove(index);
    }
}
