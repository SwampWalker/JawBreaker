package ca.tonita.jawbreaker.datasets;

import java.util.ArrayList;
import org.jfree.data.DomainInfo;
import org.jfree.data.Range;
import org.jfree.data.RangeInfo;
import org.jfree.data.xy.AbstractXYDataset;

/**
 * A MultiDataset is a collection of series, each of which has multiple
 * variables, and each variable has multiple items. For instance, a collection
 * of Equation of State series would have the variables pressure, specific
 * energy, mass density, temperature etc, and multiple points showing the
 * relations at various values of the dependent variable (depending on what one
 * thinks is the dependent variable).
 *
 * @author atonita
 */
public abstract class MultiDataset<T> extends AbstractXYDataset implements DomainInfo, RangeInfo {

    protected ArrayList<Boolean> activated = new ArrayList<Boolean>();
    protected ArrayList<T> datasets = new ArrayList<T>();
    protected int iX = 0;
    protected int iY = 1;
    protected String[] variableNames = null;

    /**
     * Adds a dataset to this MultiDataset.
     *
     * @param index The index to add this data set at.
     * @param dataset The data set to add.
     */
    public void add(int index, T dataset) {
        datasets.add(index, dataset);
        System.out.println("Setting activated.");
        activated.add(index, Boolean.TRUE);
    }

    /**
     * Removes the indexed dataset.
     *
     * @param index the index of the dataset to remove
     * @return The removed dataset.
     */
    public T remove(int index) {
        activated.remove(index);
        return datasets.remove(index);
    }

    /**
     * Returns the dataset given the series index. If some datasets are
     * disactivated, then there are more datasets than series.
     *
     * @param iSeries the index of the series to get
     * @return the indexed dataset
     */
    public T getDataset(int iSeries) {
        return datasets.get(getDatasetIndex(iSeries));
    }

    /**
     * Should return the i'th element of variable indexed by iVariable
     *
     * @param iSeries the index of the data series to get the element from
     * @param iVariable the index of the variable of the data series to get
     * @param iItem the index of the particular item to get the variable of
     * @return the variable specified
     */
    protected abstract double getItem(int iSeries, int iVariable, int iItem);

    /**
     * Gets the index into dataSets given the series number. Necessary since not
     * all dataSets will be activated.
     *
     * @param series the number of the series
     * @return the index into the dataSets.
     */
    protected int getDatasetIndex(int series) {
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

    public Range getDomainBounds(boolean includeInterval) {
        return new Range(getDomainLowerBound(includeInterval), getDomainUpperBound(includeInterval));
    }

    public double getDomainLowerBound(boolean includeInterval) {
        double lb = 0;
        int nSeries = getSeriesCount();
        if (nSeries != 0) {
            lb = Double.MAX_VALUE;
            for (int iDataSet = 0; iDataSet < nSeries; iDataSet++) {
                if (activated.get(iDataSet) && getItem(iDataSet, iX, 0) < lb) {
                    lb = getItem(iDataSet, iX, 0);
                }
            }
        }
        return lb;
    }

    /**
     * Returns the name of the domain variable.
     *
     * @return The name of the data being used for the domain.
     */
    public String getDomainName() {
        return variableNames[iX];
    }

    public double getDomainUpperBound(boolean includeInterval) {
        double ub = 1;
        int nSeries = getSeriesCount();
        if (nSeries != 0) {
            ub = Double.MIN_VALUE;
            for (int iSeries = 0; iSeries < nSeries; iSeries++) {
                int last = getItemCount(iSeries) - 1;
                if (activated.get(iSeries) && getItem(iSeries, iX, last) > ub) {
                    ub = getItem(iSeries, iX, last);
                }
            }
        }
        return ub;
    }

    public Range getRangeBounds(boolean includeInterval) {
        return new Range(getRangeLowerBound(includeInterval), getRangeUpperBound(includeInterval));
    }

    public double getRangeLowerBound(boolean includeInterval) {
        int nSeries = getSeriesCount();
        if (nSeries != 0) {
            double lb = Double.MAX_VALUE;
            for (int iDataSet = 0; iDataSet < nSeries; iDataSet++) {
                if (activated.get(iDataSet) && getItem(iDataSet, iY, 0) < lb) {
                    lb = getItem(iDataSet, iY, 0);
                }
            }
            return lb;
        }
        return 0;
    }

    /**
     * Returns the name of the Range variable.
     *
     * @return The name of the data being used for the range.
     */
    public String getRangeName() {
        return variableNames[iY];
    }

    public double getRangeUpperBound(boolean includeInterval) {
        int nSeries = getSeriesCount();
        if (nSeries != 0) {
            double ub = Double.MIN_VALUE;
            for (int iSeries = 0; iSeries < nSeries; iSeries++) {
                int last = getItemCount(iSeries) - 1;
                if (activated.get(iSeries) && getItem(iSeries, iX, last) > ub) {
                    ub = getItem(iSeries, iX, last);
                }
            }
            return ub;
        }
        return 1;
    }

    @Override
    public int getSeriesCount() {
        int count = 0;
        for (Boolean active : activated) {
            if (active) {
                count++;
            }
        }
        return count;
    }

    public String[] getVariableNames() {
        return variableNames;
    }

    public Number getX(int series, int item) {
        return getItem(series, iX, item);
    }

    public Number getY(int series, int item) {
        return getItem(series, iY, item);
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
}
