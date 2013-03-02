/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.datasets;

import java.util.ArrayList;
import jaw.breaker.equationsOfState.TabulatedHermite;
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
    private int iX = 1;
    private int iY = 2;

    /**
     * Sets the variable being used for the range.
     *
     * @param iY the index of the variable to use for the range
     */
    public void setRangeVariable(int iY) {
        if (iY < dataNames.length) {
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
        if (iX < dataNames.length) {
            this.iX = iX;
        } else {
            throw new IllegalArgumentException("Domain variable index (" + iX + ") larger than number of variables supported");
        }
    }

    public EOSDataset() {
        dataSets = new ArrayList<double[][]>();
    }
    private String[] dataNames = {"Number Density", "Pressure", "Total Energy Density"};

    public String[] getVariableNames() {
        return dataNames;
    }

    /**
     * Returns the name of the Range variable.
     *
     * @return The name of the data being used for the range.
     */
    public String getRangeName() {
        return dataNames[iY];
    }

    /**
     * Returns the name of the domain variable.
     *
     * @return The name of the data being used for the domain.
     */
    public String getDomainName() {
        return dataNames[iX];
    }

    /**
     * Adds an equation of state to this data set.
     *
     * @param eos the equation of state to add.
     */
    public void add(TabulatedHermite eos) {
        double[][] table = new double[3][];
        eos.cloneTable(table);
        dataSets.add(table);
    }

    @Override
    public int getSeriesCount() {
        return dataSets.size();
    }

    @Override
    public Comparable getSeriesKey(int series) {
        return "EOS " + series;
    }

    public int getItemCount(int series) {
        return this.dataSets.get(series)[0].length;
    }

    public Number getX(int series, int item) {
        return dataSets.get(series)[iX][item];
    }

    public Number getY(int series, int item) {
        return dataSets.get(series)[iY][item];
    }

    public double getDomainLowerBound(boolean includeInterval) {
        double lb = 0;
        if (!dataSets.isEmpty()) {
            lb = Double.MAX_VALUE;
            for (int series = 0; series < dataSets.size(); series++) {
                if (dataSets.get(series)[iX][0] < lb) {
                    lb = dataSets.get(series)[iX][0];
                }
            }
        }
        return lb;
    }

    public double getDomainUpperBound(boolean includeInterval) {
        double ub = 1;
        if (!dataSets.isEmpty()) {
            ub = Double.MIN_VALUE;
            for (int series = 0; series < dataSets.size(); series++) {
                int last = dataSets.get(series)[iX].length - 1;
                if (dataSets.get(series)[iX][last] > ub) {
                    ub = dataSets.get(series)[iX][last];
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
            for (int series = 0; series < dataSets.size(); series++) {
                if (dataSets.get(series)[iY][0] < lb) {
                    lb = dataSets.get(series)[iY][0];
                }
            }
            return lb;
        }
        return 0;
    }

    public double getRangeUpperBound(boolean includeInterval) {
        if (!dataSets.isEmpty()) {
            double ub = Double.MIN_VALUE;
            for (int series = 0; series < dataSets.size(); series++) {
                int last = dataSets.get(series)[iY].length - 1;
                if (dataSets.get(series)[iY][last] > ub) {
                    ub = dataSets.get(series)[iY][last];
                }
            }
            return ub;
        }
        return 1;
    }

    public Range getRangeBounds(boolean includeInterval) {
        return new Range(getRangeLowerBound(includeInterval), getRangeUpperBound(includeInterval));
    }
}
