/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.datasets;

import java.util.ArrayList;
import jaw.breaker.models.TOVData;
import org.jfree.data.DomainInfo;
import org.jfree.data.Range;
import org.jfree.data.RangeInfo;
import org.jfree.data.xy.AbstractXYDataset;

/**
 *
 * @author atonita
 */
public class TOVDataset extends AbstractXYDataset implements DomainInfo, RangeInfo {

    private TOVData tov = null;

    public TOVData getTov() {
        return tov;
    }

    public void setTOV(TOVData tov) {
        this.tov = tov;
    }
    private ArrayList<String> names = null;
    private int iX = 0;
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

    public TOVDataset() {
        names = new ArrayList<String>();
    }
    private String[] variableNames = {"Radius", "Number Density", "Pressure", "Total Energy Density"};

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

    @Override
    public int getSeriesCount() {
        if (tov != null) {
            return 1;
        }
        return 0;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        if (series == 0 && tov != null) {
            return tov.getIdentifier();
        }
        return null;
    }

    public int getItemCount(int series) {
        if (series == 0) {
            return tov.getNPoints();
        }
        return 0;
    }

    public Number getX(int series, int item) {
        if (series == 0) {
            return tov.getVariables(item)[iX];
        }
        return null;
    }

    public Number getY(int series, int item) {
        if (series == 0) {
            return tov.getVariables(item)[iY];
        }
        return null;
    }

    public double getDomainLowerBound(boolean includeInterval) {
        double lb = 0.0;
        if (tov != null) {
            lb = tov.getVariables(0)[iX];
        }
        return lb;
    }

    public double getDomainUpperBound(boolean includeInterval) {
        double ub = 1;
        if (tov != null) {
            ub = tov.getVariables(tov.getNPoints() - 1)[iX];
        }
        return ub;
    }

    public Range getDomainBounds(boolean includeInterval) {
        return new Range(getDomainLowerBound(includeInterval), getDomainUpperBound(includeInterval));
    }

    public double getRangeLowerBound(boolean includeInterval) {
        double lb = 0.0;
        if (tov != null) {
            lb = tov.getVariables(0)[iY];
        }
        return lb;
    }

    public double getRangeUpperBound(boolean includeInterval) {
        double ub = 1;
        if (tov != null) {
            ub = tov.getVariables(tov.getNPoints() - 1)[iY];
        }
        return ub;
    }

    public Range getRangeBounds(boolean includeInterval) {
        return new Range(getRangeLowerBound(includeInterval), getRangeUpperBound(includeInterval));
    }
}
