/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.datasets;

import java.util.ArrayList;
import ca.tonita.jawbreaker.models.TOVData;
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
    private int iX = 0;
    private int iY = 1;

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
    }
    private static final String[] variableNames = {"Radius", "Pressure", "Gravitational Mass", "Lambda", "Number Density", "Energy Density"};

    public static String[] getVariableNames() {
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
        if (tov != null && tov.getNPoints() != 0) {
            return 1;
        }
        return 0;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        if (series == 0 && tov != null && tov.getNPoints() != 0) {
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
            return getVariable(item, iX);
        }
        return null;
    }

    public Number getY(int series, int item) {
        if (series == 0) {
            return getVariable(item, iY);
        }
        return null;
    }

    public double getDomainLowerBound(boolean includeInterval) {
        double lb = 0.0;
        if (tov != null && tov.getNPoints() != 0) {
            lb = getVariable(0, iX);
        }
        return lb;
    }

    public double getDomainUpperBound(boolean includeInterval) {
        double ub = 1;
        if (tov != null && tov.getNPoints() != 0) {
            ub = getVariable(tov.getNPoints()-1, iX);
        }
        return ub;
    }

    public Range getDomainBounds(boolean includeInterval) {
        return new Range(getDomainLowerBound(includeInterval), getDomainUpperBound(includeInterval));
    }

    public double getRangeLowerBound(boolean includeInterval) {
        double lb = 0.0;
        if (tov != null && tov.getNPoints() != 0) {
            lb = getVariable(0, iY);
        }
        return lb;
    }

    public double getRangeUpperBound(boolean includeInterval) {
        double ub = 1;
        if (tov != null && tov.getNPoints() != 0) {
            ub = getVariable(tov.getNPoints()-1, iY);
        }
        return ub;
    }

    public Range getRangeBounds(boolean includeInterval) {
        return new Range(getRangeLowerBound(includeInterval), getRangeUpperBound(includeInterval));
    }

    /**
     * Returns the specified data element.
     * @param iPoint the index of the point to get the variable of
     * @param iY the variable to get, as specified by the variableNames field
     * @return the variable at the specified point
     */
    private double getVariable(int iPoint, int iY) {
        double variable;
        if (iY == 0) {
            variable = tov.getRadius(iPoint);
        } else if (iY < 4) {
            variable = tov.getVariables(iPoint)[iY - 1];
        } else {
            variable = tov.getSecondaries(iPoint)[iY - 4];
        }
        return variable;
    }
}
