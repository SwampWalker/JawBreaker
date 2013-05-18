/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.datasets;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class EOSDataset extends MultiDataset<EOSBean> {

    private final static String[] eosVariableNames = new String[]{"Number Density", "Pressure", "Total Energy Density", "Total Energy Density Derivative"};
    private ArrayList<String> names = null;

    public EOSDataset() {
        variableNames = eosVariableNames;
        iX = 1;
        iY = 2;
    }

    /**
     * Adds an equation of state to this data set.
     *
     * @param index the index into the datasets to add the table
     * @param eos the equation of state to add.
     */
    public void add(int index, TabulatedHermite eos) {
        double[][] table = new double[4][];
        eos.cloneTable(table);
        EOSBean bean = new EOSBean();
        bean.setData(table);
        bean.setName(eos.getIdentifier());
        super.add(index, bean);
    }

    @Override
    protected double getItem(int iSeries, int iVariable, int iItem) {
        return getDataset(iSeries).getData()[iVariable][iItem];
    }

    @Override
    public int getItemCount(int series) {
        return getDataset(series).getData()[0].length;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        return getDataset(series).getName();
    }
}
