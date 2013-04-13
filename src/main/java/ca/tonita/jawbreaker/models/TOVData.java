package ca.tonita.jawbreaker.models;

import java.util.ArrayList;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;

/**
 *
 * @author atonita
 */
public class TOVData {

    private ArrayList<double[]> variables;
    private ArrayList<double[]> secondaries;
    private ArrayList<Double> radii;
    private String name = "Unset name";

    public ArrayList<double[]> getVariables() {
        return variables;
    }

    public void setVariables(ArrayList<double[]> variables) {
        this.variables = variables;
        this.secondaries = null;
    }

    public ArrayList<Double> getRadii() {
        return radii;
    }

    public void setRadii(ArrayList<Double> radii) {
        this.radii = radii;
    }

    public TOVData() {
        variables = new ArrayList<double[]>();
        radii = new ArrayList<Double>();
    }

    public double getRadius(int i) {
        return radii.get(i);
    }

    public double getMass(int i) {
        return variables.get(i)[1];
    }

    public double getPressure(int i) {
        return variables.get(i)[0];
    }

    public double getLambda(int i) {
        return variables.get(i)[2];
    }

    public double getTotalMass() {
        return variables.get(variables.size() - 1)[1];
    }

    public String getIdentifier() {
        return name;
    }

    public int getNPoints() {
        if (radii == null) {
            return 0;
        }
        return radii.size();
    }

    public double[] getVariables(int i) {
        return variables.get(i);
    }
    
    public void computeSecondaries(TabulatedHermite eos) {
        secondaries = new ArrayList<double[]>();
        for (double[] y : variables) {
            double[] z = new double[2];
            z[0] = eos.numberDensity(y[0]);
            z[1] = eos.energyDensity(y[0]);
            secondaries.add(z);
        }
    }
    
    public double[] getSecondaries(int i) {
        if (secondaries == null) {
            throw new UnsupportedOperationException("Secondaries not computed.");
        }
        return secondaries.get(i);
    }
}