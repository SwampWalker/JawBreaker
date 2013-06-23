package ca.tonita.jawbreaker.models;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class TOVData {

    private ArrayList<double[]> variables;
    private ArrayList<double[]> secondaries;
    private ArrayList<Double> radii;
    private String name = "Unset name";
    private double conservedMass;

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
    
    public double getGravitationalMass() {
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
            double[] z = new double[3];
            z[0] = eos.numberDensity(y[0]);
            z[1] = eos.energyDensity(y[0]);
            z[2] = z[0]*eos.getParticleMass();
            secondaries.add(z);
        }
    }

    public double[] getSecondaries(int i) {
        if (secondaries == null) {
            throw new UnsupportedOperationException("Secondaries not computed.");
        }
        return secondaries.get(i);
    }

    /**
     * Returns the radius of the star (actually the last radius computed).
     *
     * @return the last radius where the star was computed.
     */
    public double getRadius() {
        return radii.get(radii.size() - 1);
    }

    public double getConservedMass() {
        return conservedMass;
    }

    /**
     * Sets the conserved mass.
     * @param mass the mass to set
     */
    public void setConservedMass(double mass) {
        conservedMass = mass;
    }

    public double getCoreRadius() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public double getCoreRestMass() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
