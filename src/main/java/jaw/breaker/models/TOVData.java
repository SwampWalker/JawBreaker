package jaw.breaker.models;

import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class TOVData {
    private ArrayList<double[]> variables;
    private ArrayList<Double> radii;
    private String name = "Unset name";

    public ArrayList<double[]> getVariables() {
        return variables;
    }

    public void setVariables(ArrayList<double[]> variables) {
        this.variables = variables;
    }

    public ArrayList<Double> getRadii() {
        return radii;
    }

    public void setRadii(ArrayList<Double> radii) {
        this.radii = radii;
    }
    
    public TOVData() {
        
    }
    
    public double getRadius(int i) {
        return radii.get(i);
    }
    
    public double getMass(int i) {
        return variables.get(i)[1];
    }
    
    public double getPressure (int i) {
        return variables.get(i)[0];
    }
    
    public double getLambda(int i) {
        return variables.get(i)[2];
    }
    
    public double getTotalMass() {
        return variables.get(variables.size()-1)[1];
    }

    public String getIdentifier() {
        return name;
    }
    
    public int getNPoints() {
        return radii.size();
    }

    public double[] getVariables(int i) {
        return variables.get(i);
    }
}
