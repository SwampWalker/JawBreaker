package ca.tonita.physics.gr.hydro;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.math.numerical.RK4;
import ca.tonita.physics.gr.hydro.TOVEquations;
import ca.tonita.physics.gr.hydro.TOVIndex;
import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class TOVData {

    private ArrayList<double[]> variables;
    private ArrayList<double[]> secondaries;
    private double coreRadius;
    private double[] coreVariables;
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
        return variables.get(i)[TOVIndex.MASS];
    }

    public double getPressure(int i) {
        return variables.get(i)[TOVIndex.PRESSURE];
    }

    public double getLambda(int i) {
        return variables.get(i)[TOVIndex.LAMBDA];
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
        if (secondaries == null) {
            secondaries = new ArrayList<double[]>();
            for (double[] y : variables) {
                double[] z = new double[3];
                z[0] = eos.numberDensity(y[0]);
                z[1] = eos.energyDensity(y[0]);
                z[2] = z[0] * eos.getParticleMass();
                secondaries.add(z);
            }

            // Compute the core quantities.
            int iLower = 0;
            double coreP = eos.getEdgePressure();
            if (variables.get(0)[0] < coreP) {
                coreRadius = 0;
            } else {
                while (variables.get(iLower + 1)[0] > coreP) {
                    iLower++;
                }
                TOVEquations eqns = new TOVEquations(eos);
                double dpdr = eqns.dpdr(radii.get(iLower), variables.get(iLower));
                double h;
                if (iLower != 0) {
                    h = (coreP - variables.get(iLower)[0]) / dpdr;
                } else {
                    h = 0.0001;
                }
                coreVariables = RK4.step(variables.get(iLower), radii.get(iLower), eqns, h);
                coreRadius = radii.get(iLower) + h;
                int i = 0;
                while (Math.abs(coreVariables[0] - coreP) > 1.0E-12 && i < 10) {
                    i++;
                    dpdr = eqns.dpdr(coreRadius, coreVariables);
                    h = (coreP - coreVariables[0]) / dpdr;
                    coreVariables = RK4.step(coreVariables, coreRadius, eqns, h);
                    coreRadius += h;
                }
            }
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
        return variables.get(variables.size()-1)[TOVIndex.RESTMASS];
    }

    public double getCoreRadius() {
        return coreRadius;
    }

    public double getCoreRestMass() {
        return coreVariables[TOVIndex.RESTMASS];
    }

    public double getCoreMass() {
        return coreVariables[TOVIndex.MASS];
    }
}
