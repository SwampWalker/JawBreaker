package ca.tonita.physics.gr.hydro;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.math.numerical.RK4;
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
    private TabulatedHermite eos;

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

    public TOVData(TabulatedHermite eos) {
        variables = new ArrayList<double[]>();
        radii = new ArrayList<Double>();
        this.eos = eos;
    }

    public TabulatedHermite getEos() {
        return eos;
    }

    public void setEos(TabulatedHermite eos) {
        this.eos = eos;
    }

    public double getRadius(int i) {
        return radii.get(i);
    }

    /**
     * Returns an areal radius such that the rest mass inside of that radius is
     * the prescribed value to the tolerance supplied. If the rest mass is
     * greater than the total rest mass of the star, the radius is returned.
     *
     * @param restMass The rest mass to get.
     * @param tolerance The tolerance to use to compute the radius.
     * @return the radius enclosing the supplied rest mass
     */
    public double getRadiusByRestMass(double restMass, double tolerance) {
        if (restMass > getRestMass()) {
            return getRadius();
        }
        int iLower = 0;
        while (variables.get(iLower + 1)[TOVIndex.RESTMASS] < restMass) {
            iLower++;
        }
        TOVEquations eqns = new TOVEquations(eos);
        double r = radii.get(iLower);
        double[] currentVariables = variables.get(iLower);
        // TODO: parameterise this numeric parameter.
        int iStep = 0;
        while (Math.abs(restMass - currentVariables[TOVIndex.RESTMASS]) > tolerance && iStep < 30) {
            double dm0dr = eqns.dm0dr(r, currentVariables);
            double h = (restMass - currentVariables[TOVIndex.RESTMASS]) / dm0dr;
            currentVariables = RK4.step(currentVariables, r, eqns, h);
            r += h;
            iStep++;
        }
        return r;
    }

    public double getMass(int i) {
        return variables.get(i)[TOVIndex.MASS];
    }

    /**
     * Computes the radius of the sphere containing the prescribed rest mass
     *
     * @param restMass The rest mass to be enclosed by the sphere of a given
     * radius.
     * @return the radius enclosing the rest mass prescribed.
     */
    public double enclosingRadius(double restMass) {
        int iLower = 0;
        if (restMass > variables.get(variables.size() - 1)[TOVIndex.RESTMASS]) {
            throw new IllegalArgumentException("Cannot find radius enclosing mass greater than that of the TOV.");
        }
        while (variables.get(iLower + 1)[TOVIndex.RESTMASS] > restMass) {
            iLower++;
        }
        TOVEquations eqns = new TOVEquations(eos);
        double r = radii.get(iLower);
        double dm0dr = eqns.dm0dr(r, variables.get(iLower));
        double h;
        if (iLower != 0) {
            h = (restMass - variables.get(iLower)[TOVIndex.RESTMASS]) / dm0dr;
        } else {
            h = 0.001 * radii.get(0);
            if (h > 0.0001) {
                h = 0.0001;
            }
        }
        double[] enclosingVariables = RK4.step(variables.get(iLower), radii.get(iLower), eqns, h);
        int i = 0;
        while (Math.abs(enclosingVariables[TOVIndex.RESTMASS] - restMass) > 1.0E-12 && i < 10) {
            i++;
            dm0dr = eqns.dpdr(r, enclosingVariables);
            h = (restMass - variables.get(iLower)[TOVIndex.RESTMASS]) / dm0dr;
            enclosingVariables = RK4.step(enclosingVariables, r, eqns, h);
            r += h;
        }
        return r;
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

    /**
     * Returns the value of the TOV variables at the prescribed radius.
     *
     * @param r The radius to get the variables at.
     * @return the variables at that radius.
     */
    public double[] getVariables(double r) {
        int iLower = 0;
        while (radii.get(iLower + 1) < r) {
            iLower++;
        }
        TOVEquations eqns = new TOVEquations(eos);
        return RK4.step(variables.get(iLower), radii.get(iLower), eqns, r - radii.get(iLower));
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
            computeSecondaries(eos);
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

    public double getRestMass() {
        return variables.get(variables.size() - 1)[TOVIndex.RESTMASS];
    }

    /**
     * Returns the radius of the core of the Neutron star.
     *
     * @return the core radius
     */
    public double getCoreRadius() {
        return coreRadius;
    }

    /**
     * Gets the rest mass of the core.
     *
     * @return the core's rest mass
     */
    public double getCoreRestMass() {
        return coreVariables[TOVIndex.RESTMASS];
    }

    /**
     * Returns the value of the mass potential at the surface of the core (the
     * core-crust boundary).
     *
     * @return the mass potential of the core surface
     */
    public double getCoreMass() {
        return coreVariables[TOVIndex.MASS];
    }

    /**
     * Gets all the core variables at the core radius.
     *
     * @return the core variables
     */
    public double[] getCoreVariables() {
        return coreVariables;
    }
}
