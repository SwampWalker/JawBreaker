/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.math.numerical;

import java.util.ArrayList;
import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class RK4Test extends TestCase implements QuasiLinearODESystem {

    public RK4Test(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of step method, of class RK4.
     */
    public void testStep() {
        System.out.println("step");
        double[] y0 = {0, 0};
        double t = 0.0;
        QuasiLinearODESystem ode = this;
        double h = 0.1;
        double[] y1 = {0.5382550, 0.3196263};
        double[] w1 = RK4.step(y0, t, ode, h);

        // Test
        boolean withinTolerance = true;
        for (int i = 0; i < 2; i++) {
            if (Math.abs(y1[i] - w1[i]) > 1.0e-6) {
                withinTolerance = false;
            }
        }

        if (!withinTolerance) {
            fail("At least one of the evolved steps did not match the solution to within tolerance");
        }
    }

    /**
     * From Example 1 of Section 5.9 of Burden and Faires (pg 316-318)
     *
     * @param t time coordinate
     * @param I the value of the current at time t
     * @return the negative derivative of current, the right hand side
     */
    public double[] rightHandSide(double t, double[] I) {
        double[] f = new double[2];
        f[0] = -4 * I[0] + 3 * I[1] + 6;
        f[1] = -2.4 * I[0] + 1.6 * I[1] + 3.6;
        return f;
    }

    /**
     * Test of evolve method, of class RK4.
     */
    public void testEvolve() {
        System.out.println("evolve");
        ArrayList<double[]> y = new ArrayList<double[]>();
        y.add(new double[]{0, 0});
        ArrayList<Double> t = new ArrayList<Double>();
        t.add(0.);
        QuasiLinearODESystem ode = this;
        double h = 0.1;
        int outputEvery = 1;
        double tMax = 0.5;
        int maxSteps = 5;
        RK4.evolve(y, t, ode, h, outputEvery, tMax, maxSteps);

        if (t.size() != y.size() && y.size() != 6) {
            fail("Incorrect number of steps taken.");
        }

        double[][] yBF = {{0, 0},
            {0.5382550, 0.3196263},
            {0.9684983, 0.5687817},
            {1.310717, 0.7607328},
            {1.581263, 0.9063208},
            {1.793505, 1.014402}
        };
        boolean failed = false;
        for (int i = 0; i < 6 && !failed; i++) {
            double[] yi = y.get(i);
            if (Math.abs(yi[0] - yBF[i][0]) > 1.0E-5 || Math.abs(yi[1] - yBF[i][1]) > 1.0e-5) {
                System.out.println("t=" + t.get(i) + " : " + yi[0] + " " + yi[1]);
                failed = true;
            }
        }
        if (failed) {
            fail("The evolution failed.");
        }
    }
}
