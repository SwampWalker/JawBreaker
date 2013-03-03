/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.math.numerical;

import java.util.ArrayList;

/**
 *
 * @author atonita
 */
public class RK4 {

    /**
     * To avoid constant unnecessary divisions.
     */
    private static final double sixth = 1. / 6.;

    /**
     * Takes a single step using RK4.
     *
     * @param y0 The initial condition.
     * @param t The current value of the independent variable.
     * @param ode The <code>QuasiLinearODESystem</code> being approximately
     * solved.
     * @param h The step size to take.
     * @return The approximate value of y at t + h.
     */
    public static double[] step(double[] y0, double t, QuasiLinearODESystem ode, double h) {
        int n = y0.length;
        double[] k1 = ode.rightHandSide(t, y0);
        double[] yNext = new double[n];
        for (int i = 0; i < n; i++) {
            yNext[i] = y0[i] + 0.5 * h * k1[i];
        }
        double[] k2 = ode.rightHandSide(t + 0.5 * h, yNext);
        for (int i = 0; i < n; i++) {
            yNext[i] = y0[i] + 0.5 * h * k2[i];
        }
        double[] k3 = ode.rightHandSide(t + 0.5 * h, yNext);
        for (int i = 0; i < n; i++) {
            yNext[i] = y0[i] + h * k3[i];
        }
        double[] k4 = ode.rightHandSide(t + h, yNext);
        for (int i = 0; i < n; i++) {
            yNext[i] = y0[i] + h * sixth * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }

        return yNext;
    }

    /**
     * Evolves the ODE system. New values are added to the
     * <code>ArrayList</code> y.
     *
     * @param y The values. Should include y0 as the last element.
     * @param t The values of t for the corresponding y. Should contain t0.
     * @param ode The ode system to evolve.
     * @param h The step size to use.
     * @param outputEvery The number of steps to evolve before outputting. A
     * value of 1 meaning output every step.
     * @param tMax The maximum time to evolve to.
     * @param maxSteps The maximum number of steps to take from the last value
     * of t in yOfT.
     */
    public static void evolve(ArrayList<double[]> y, ArrayList<Double> t, QuasiLinearODESystem ode, double h, int outputEvery, double tMax, int maxSteps) {
        double t0 = t.get(y.size() - 1);
        double[] y0 = new double[y.get(0).length];
        System.arraycopy(y.get(y.size() - 1), 0, y0, 0, y0.length);
        for (int i = 0; i < maxSteps && t0 + h <= tMax; i++) {
            y0 = step(y0, t0, ode, h);
            t0 = t0 + h;
            if (outputEvery > 0 && i % outputEvery == 0) {
                y.add(y0);
                t.add(t0);
            }
        }
    }
}
