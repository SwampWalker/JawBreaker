package ca.tonita.math.numerical;

/**
 *
 * @author atonita
 */
public class EvolutionTerminator {
    
    /**
     * The maximum value of the independent variable to allow.
     */
    protected double tMax;
    
    /**
     * The maximum number of steps to allow.
     */
    protected int maxSteps;

    public EvolutionTerminator(double tMax, int maxSteps) {
        this.tMax = tMax;
        this.maxSteps = maxSteps;
    }

    public double getTMax() {
        return tMax;
    }

    public void setTMax(double tMax) {
        this.tMax = tMax;
    }

    public int getMaxSteps() {
        return maxSteps;
    }

    public void setMaxSteps(int maxSteps) {
        this.maxSteps = maxSteps;
    }

    /**
     * Determines whether to terminate evolution.
     *
     * @param y the variables being evolved
     * @param t the independent variable
     * @param steps the number of steps in the evolution so far.
     * @return whether to terminate the evolution or not
     */
    public boolean terminate(double[] y, double t, int steps) {
        if (steps >= maxSteps) {
            return true;
        } else if (t >= tMax) {
            return true;
        }
        return false;
    }
}
