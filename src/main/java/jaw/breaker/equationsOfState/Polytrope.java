/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.equationsOfState;

/**
 *
 * @author tonita
 */
public class Polytrope {
    private double kappa;
    private double gamma;
    private double particleMass;
    private double inverseGammaMinusOne;

    /**
     * Constructs a polytropic equation of state. Here, c the speed of light and
     * G the gravitational constant are taken equal to 1 and unitless. That is
     * <code>c = G = 1</code>. The relation between number density and pressure
     * is <code>P = kappa*(particleMass*n)^gamma</code>.
     * 
     * @param kappa The polytropic constant.
     * @param gamma The adiabatic index.
     * @param particleMass The mass of a particle, whatever the particle to be
     * considered is.
     */
    public Polytrope(double kappa, double gamma, double particleMass) {
        this.kappa = kappa;
        this.gamma = gamma;
        this.particleMass = particleMass;
        inverseGammaMinusOne = 1./(gamma - 1);
    }
    
    /**
     * The relation is <code>P = kappa*(particleMass*n)^gamma</code>.
     * @param numberDensity The number density of particles.
     * @return The pressure.
     */
    public double pressure(double numberDensity) {
        return kappa*Math.pow(particleMass*numberDensity, gamma);
    }
    
    /**
     * Returns the total energy density as a function of number density.
     * @param numberDensity The number density.
     * @return the total energy density (rest mass energy plus internal energy)
     */
    public double energyDensity(double numberDensity) {
        // P = rho*epsilon*(gamma-1)
        double epsilon = kappa*Math.pow(particleMass*numberDensity, gamma-1)*inverseGammaMinusOne;
        return particleMass*numberDensity*(1+epsilon);
    }
    
    /**
     * Essentially the internal energy, divided by the number of particles. Also
     * equal to the specific internal energy multiplied by the particle mass in
     * units where <code>c=G=1</code>.
     * 
     * @param numberDensity the number density of particles
     * @return the energy per particle excluding rest mass
     */
    public double energyPerParticle(double numberDensity) {
        double epsilon = kappa*Math.pow(particleMass*numberDensity, gamma-1)*inverseGammaMinusOne;
        return epsilon * particleMass;
    }
}
