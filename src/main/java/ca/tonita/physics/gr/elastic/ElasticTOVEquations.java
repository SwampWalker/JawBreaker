package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;

/**
 *
 * @author atonita
 */
public final class ElasticTOVEquations {

    private final SphericalBodyManifoldRK4 body;

    /**
     * Constructs the ElasticTOVEquations. It needs a body manifold and an
     * equation of state.
     *
     * @param body
     * @param eos
     */
    public ElasticTOVEquations(SphericalBodyManifoldRK4 body, TabulatedHermite eos) {
        this.body = body;
    }

    /**
     * Returns the trace of the strain tensor.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the normalized trace of the strain
     */
    public double trace(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double GrrOvergrr = (1 - 2 * m / r) / (1 - 2 * bodyVars.getMassPotential() / xi);
        return GrrOvergrr * xiPrime * xiPrime + 2 * xi * xi / (r * r) - 3;
    }

    /**
     * Returns the derivative of the trace of the strain tensor with respect to
     * the configuration.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the normalized trace of the strain
     */
    protected double dtracedxi(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        return (1 - 2 * m / r) / Math.pow(1 - 2 * bodyVars.getMassPotential() / xi, 2)
                * (2 / xi * bodyVars.getdMassPotential() - 2 * bodyVars.getMassPotential() / (xi * xi))
                + 4 * xi / (r * r);
    }

    /**
     * Returns the derivative of the trace of the strain tensor with respect to
     * the configuration gradient.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the normalized trace of the strain
     */
    protected double dtracedxiPrime(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double GrrOvergrr = (1 - 2 * m / r) / (1 - 2 * bodyVars.getMassPotential() / xi);
        return 2 * GrrOvergrr * xiPrime;
    }

    /**
     * Returns the derivative of the trace of the strain tensor with respect to
     * the mass potential.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the normalized trace of the strain
     */
    protected double dtracedmass(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        return - 2 * xiPrime * xiPrime / r / (1 - 2 * bodyVars.getMassPotential() / xi);
    }

    /**
     * Returns the isotropic pressure component.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the isotropic component of the pressure.
     */
    protected double pressureIsotropic(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        return bodyVars.getPressure() + bodyVars.getLameLambda() * 0.25 * trace(r, bodyVars, xi, xiPrime, m);
    }

    /**
     * Returns the volume contraction factor: the ratio of the new number
     * density to the old number density.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the volume contraction factor.
     */
    protected double volumeContractionFactor(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double xiOverR = xi / r;
        double M = bodyVars.getMassPotential();
        return xiOverR * xiOverR * xiPrime * Math.sqrt((1. - 2 * m / r) / (1. - 2 * M / xi));
    }

    /**
     * Returns the tangential pressure. This is the (\theta,\theta) or
     * (\phi,\phi) component of the mixed pressure tensor.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the pressure, tangential
     */
    public double pressureTangential(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double chi = volumeContractionFactor(r, bodyVars, xi, xiPrime, m);
        double P_I = pressureIsotropic(r, bodyVars, xi, xiPrime, m);
        double mu = bodyVars.getShearModulus();
        double xiOverR = xi / r;
        return chi * (P_I + mu * 0.5 * (xiOverR * xiOverR * xiPrime * xiPrime - 1.)) * xiOverR * xiOverR * xiPrime * xiPrime;
    }

    /**
     * Returns the radial pressure. This is the (r,r) component of the mixed
     * pressure tensor.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the pressure, radial
     */
    public double pressureRadial(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double chi = volumeContractionFactor(r, bodyVars, xi, xiPrime, m);
        double P_I = pressureIsotropic(r, bodyVars, xi, xiPrime, m);
        double mu = bodyVars.getShearModulus();
        double M = bodyVars.getMassPotential();
        double ratio = (1. - 2 * m / r) / (1. - 2 * M / xi) * xiPrime * xiPrime;
        return chi * (P_I + mu * 0.5 * (ratio - 1.)) * ratio;
    }
}
