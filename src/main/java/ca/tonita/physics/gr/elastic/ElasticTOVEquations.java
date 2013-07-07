package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;

/**
 * This particular set of the ElasticEquations assumes a specific equation of
 * state in terms of the right Cauchy-Green tensor. See the notes for the
 * precise formulation of the elastic equation of state assumed.
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
     */
    public ElasticTOVEquations(SphericalBodyManifoldRK4 body) {
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
     * Returns the trace of the square of the stress tensor minus the delta.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the square of the normalized stress
     */
    public double traceOfSquare(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double GrrOvergrr = (1 - 2 * m / r) / (1 - 2 * bodyVars.getMassPotential() / xi);
        return Math.pow(GrrOvergrr * xiPrime * xiPrime - 1, 2) + 2 * Math.pow(xi * xi / (r * r) - 1, 2);
    }
    
    public double energyPerParticle(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double trace = trace(r, bodyVars, xi, xiPrime, m);
        double traceOfSquare = traceOfSquare(r, bodyVars, xi, xiPrime, m);
        return body.particleMass() * (1 + bodyVars.getEnergyPerParticle()) + bodyVars.getPressure() * trace / bodyVars.getNumberDensity()
                + bodyVars.getLameLambda() * 0.125 / bodyVars.getNumberDensity() * trace * trace
                + bodyVars.getShearModulus() * 0.25 / bodyVars.getNumberDensity() * traceOfSquare;
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
    protected double dtrace_dxiPrime(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
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
     * Returns the isotropic pressure component.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @return the isotropic component of the pressure.
     */
    protected double dpressureIsotropic_dxiPrime(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        return 0.25 * bodyVars.getLameLambda() * dtrace_dxiPrime(r, bodyVars, xi, xiPrime, m);
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
    protected double dvolumeContractionFactor_dxiPrime(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double xiOverR = xi / r;
        double M = bodyVars.getMassPotential();
        return xiOverR * xiOverR * Math.sqrt((1. - 2 * m / r) / (1. - 2 * M / xi));
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
    public double dpressureRadial_dxiPrime(double r, SphericalElasticBean bodyVars, double xi, double xiPrime, double m) {
        double chi = volumeContractionFactor(r, bodyVars, xi, xiPrime, m);
        double dchi = dvolumeContractionFactor_dxiPrime(r, bodyVars, xi, xiPrime, m);
        double P_I = pressureIsotropic(r, bodyVars, xi, xiPrime, m);
        double dP_I = dpressureIsotropic_dxiPrime(r, bodyVars, xi, xiPrime, m);
        double mu = bodyVars.getShearModulus();
        double M = bodyVars.getMassPotential();
        double dratioHalf = (1. - 2 * m / r) / (1. - 2 * M / xi) * xiPrime;
        double ratio = dratioHalf * xiPrime;
        double term1 = (2 * chi + xiPrime * dchi) * (P_I + mu * 0.5 * (ratio - 1.)) * ratio;
        double term2 = chi * (dP_I + mu * dratioHalf) * ratio * xiPrime * xiPrime;
        return term1 + term2;
    }

    /**
     * Returns the xiPrime by inverting the constituent equation of the radial
     * pressure. The inversion of the nonlinear constituent equation is
     * performed via using Newton-Raphson iterations.
     *
     * @param r the radius
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @param pRadial The radial pressure.
     * @return the pressure, radial
     */
    public double xiPrime(double r, SphericalElasticBean bodyVars, double xi, double m, double pRadial) {
        double xiPrime = 1;
        double tolerance = 1.0e-15;
        double pressureGuess = pressureRadial(r, bodyVars, xi, xiPrime, m);
        int i = 0;
        while (Math.abs(pressureGuess - pRadial) > tolerance && i < 25) {
            i++;
            xiPrime = xiPrime + (pRadial - pressureGuess) / dpressureRadial_dxiPrime(r, bodyVars, xi, xiPrime, m);
            pressureGuess = this.pressureRadial(r, bodyVars, xi, xiPrime, m);
        }
        return xiPrime;
    }

    /**
     * Returns the derivative of mass with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @param pRadial The radial pressure.
     * @return the derivative of mass wrt r
     */
    public double dmdr(double r, SphericalElasticBean bodyVars, double xi, double m, double pRadial) {
        double rho = volumeContractionFactor(r, bodyVars, xi, xi, m) * bodyVars.getNumberDensity() * energyPerParticle(r, bodyVars, xi, xi, m);
        return 4 * Math.PI * rho * r * r;
    }

    /**
     * Returns the derivative of lambda with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @param pRadial The radial pressure.
     * @return the derivative of lambda wrt r
     */
    public double dlambdadr(double r, SphericalElasticBean bodyVars, double xi, double m, double pRadial) {
        if (r == 0.0) {
            return 0;
        }
        return (m + 4 * Math.PI * Math.pow(r, 3) * pRadial)
                / (r * (r - 2 * m));
    }

    /**
     * Returns the derivative of pressure with respect to r.
     *
     * @param r The areal radial coordinate r.
     * @param bodyVars the quantities from the body manifold
     * @param xi the configuration
     * @param xiPrime the configuration gradient
     * @param m the mass potential
     * @param pRadial The radial pressure.
     * @param xiPrime The configuration gradient.
     * @return the derivative of pressure wrt r
     */
    public double dpdr(double r, SphericalElasticBean bodyVars, double xi, double m, double pRadial, double xiPrime) {
        if (r == 0.0) {
            return 0;
        }
        double pTangential = pressureTangential(r, bodyVars, xi, xiPrime, m);
        double dlambdadr = (m + 4 * Math.PI * Math.pow(r, 3) * pRadial)
                / (r * (r - 2 * m));
        double rho = volumeContractionFactor(r, bodyVars, xi, xi, m) * bodyVars.getNumberDensity() * energyPerParticle(r, bodyVars, xi, xi, m);
        return -(pRadial + rho) * dlambdadr + 2. / r * (pTangential - pRadial);
    }

    /**
     * Returns the value of xi' of the interior surface of the elastic crust for
     * which the tangential pressure equals the radial pressure.
     *
     * @param r The areal radius of the core's inner surface in the spacetime
     * manifold.
     * @param m The mass potential of the space time at the core's inner
     * surface.
     * @return The necessary value of xi' to ensure isotropy.
     */
    public double xiPrimeCoreIsotropic(double r, double m) {
        double R = body.getBackground().getCoreRadius();
        double M = body.getBackground().getCoreMass();
        return Math.sqrt((R * (R - 2 * M)) / (r * (r - 2 * m)));
    }
}
