package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.gr.hydro.TOVData;
import ca.tonita.math.numerical.QuasiLinearFirstOrderODESystem;
import ca.tonita.math.numerical.RK4;
import ca.tonita.physics.gr.hydro.TOVEquations;

/**
 * Implements the equations of the body manifold given an RK4 solution as the
 * background.
 *
 * @author atonita
 */
public class SphericalBodyManifoldRK4 implements QuasiLinearFirstOrderODESystem {

    private final TOVData background;
    private final double radius;
    private final double inverseStepSize;
    private final TabulatedHermite eos;
    private final TOVEquations eqns;

    /**
     * Creates the body manifold.
     *
     * @param solution The TOV solution to use as background. The points must be
     * equally spaced.
     */
    public SphericalBodyManifoldRK4(TOVData solution, TabulatedHermite eos) {
        this.background = solution;
        this.eos = eos;
        radius = solution.getRadius();
        inverseStepSize = 1./solution.getRadius(1);
        eqns = new TOVEquations(eos);
    }

    /**
     * Returns all the quantities of interest at the given radius.
     *
     * @param xi The radius to compute elastic terms at.
     * @return a <code>SphericalElasticBean</code> containing the elastic values
     */
    public SphericalElasticBean getQuantities(double xi) {
        if (xi > radius) {
            throw new UnsupportedOperationException("Radius out of bounds: maximum radius is " + radius + " received " + xi);
        } else if (xi < 0) {
            throw new UnsupportedOperationException("Radius out of bounds: must be greater than zero, received " + xi);
        }
        SphericalElasticBean bean = new SphericalElasticBean();
        int iLowerBound = (int)(xi * inverseStepSize);
        if (background.getRadius(iLowerBound) > xi) {
            iLowerBound--;
        }
        double h = xi - background.getRadius(iLowerBound);
        double[] newPoint = RK4.step(background.getVariables(iLowerBound), background.getRadius(iLowerBound), eqns, h);
        bean.setPressure(newPoint[0]);
        bean.setDpressure(eqns.dmdr(xi, newPoint));
        bean.setMassPotential(newPoint[1]);
        bean.setdMassPotential(eqns.dmdr(xi, newPoint));
        
        bean.setNumberDensity(eos.numberDensity(bean.getPressure()));
        bean.setEnergyPerParticle(eos.energyPerParticle(bean.getPressure()));
        bean.setDnumberDensity(eos.dnumberDensity(bean.getPressure()));
        bean.setDenergyPerParticle(eos.denergyPerParticle(bean.getPressure()));
        
        bean.setLameLambda(eos.lambda(bean.getPressure()));
        bean.setShearModulus(eos.shearModulus(bean.getPressure()));
        return bean;
    }

    /**
     * Returns the particle mass of the average particle in the body manifold.
     * @return the particle mass
     */
    public double particleMass() {
        return eos.getParticleMass();
    }
    
    public TOVData getBackground() {
        return background;
    }

    public double[] rightHandSide(double t, double[] y) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
