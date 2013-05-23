package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.TOVData;
import ca.tonita.math.numerical.QuasiLinearFirstOrderODESystem;
import ca.tonita.math.numerical.RK4;
import ca.tonita.physics.gr.TOVEquations;

/**
 * Implements the equations of the body manifold given an RK4 solution as the
 * background.
 *
 * @author atonita
 */
public class SphericalBodyManifoldRK4 {

    private final TOVData solution;
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
        this.solution = solution;
        this.eos = eos;
        radius = solution.getRadius();
        inverseStepSize = 1./solution.getRadius(1);
        eqns = new TOVEquations(eos);
    }

    /**
     * Returns all the quantities of interest at the given radius.
     *
     * @param r The radius to compute elastic terms at.
     * @return a <code>SphericalElasticBean</code> containing the elastic values
     */
    public SphericalElasticBean getQuantities(double r) {
        if (r > radius) {
            throw new UnsupportedOperationException("Radius out of bounds: maximum radius is " + radius + " received " + r);
        } else if (r < 0) {
            throw new UnsupportedOperationException("Radius out of bounds: must be greater than zero, received " + r);
        }
        SphericalElasticBean bean = new SphericalElasticBean();
        int iLowerBound = (int)(r * inverseStepSize);
        if (solution.getRadius(iLowerBound) > r) {
            iLowerBound--;
        }
        double h = r - solution.getRadius(iLowerBound);
        double[] newPoint = RK4.step(solution.getVariables(iLowerBound), solution.getRadius(iLowerBound), eqns, h);
        bean.setPressure(newPoint[0]);
        bean.setDpressure(eqns.dmdr(r, newPoint));
        bean.setMassPotential(newPoint[1]);
        bean.setdMassPotential(eqns.dmdr(r, newPoint));
        
        bean.setNumberDensity(eos.numberDensity(bean.getPressure()));
        bean.setEnergyPerParticle(eos.energyPerParticle(bean.getPressure()));
        bean.setDnumberDensity(eos.dnumberDensity(bean.getPressure()));
        bean.setDenergyPerParticle(eos.denergyPerParticle(bean.getPressure()));
        return bean;
    }
}
