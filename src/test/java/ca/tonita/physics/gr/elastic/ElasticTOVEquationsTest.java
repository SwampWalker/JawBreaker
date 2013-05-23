package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.Polytrope;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.TOVData;
import ca.tonita.physics.gr.TOVBuilder;
import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class ElasticTOVEquationsTest extends TestCase {

    // Equation of state.
    private double kappa = 100;
    private double gamma = 2;
    private double particleMass = 1;
    private Polytrope poly;
    // Tabulated equation of state.
    private int nPoints = 10000;
    private double logNMin = -10;
    private double logNMax = -2;
    private TabulatedHermite eos;
    // Background model.
    private double centralPressure = 1.0e-4;
    private double stepSize = 1.0e-3;
    private double terminationPressure = 1.0e-10;
    private TOVData background = new TOVData();
    // Body manifold.
    private SphericalBodyManifoldRK4 body;
    private ElasticTOVEquations instance;

    public ElasticTOVEquationsTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();

        poly = new Polytrope(kappa, gamma, particleMass);
        double[] logn = new double[nPoints];
        double[] logp = new double[nPoints];
        double[] energyPerParticle = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            logn[i] = logNMin + i * (logNMax - logNMin) / (nPoints - 1);
            double n = Math.pow(10, logn[i]);
            logp[i] = Math.log10(poly.pressure(n));
            energyPerParticle[i] = poly.energyPerParticle(n);
        }
        eos = new TabulatedHermite(logn, logp, energyPerParticle, particleMass);
        TOVBuilder.evolve(background, eos, centralPressure, stepSize, 1, terminationPressure);
        body = new SphericalBodyManifoldRK4(background, eos);
        instance = new ElasticTOVEquations(body, eos);
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of dtracedxi method, of class ElasticTOVEquations.
     */
    public void testDtracedxi() {
        System.out.println("dtracedxi");
        double r = background.getRadius()*Math.random();
        SphericalElasticBean bodyVars = body.getQuantities(r);
        double h = 1.0e-6;
        double xi = r;
        double xiPrime = 1;
        double m = bodyVars.getMassPotential();
        SphericalElasticBean bodyVarsp = body.getQuantities(xi + h);
        SphericalElasticBean bodyVarsm = body.getQuantities(xi - h);
        double Hp = instance.trace(r, bodyVarsp, xi + h, xiPrime, m);
        double Hm = instance.trace(r, bodyVarsm, xi - h, xiPrime, m);
        double expResult = (Hp - Hm)/(2*h);
        double result = instance.dtracedxi(r, bodyVars, xi, xiPrime, m);
        assertEquals(expResult, result, 1.0e-6);
    }

    /**
     * Test of dtracedxiPrime method, of class ElasticTOVEquations.
     */
    public void testDtracedxiPrime() {
        System.out.println("dtracedxiPrime");
        double r = background.getRadius()*Math.random();
        SphericalElasticBean bodyVars = body.getQuantities(r);
        double xi = r;
        double xiPrime = 1;
        double m = bodyVars.getMassPotential();
        double h = 1.0e-6;
        double Hp = instance.trace(r, bodyVars, xi, xiPrime + h, m);
        double Hm = instance.trace(r, bodyVars, xi, xiPrime - h, m);
        double expResult = (Hp - Hm)/(2*h);
        double result = instance.dtracedxiPrime(r, bodyVars, xi, xiPrime, m);
        assertEquals(expResult, result, 1.0e-6);
    }

    /**
     * Test of dtracedmass method, of class ElasticTOVEquations.
     */
    public void testTrace() {
        System.out.println("trace");
        double r = background.getRadius()*Math.random();
        SphericalElasticBean bodyVars = body.getQuantities(r);
        double xi = r;
        double xiPrime = 1;
        double m = bodyVars.getMassPotential();
        double expResult = 0.0;
        double result = instance.trace(r, bodyVars, xi, xiPrime, m);
        assertEquals(expResult, result, 1.0e-4);
    }

    /**
     * Test of trace method, of class ElasticTOVEquations.
     */
    public void testDtracedmass() {
        System.out.println("dtracedmass");
        double r = background.getRadius()*Math.random();
        SphericalElasticBean bodyVars = body.getQuantities(r);
        double xi = r;
        double xiPrime = 1;
        double m = bodyVars.getMassPotential();
        double h = 1.0e-6;
        double Hp = instance.trace(r, bodyVars, xi, xiPrime, m + h);
        double Hm = instance.trace(r, bodyVars, xi, xiPrime, m - h);
        double expResult = (Hp - Hm)/(2*h);
        double result = instance.dtracedmass(r, bodyVars, xi, xiPrime, m);
        assertEquals(expResult, result, 1.0e-6);
    }
}
