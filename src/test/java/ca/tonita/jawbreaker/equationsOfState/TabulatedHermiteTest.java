package ca.tonita.jawbreaker.equationsOfState;

import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class TabulatedHermiteTest extends TestCase {

    TabulatedHermite instance;
    private Polytrope poly;
    private double logNMin = -5;
    private double logNMax = -1;

    public TabulatedHermiteTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        double kappa = 100;
        double gamma = 2;
        double particleMass = 1;
        int nPoints = 10000;
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
        instance = new TabulatedHermite(logn, logp, energyPerParticle, particleMass);
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of getIdentifier method, of class TabulatedHermite.
     */
    public void testGetIdentifier() {
        System.out.println("getIdentifier");
        String expResult = "fixed";
        instance.setIdentifier(expResult);
        String result = instance.getIdentifier();
        assertEquals(expResult, result);
    }

    /**
     * Test of setIdentifier method, of class TabulatedHermite.
     */
    public void testSetIdentifier() {
        System.out.println("setIdentifier");
        String identifier = "";
        instance.setIdentifier(identifier);
    }

    /**
     * Test of toString method, of class TabulatedHermite.
     */
    public void testToString() {
        System.out.println("toString");
        String expResult = "expResult";
        instance.setIdentifier(expResult);
        String result = instance.toString();
        assertEquals(expResult, result);
    }

    /**
     * Test of energyDensity method, of class TabulatedHermite.
     */
    public void testEnergyDensity() {
        System.out.println("energyDensity");
        double[] n = getN();
        for (int i = 0; i < 3; i++) {
            double expResult = poly.energyDensity(n[i]);
            double result = instance.energyDensity(poly.pressure(n[i]));
            assertEquals(expResult, result, 1.0e-2 * expResult);
        }
    }

    /**
     * Test of denergyDensity method, of class TabulatedHermite.
     */
    public void testDenergyDensity() {
        System.out.println("denergyDensity");
        double[] p = getP();
        double h = 1.0e-6;
        for (int i = 0; i < 3; i++) {
            double ep = instance.energyDensity(p[i] * (1 + h));
            double em = instance.energyDensity(p[i] * (1 - h));
            double de = instance.denergyDensity(p[i]);
            assertEquals((ep - em) * 0.5 / p[i] / h, de, de * h);
        }
    }

    /**
     * Test of cloneTable method, of class TabulatedHermite.
     */
    public void testCloneTable() {
        System.out.println("cloneTable");
        double[][] table = new double[4][];
        instance.cloneTable(table);
    }

    /**
     * Test of numberDensity method, of class TabulatedHermite.
     */
    public void testNumberDensity() {
        System.out.println("numberDensity");
        double[] n = getN();
        for (int i = 0; i < 3; i++) {
            double result = instance.numberDensity(poly.pressure(n[i]));
            assertEquals(n[i], result, 1.0e-3 * n[i]);
        }
    }

    /**
     * Test of dnumberDensity method, of class TabulatedHermite.
     */
    public void testDnumberDensity() {
        System.out.println("dnumberDensity");
        double[] p = getP();
        double h = 1.0e-6;
        for (int i = 0; i < 3; i++) {
            double np = instance.numberDensity(p[i] * (1 + h));
            double nm = instance.numberDensity(p[i] * (1 - h));
            double dn = instance.dnumberDensity(p[i]);
            assertEquals((np - nm) * 0.5 / p[i] / h, dn, dn * h);
        }
    }

    /**
     * Test of getParticleMass method, of class TabulatedHermite.
     */
    public void testGetParticleMass() {
        System.out.println("getParticleMass");
        double expResult = 1.0;
        double result = instance.getParticleMass();
        assertEquals(expResult, result, 1.0e-6);
    }

    /**
     * Test of energyPerParticle method, of class TabulatedHermite.
     */
    public void testEnergyPerParticle() {
        System.out.println("energyPerParticle");
        double[] p = getP();
        for (int i = 0; i < 3; i++) {
            double expResult = instance.energyDensity(p[i]);
            double result = (instance.getParticleMass() + instance.energyPerParticle(p[i])) * instance.numberDensity(p[i]);
            assertEquals(expResult, result, 1.0e-4);
        }
    }

    /**
     * Test of denergyPerParticle method, of class TabulatedHermite.
     */
    public void testDEnergyPerParticle() {
        System.out.println("dEnergyPerParticle");
        double[] p = getP();
        double h = 1.0e-6;
        for (int i = 0; i < 3; i++) {
            double ep = instance.energyPerParticle(p[i]*(1 + h));
            double em = instance.energyPerParticle(p[i]*(1 - h));
            double expResult = (ep - em)*0.5/p[i]/h;
            double result = instance.denergyPerParticle(p[i]);
            assertEquals(expResult, result, result*h);
        }
    }

    /**
     * Test of getEdgePressure method, of class TabulatedHermite.
     */
    public void testGetEdgePressure() {
        System.out.println("getEdgePressure");
        double result = instance.getEdgePressure();
        double numberDensity = instance.numberDensity(result);
        assertEquals(2.1707e-4, numberDensity, 1.0e-8);
    }

    /**
     * Test of getEdgeDensity method, of class TabulatedHermite.
     */
    public void testGetEdgeDensity() {
        System.out.println("getEdgeDensity");
        double result = instance.getEdgeDensity();
        assertEquals(2.1707e-4, result, 1.0e-8);
    }

    private double[] getN() {
        double min = Math.pow(10, logNMin);
        double max = Math.pow(10, logNMax);
        return new double[]{
                    0.5 * min,
                    min + (max - min) * Math.random(),
                    2 * max
                };
    }

    private double[] getP() {
        double[] n = getN();
        double[] p = new double[3];
        for (int i = 0; i < 3; i++) {
            p[i] = poly.pressure(n[i]);
        }
        return p;
    }
}
