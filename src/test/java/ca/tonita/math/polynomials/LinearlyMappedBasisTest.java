package ca.tonita.math.polynomials;

import ca.tonita.math.linearalgebra.LinearAlgebra;
import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class LinearlyMappedBasisTest extends TestCase {
    private double tolerance = 1.0e-10;

    public LinearlyMappedBasisTest(String testName) {
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
     * Test of getAbscissas method, of class LinearlyMappedBasis.
     */
    public void testGetAbscissas() {
        System.out.println("getAbscissas");
        int rank = 6;
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setRank(rank);
        double[] expResult = new double[6];
        for (int i = 0; i < rank; i++) {
            expResult[i] = 0.5 * (Math.cos(-Math.PI + Math.PI * i / (6 - 1)) + 1);
        }
        double[] result = instance.getAbscissas();
        for (int i = 0; i < rank; i++) {
            if (Math.abs(expResult[i] - result[i]) > tolerance) {
                fail("At point " + i + " expected " + expResult[i] + " received " + result[i] + ".");
            }
        }
    }

    /**
     * Test of getDifferentiationMatrix method, of class LinearlyMappedBasis.
     */
    public void testGetDifferentiationMatrix() {
        System.out.println("getDifferentiationMatrix");
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        int rank = 20;
        instance.setRank(rank);
        double[] x = instance.getAbscissas();
        double[] u = new double[rank];
        // Derivative of e^x is e^x.
        for (int i = 0; i < rank; i++) {
            u[i] = Math.exp(x[i]);
        }
        double[][] result = instance.getDifferentiationMatrix();
        double[] du = LinearAlgebra.matrixVectorMultiply(result, u);
        for (int i = 0; i < rank; i++) {
            if (Math.abs(du[i] - u[i]) > tolerance) {
                fail("Failure in derivative at " + i + " expected " + u[i] + " computed " + du[i] + ".");
            }
        }
    }

    /**
     * Test of getRank method, of class LinearlyMappedBasis.
     */
    public void testGetRank() {
        System.out.println("getRank");
        int rank = 6;
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setRank(rank);
        int expResult = rank;
        int result = instance.getRank();
        assertEquals(expResult, result);
    }

    /**
     * Test of integrate method, of class LinearlyMappedBasis.
     */
    public void testIntegrate() {
        System.out.println("integrate");
        int rank = 30;
        LinearlyMappedBasis instance =  new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setRank(rank);
        double[] integrand = new double[rank];
        double[] x = instance.getAbscissas();
        for (int i = 0; i < rank; i++) {
            integrand[i] = Math.exp(x[i]);
        }
        double expResult = Math.exp(1) - Math.exp(0);
        double result = instance.integrate(integrand);
        assertEquals(expResult, result, 0.01);
    }

    /**
     * Test of setRank method, of class LinearlyMappedBasis.
     */
    public void testSetRank() {
        System.out.println("setRank");
        int rank = 5;
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setRank(rank);
    }

    /**
     * Test of setDomain method, of class LinearlyMappedBasis.
     */
    public void testSetDomain() {
        System.out.println("setDomain");
        double[] domain = {Math.random(), Math.random()};
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setDomain(domain);
        assertEquals(domain[0]+domain[1], instance.getDomain()[0] + instance.getDomain()[1], tolerance);
    }

    /**
     * Test of setLeft method, of class LinearlyMappedBasis.
     */
    public void testSetLeft() {
        System.out.println("setLeft");
        double left = Math.random();
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setLeft(left);
        assertEquals(left,instance.getDomain()[0], tolerance);
    }

    /**
     * Test of setRight method, of class LinearlyMappedBasis.
     */
    public void testSetRight() {
        System.out.println("setRight");
        double right = Math.random();
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        instance.setRight(right);
        assertEquals(right,instance.getDomain()[1], tolerance);
    }

    /**
     * Test of getDomain method, of class LinearlyMappedBasis.
     */
    public void testGetDomain() {
        System.out.println("getDomain");
        LinearlyMappedBasis instance = new LinearlyMappedBasis(new ChebyshevExtrema());
        double[] expResult = new double[]{Math.random(), Math.random()};
        instance.setDomain(expResult);
        double[] result = instance.getDomain();
        assertEquals(expResult[0]+expResult[1], result[0] + result[1], tolerance);
    }
}
