/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.math.polynomials;

import ca.tonita.math.numerical.LinearAlgebra;
import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class ChebyshevExtremaTest extends TestCase {
    
    private static double tolerance = 1.0e-10;

    public ChebyshevExtremaTest(String testName) {
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
     * Test of setRank method, of class ChebyshevExtrema.
     */
    public void testSetRank() {
        System.out.println("setRank");
        int rank = 0;
        ChebyshevExtrema instance = new ChebyshevExtrema();

        // Zero rank not allowed.
        boolean rankZeroNotAllowed = false;
        try {
            instance.setRank(rank);
        } catch (IllegalArgumentException e) {
            rankZeroNotAllowed = true;
        }
        assertEquals(rankZeroNotAllowed, true);

    }

    /**
     * Test of getAbscissas method, of class ChebyshevExtrema.
     */
    public void testGetAbscissas_0args() {
        System.out.println("getAbscissas");
        int rank = 3;
        ChebyshevExtrema instance = new ChebyshevExtrema();
        double[] expResult = {-1, 0, 1};
        double[] result = instance.getAbscissas(rank);
        boolean allAbscissasEqual = true;
        for (int i = 0; i < rank; i++) {
            if (Math.abs(result[i] - expResult[i]) > 1.0e-15) {
                allAbscissasEqual = false;
            }
        }
        assertEquals(allAbscissasEqual, true);
    }

    /**
     * Test of getAbscissas method, of class ChebyshevExtrema.
     */
    public void testGetAbscissas_int() {
        System.out.println("getAbscissas");
        // The 0 args version test implicitly tests this function.
    }

    /**
     * Test of getValuesToCoefficientsMatrix method, of class ChebyshevExtrema.
     */
    public void testGetValuesToCoefficientsMatrix() {
        System.out.println("getValuesToCoefficientsMatrix");
        ChebyshevExtrema instance = new ChebyshevExtrema();
        int rank = 6;
        instance.setRank(rank);
        double[] x = instance.getAbscissas();
        double[][] M = instance.getValuesToCoefficientsMatrix();
        // We won't analyse M directly. Instead we will ensure that it can
        // distinguish the basis functions directly...
        for (int i = 0; i < rank; i++) {
            double[] u = new double[rank];
            for (int j = 0; j < rank; j++) {
                u[j] = instance.function(i, x[j]);
                System.out.println(u[j] + " " + x[j]);
            }
            System.out.print("\n");
            double[] a = LinearAlgebra.matrixVectorMultiply(M, u);
            for (int j = 0; j < rank; j++) {
                if (j == i) {
                    if (Math.abs(a[j] - 1) > tolerance) {
                        fail("Expected to recover coefficient 1 for basis function " + i + " received coefficient " + a[j] + ".");
                    }
                } else {
                    if (Math.abs(a[j]) > tolerance) {
                        fail("Expected to recover coefficient 0 for basis function " + i + " in position " + j + " received coefficient " + a[j] + ".");
                    }
                }
            }
        }
    }

    /**
     * Test of function method, of class ChebyshevExtrema.
     */
    public void testFunction() {
        System.out.println("function");
        int n = 5;
        double[] x = ChebyshevExtrema.getAbscissas(n);
        double[] expResult = new double[]{1,-1,1,-1,1};
        for (int i = 0; i < n; i++) {
            double result = ChebyshevExtrema.function(n-1, x[i]);
            if (Math.abs(result - expResult[i]) > tolerance) {
                    fail("Failure at extrema " + i + ", expected "  + expResult[i] + " received " + result + ".");
            }
        }
    }
}
