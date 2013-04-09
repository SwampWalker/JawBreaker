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
            }
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

    /**
     * Test of getCoefficientsToValuesMatrix method, of class ChebyshevExtrema.
     */
    public void testGetCoefficientsToValuesMatrix() {
        System.out.println("getCoefficientsToValuesMatrix");
        ChebyshevExtrema instance = new ChebyshevExtrema();
        int rank = 6;
        instance.setRank(rank);
        double[][] result = instance.getCoefficientsToValuesMatrix();
        double[] x = instance.getAbscissas();
        for (int i = 0; i < rank; i++) {
            double[] coeffs = new double[rank];
            double[] u = new double[rank];
            for (int j = 0; j < rank; j++) {
                coeffs[j] = 0;
                u[j] = ChebyshevExtrema.function(i, x[j]);
            }
            coeffs[i] = 1.0;
            double[] u2 = LinearAlgebra.matrixVectorMultiply(result, coeffs);
            for (int j = 0; j < rank; j++) {
                if (Math.abs(u2[j] - u[j]) > tolerance) {
                    for (int k = 0; k < rank; k++) {
                        System.out.println(x[k] + "\t" + u[k] + "\t" + u2[k]);
                    }
                    fail("Tolerance failure in coefficient " + i + " at abscissas " + j + ".");
                }
            }
        }
    }

    /**
     * Test of getCoefficientDifferentiationMatrix method, of class ChebyshevExtrema.
     */
    public void testGetCoefficientDifferentiationMatrix() {
        System.out.println("getCoefficientDifferentiationMatrix");
        ChebyshevExtrema instance = new ChebyshevExtrema();
        int rank = 5;
        instance.setRank(rank);
        double[][] result = instance.getCoefficientDifferentiationMatrix();
        for (int i = 0; i < rank; i++) {
            for (int j = 0; j <= i; j++) {
                if (Math.abs(result[i][j]) > tolerance) {
                    fail("Non zero entry found in position (" + i + "," + j + ").");
                }
                if (i < rank - 1 && i > 0) {
                    if (Math.abs(result[i][i+1] - (i+1)*2) > tolerance) {
                        fail("Super diagonal entry " + i + " not equal to " + (2*i+2) + ", received " + result[i][i+1]);
                    }
                }
            }
        }
        // We will test the rest implicitly with the testGetDifferentiationMatrix() test.
    }

    /**
     * Test of getDifferentiationMatrix method, of class ChebyshevExtrema.
     */
    public void testGetDifferentiationMatrix() {
        System.out.println("getDifferentiationMatrix");
        ChebyshevExtrema instance = new ChebyshevExtrema();
        instance.setRank(20);
        double[][] result = instance.getDifferentiationMatrix();
        double[] x = instance.getAbscissas();
        double[] u = new double[instance.getRank()];
        for (int i = 0; i < x.length; i++) {
            u[i] = Math.sin(x[i]);
        }
        double[] du = LinearAlgebra.matrixVectorMultiply(result, u);
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(du[i] - Math.cos(x[i])) > tolerance) {
                for (int j = 0; j < x.length; j++) {
                    System.out.println(i + ": " + x[j] + " " + du[j] + " " + Math.cos(x[j]) );
                }
                fail("Tolerance failure at point " + i + ".");
            }
        }
    }

    /**
     * Test of getRank method, of class ChebyshevExtrema.
     */
    public void testGetRank() {
        System.out.println("getRank");
        ChebyshevExtrema instance = new ChebyshevExtrema();
        int expResult = 5;
        instance.setRank(expResult);
        int result = instance.getRank();
        assertEquals(expResult, result);
    }

    /**
     * Test of integrate method, of class ChebyshevExtrema.
     */
    public void testIntegrate() {
        System.out.println("integrate");
        int rank = 20;
        double[] integrand = new double[rank];
        ChebyshevExtrema instance = new ChebyshevExtrema();
        instance.setRank(rank);
        double[] x = instance.getAbscissas();
        for (int i = 0; i < rank; i++) {
            integrand[i] = df(x[i]);
        }
        double expResult = f(1) - f(-1);
        double result = instance.integrate(integrand);
        assertEquals(expResult, result, tolerance);
    }
    
    private double f(double x) {
        return 1.3*Math.pow(x, 12);
    }
    
    private double df(double x) {
        return 1.3*12*Math.pow(x, 11);
    }
}
