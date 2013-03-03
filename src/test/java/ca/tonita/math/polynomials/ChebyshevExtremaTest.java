/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.math.polynomials;

import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class ChebyshevExtremaTest extends TestCase {

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
    public void testGetAbscissas() {
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
}
