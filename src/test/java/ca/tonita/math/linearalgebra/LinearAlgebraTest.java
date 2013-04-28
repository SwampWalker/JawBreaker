package ca.tonita.math.linearalgebra;

import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class LinearAlgebraTest extends TestCase {
    
    private static double tolerance = 1.0e-10;
    
    public LinearAlgebraTest(String testName) {
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
     * Test of factorizeLU method, of class LinearAlgebra.
     */
    public void testFactorizeLU() {
        System.out.println("factorizeLU");
        double[][] A = new double[4][];
        A[0] = new double[]{2,1,1,0};
        A[1] = new double[]{4,3,3,1};
        A[2] = new double[]{8,7,9,5};
        A[3] = new double[]{6,7,9,8};
        int[] expResult = {2,3,1,0};
        int[] result = LinearAlgebra.LUFactor(A);
        double[][] LU = new double[4][];
        LU[0] = new double[]{8,7,9,5};
        LU[1] = new double[]{3./4,7./4,9./4,17./4};
        LU[2] = new double[]{1./2,-2./7,-6./7,-2./7};
        LU[3] = new double[]{1./4,-3./7,1./3,2./3};
        for (int i = 0; i < 4; i++) {
            // Test pivot.
            if (result[i] != expResult[i]) {
                fail("Incorrect pivot in row " + i + ", expected " + expResult[i] + " received " + result[i] + ".");
            }
            // Test rows of LU.
            for (int j = 0; j < 4; j++) {
                if (Math.pow(A[i][j] - LU[i][j], 2) > tolerance*tolerance) {
                    fail("Incorrect value of LU factors at position (" + i + "," + j + "). Expected " + LU[i][j] + " received " + A[i][j] + ".");
                }
            }
        }
    }

    /**
     * Test of solve method, of class LinearAlgebra.
     */
    public void testSolve() {
        System.out.println("solve");
        double[][] A = new double[4][];
        A[0] = new double[]{2,1,1,0};
        A[1] = new double[]{4,3,3,1};
        A[2] = new double[]{8,7,9,5};
        A[3] = new double[]{6,7,9,8};
        double[] expResult = {1,2,3,4};
        double[] b = new double[]{2+2+3,4+3*2+3*3+4,8+7*2+9*3+5*4,6+7*2+9*3+8*4};
        int result = LinearAlgebra.solve(A, b);
        for (int i = 0; i < 4; i++) {
            if (Math.pow(b[i] - expResult[i], 2) > tolerance*tolerance) {
                fail("Incorrect value at row " + i + ", expected " + expResult[i] + " receive " + b[i] + ".");
            }
        }
    }

    /**
     * Test of matrixMultiply method, of class LinearAlgebra.
     */
    public void testMatrixMultiply() {
        System.out.println("matrixMultiply");
        double[][] A = new double[2][];
        A[0] = new double[]{1,2,3};
        A[1] = new double[]{4,5,6};
        double[][] B = new double[3][];
        B[0] = new double[]{7,0};
        B[1] = new double[]{8,1};
        B[2] = new double[]{9,2};
        double[][] expResult = new double[2][];
        expResult[0] = new double[]{50, 8};
        expResult[1] = new double[]{122,17};
        double[][] result = LinearAlgebra.matrixMultiply(A, B);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if (Math.pow(expResult[i][j] - result[i][j], 2) > tolerance*tolerance) {
                    fail("Bad value at (" + i + "," + j + "), expected " + expResult[i][j] + " received " + result[i][j] + ".");
                }
            }
        }
    }

    /**
     * Test of matrixVectorMultiply method, of class LinearAlgebra.
     */
    public void testMatrixVectorMultiply() {
        System.out.println("matrixVectorMultiply");
        double[][] A = new double[4][];
        A[0] = new double[]{2,1,1,0};
        A[1] = new double[]{4,3,3,1};
        A[2] = new double[]{8,7,9,5};
        A[3] = new double[]{6,7,9,8};
        double[] x = {1,2,3,4};
        double[] expResult = new double[]{2+2+3,4+3*2+3*3+4,8+7*2+9*3+5*4,6+7*2+9*3+8*4};
        double[] result = LinearAlgebra.matrixVectorMultiply(A, x);
        for (int i = 0; i < 4; i++) {
            if (Math.pow(result[i] - expResult[i], 2) > tolerance*tolerance) {
                fail("Bad value at row " + i + ", expected " + expResult[i] + " received " + result[i] + ".");
            }
        }
    }
}
