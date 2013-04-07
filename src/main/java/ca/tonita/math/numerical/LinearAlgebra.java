package ca.tonita.math.numerical;

/**
 *
 * @author atonita
 */
public class LinearAlgebra {

    /**
     * Performs an LU factorization with partial pivoting. On exit, A is
     * overwritten with the LU decomposition of A. The pivot matrix is returned
     * in a permutation vector.
     *
     * @param A The matrix to factorize.
     * @return The permutation vector.
     */
    public static int[] factorizeLU(double[][] A) {
        int N = A.length;
        int[] pivots = new int[N];
        for (int i = 0; i < N; i++) {
            pivots[i] = i;
        }
        for (int iRow = 0; iRow < N - 1; iRow++) {
            // Find pivot.
            int iSwap = iRow;
            double pivotValue = A[iSwap][iSwap] * A[iSwap][iSwap];
            for (int iPivot = iRow + 1; iPivot < N; iPivot++) {
                double potentialPivotValue = A[iPivot][iRow] * A[iPivot][iRow];
                if (potentialPivotValue > pivotValue) {
                    pivotValue = potentialPivotValue;
                    iSwap = iPivot;
                }
            }
            // Swap rows.
            double[] placeHolder = A[iSwap];
            A[iSwap] = A[iRow];
            A[iRow] = placeHolder;
            int c = pivots[iSwap];
            pivots[iSwap] = pivots[iRow];
            pivots[iRow] = c;

            // Factorize.
            for (int iRow2 = iRow + 1; iRow2 < N; iRow2++) {
                // L part.
                A[iRow2][iRow] /= A[iRow][iRow];
                // U part.
                for (int iColumn = iRow + 1; iColumn < N; iColumn++) {
                    A[iRow2][iColumn] -= A[iRow2][iRow] * A[iRow][iColumn];
                }
            }
        }
        return pivots;
    }

    /**
     * Solves the equation Ax = b. On exit, A is replaced with it's pivoted LU
     * factorization which is essentially useless without the pivots...
     *
     * @param A The matrix.
     * @param b The right hand side vector.
     * @return The solution vector x.
     */
    public static double[] solve(double[][] A, double[] b) {
        // This turns Ax = b -> PAx = LUx = Pb
        int[] pivots = factorizeLU(A);
        // Pivot b, that is compute Pb.
        int N = b.length;
        double[] x = new double[N];
        for (int i = 0; i < N; i++) {
            x[i] = b[pivots[i]];
        }

        // Forward substitute to find y: LUx = Ly = Pb
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                x[i] -= x[j]*A[i][j];
            }
        }
        
        // Now, variable x stores y, we back substitute for actual x: Ux = y
        for (int i = N-1; i >= 0; i--) {
            for (int j = i + 1; j < N; j++) {
                x[i] -= x[j]*A[i][j];
            }
            x[i] /= A[i][i];
        }
        return x;
    }
    
    /**
     * Multiplies two matrices together, returning the result: AB = C.
     * @param A The left hand matrix.
     * @param B The right hand matrix.
     * @return The product.
     */
    public static double[][] matrixMultiply(double[][] A, double[][] B) {
        int N = A.length;
        int M = B[0].length;
        int O = B.length;
        double[][] C = new double[N][M];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                C[i][j] = 0;
                for (int k = 0; k < O; k++) {
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }
        return C;
    }
    
    /**
     * Multiplies a vector into a matrix: Ax = b. Returns b.
     * @param A The matrix to multiply into.
     * @param x The vector to multiply with.
     * @return The result.
     */
    public static double[] matrixVectorMultiply(double[][] A, double[] x) {
        int N = A.length;
        int M = x.length;
        double[] b = new double[N];
        for (int i = 0; i < N; i++) {
            b[i] = 0;
            for (int j = 0; j < M; j++) {
                b[i] += A[i][j]*x[j];
            }
        }
        return b;
    }
}
