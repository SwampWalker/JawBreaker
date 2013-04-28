package ca.tonita.math.linearalgebra;

public final class LinearAlgebra {

    /**
     * Computes the LU factorization of A, overwriting the contents. The lower
     * sub-diagonal contains L, the upper triangle contains U.
     *
     * This function assumes that the input matrix is square.
     *
     * @param A The square matrix to factor.
     * @return error code, 0 if everything was ok.
     */
    public static int LUFactorNoPivot(double[][] A) {
        for (int k = 0; k < A.length; k++) {
            for (int j = k + 1; j < A.length; j++) {
                A[j][k] = A[j][k] / A[k][k];
                for (int i = k + 1; i < A.length; i++) {
                    A[j][i] -= A[j][k] * A[k][i];
                }
            }
        }
        return 0;
    }

    /**
     * Factors A into LU form with pivoting. The pivot vector is returned.
     *
     * @param A The square matrix A.
     * @return The pivot vector.
     */
    public static int[] LUFactor(double[][] A) {
        int N = A.length;
        int[] pivots = new int[N];
        for (int i = 0; i < N; i++) {
            pivots[i] = i;
        }
        for (int k = 0; k < N; k++) {
            // Find the pivot.
            int iMax = k;
            double max = Math.abs(A[k][k]);
            for (int i = k; i < N; i++) {
                if (Math.abs(A[i][k]) > max) {
                    iMax = i;
                    max = Math.abs(A[i][k]);
                }
            }

            // Swap rows.
            double[] rowMax = A[iMax];
            A[iMax] = A[k];
            A[k] = rowMax;
            // Note this in pivots.
            int row = pivots[k];
            pivots[k] = pivots[iMax];
            pivots[iMax] = row;

            // Row reduce
            for (int j = k + 1; j < N; j++) {
                A[j][k] = A[j][k] / A[k][k];
                for (int i = k + 1; i < N; i++) {
                    A[j][i] -= A[j][k] * A[k][i];
                }
            }
        }

        return pivots;
    }

    /**
     * Solves the system Ax = b, overwriting A and b with the LU factorization
     * of A and the solution x respectively.
     *
     * @param A The square matrix, A.
     * @param b The right hand side vector.
     * @return The status, 0 indicating success.
     */
    public static int solve(double[][] A, double[] b) {
        int[] pivots = LUFactor(A);

        // now we need to change the ordering of b.
        double[] bCopy = b.clone();
        for (int i = 0; i < b.length; i++) {
            b[i] = bCopy[pivots[i]];
        }

        forwardSubstitution(A, b);
        backwardSubstitution(A, b);

        return 0;
    }

    /**
     * Forwards substitution. Solves the system Lx = b, where L is lower
     * triangular.
     *
     * @param L A lower triangular matrix (should still be double[N][N]).
     * @param b The right hand side.
     */
    private static void forwardSubstitution(double[][] L, double[] b) {
        for (int iRow = 0; iRow < L.length; iRow++) {
            for (int iCol = 0; iCol < iRow; iCol++) {
                b[iRow] -= b[iCol] * L[iRow][iCol];
            }
        }
    }

    /**
     * Backwards substitution. Solves the system Ux = b, where U is upper
     * triangular.
     *
     * @param U An upper triangular matrix (should still be double[N][N]).
     * @param b The right hand side.
     */
    private static void backwardSubstitution(double[][] U, double[] b) {
        for (int iRow = U.length - 1; iRow >= 0; iRow--) {
            for (int iCol = U.length - 1; iCol >= iRow + 1; iCol--) {
                b[iRow] -= b[iCol] * U[iRow][iCol];
            }
            b[iRow] /= U[iRow][iRow];
        }
    }

    public static void printMatrix(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a.length; j++) {
                System.out.print(a[i][j] + " ");
            }
            System.out.print("\n");
        }
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
    
    /**
     * Multiplies a vector into a matrix: Ax = b. Returns b.
     * @param A The matrix to multiply into.
     * @param x The vector to multiply with.
     * @param b The result.
     */
    public static void matrixVectorMultiply(double[][] A, double[] x, double[] b) {
        int N = A.length;
        int M = x.length;
        for (int i = 0; i < N; i++) {
            b[i] = 0;
            for (int j = 0; j < M; j++) {
                b[i] += A[i][j]*x[j];
            }
        }
    }
}
