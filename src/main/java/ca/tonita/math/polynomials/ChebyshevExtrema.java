package ca.tonita.math.polynomials;

/**
 * The basis of Chebyshev polynomials using the extrema as the collocation
 * points. Because of some limits, this object represents a truncated basis.
 *
 * @author atonita
 */
public class ChebyshevExtrema {

    /**
     * The locations of the collocation points.
     */
    private double[] abscissas;
    /**
     * The matrix M that converts the vector of values u to coefficients a:
     * Mu=a. See Boyd, Chebyshev and Fourier spectral methods, 5.37.
     */
    private double[][] valuesToCoefficients;
    /**
     * The matrix F that converts the vector of coefficients a to the vector of
     * values at the abscissas u. This is the inverse of the matrix M above.
     */
    private double[][] coefficientsToValues;
    /**
     * The maximum rank of this truncated basis.
     */
    private int rank;

    /**
     * Creates a new ChebyshevExtrema basis.
     */
    public ChebyshevExtrema() {
        rank = 1;
        abscissas = new double[]{0};
    }

    /**
     * Sets the maximum rank of the truncated basis.
     *
     * @param rank the maximum rank of this truncated basis.
     */
    public void setRank(int rank) {
        if (rank <= 0) {
            throw new IllegalArgumentException("Rank must be greater than 0.");
        }
        this.rank = rank;
        this.abscissas = null;
    }

    /**
     * Returns the abscissas for the specified rank.
     *
     * @param rank the desired number of abscissas
     * @return the abscissas for the given rank
     */
    public static double[] getAbscissas(int rank) {
        double[] x = new double[rank];
        double inverseRankMinusOne = 1. / (rank - 1);
        for (int i = 0; i < rank; i++) {
            x[rank - i - 1] = Math.cos(Math.PI * i * inverseRankMinusOne);
        }
        return x;
    }

    /**
     * Returns the abscissas for this basis
     *
     * @return the abscissas
     */
    public double[] getAbscissas() {
        if (abscissas == null) {
            abscissas = getAbscissas(rank);
        }
        return abscissas;
    }

    /**
     * The matrix M that converts the vector of values u to coefficients a:
     * Mu=a. See Boyd, Chebyshev and Fourier spectral methods, eq 5.37 pg 124.
     * The coefficients are also implied by 4.49 and 4.50.
     *
     * @return The matrix which changes basis: function values -> Chebyshev
     * coefficients
     */
    public double[][] getValuesToCoefficientsMatrix() {
        if (valuesToCoefficients == null) {
            double[] x = getAbscissas();
            valuesToCoefficients = new double[rank][rank];
            double invRankMinus1 = 1. / (rank - 1);
            // Each row is multiplied by the weight over the inner product of the basis function.. Pi/(N-1) and Pi/2 respectively (with a couple exceptions...)
            for (int i = 0; i < rank; i++) {
                valuesToCoefficients[0][i] = 2 * invRankMinus1;
                valuesToCoefficients[1][i] = x[i] * 2 * invRankMinus1;
            }
            for (int iRow = 2; iRow < rank; iRow++) {
                for (int iCol = 0; iCol < rank; iCol++) {
                    valuesToCoefficients[iRow][iCol] = 2 * x[iCol] * valuesToCoefficients[iRow - 1][iCol] - valuesToCoefficients[iRow - 2][iCol];
                }
            }
            // The exceptions.
            for (int i = 0; i < rank; i++) {
                valuesToCoefficients[0][i] *= 0.5; // Inner product for first row is twice as much.
                valuesToCoefficients[i][0] *= 0.5; // Weight on end point (not an extrema) is halved.
                valuesToCoefficients[i][rank - 1] *= 0.5; // As other end point because trapezoid rule...
                valuesToCoefficients[rank - 1][i] *= 0.5; // I have no explanation for this. :( It is implied by Boyd's equation 4.49. Likely an artifact of the truncation.
                // This last factor implies that Chebyshev-Lobatto Quadrature is exact only up to rank-2 order polynomials. =/
            }
        }
        return valuesToCoefficients;
    }
    
    public double[][] getCoefficientsToValuesMatrix() {
        if (coefficientsToValues == null) {
            coefficientsToValues = new double[rank][rank];
            double[] x = getAbscissas();
            for (int iRow = 0; iRow < rank; iRow++ ) {
                coefficientsToValues[iRow][0] = 1;
                coefficientsToValues[iRow][1] = x[iRow];
                for (int iCol = 2; iCol < rank; iCol++) {
                    coefficientsToValues[iRow][iCol] = 2*x[iRow]*coefficientsToValues[iRow][iCol-1] - coefficientsToValues[iRow][iCol - 2];
                }
            }
        }
        return coefficientsToValues;
    }

    /**
     * The Chebyshev polynomials follow the 2-part linear recursion relation:
     * T_n = 2xT_{n-1} - T_{n-2}.<br> <br> This function turns that recursion
     * relation into an iteration.
     *
     * @param n The order of the basis function to compute.
     * @param x The value of the domain variable to evalute the basis function
     * at.
     * @return The basis function T_n evaluated at x.
     */
    public static double function(int n, double x) {
        double t0 = 1;
        double t1 = x;
        if (n <= 0) {
            return 1;
        } else {
            for (int i = 1; i < n; i++) {
                double oldt1 = t1;
                t1 = 2 * x * t1 - t0;
                t0 = oldt1;
            }
        }
        return t1;
    }
}
