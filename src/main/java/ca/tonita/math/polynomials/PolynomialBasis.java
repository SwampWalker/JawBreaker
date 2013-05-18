package ca.tonita.math.polynomials;

/**
 *
 * @author atonita
 */
public interface PolynomialBasis {

    /**
     * Returns the abscissas for this basis
     *
     * @return the abscissas
     */
    public double[] getAbscissas();

    /**
     * Returns the abscissas for this basis for the specified rank.
     *
     * @param n The number of the abscissas desired abscissas.
     * @return the abscissas
     */
    public double[] getAbscissas(int n);

    /**
     * Returns the differentiation matrix D such that for a function u, Du is
     * the approximation of the derivative.
     *
     * @return the differentiation operator
     */
    public double[][] getDifferentiationMatrix();

    /**
     * Returns the rank of the truncated basis.
     *
     * @return the rank
     */
    public int getRank();

    /**
     * Computes the integral of a function using Guassian quadrature.
     *
     * @param integrand The integrand to integrate, must be at least length
     * rank.
     * @return the integrated function.
     */
    public double integrate(double[] integrand);

    /**
     * Sets the maximum rank of the truncated basis.
     *
     * @param rank the maximum rank of this truncated basis.
     */
    public void setRank(int rank);

    /**
     * Returns a two element array with the left and right end points of the
     * domain over which the basis is orthogonal with respect to a given weight.
     *
     * @return the domain of the polynomial basis
     */
    public double[] getDomain();
    
    /**
     * Interpolates the function onto the given points.
     * @param function the function to interpolate
     * @param x the points to interpolate the function onto
     * @return the interpolated values
     */
    public double[] interpolate(double[] function, double[] x);
}
