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
    double[] getAbscissas();

    /**
     * Returns the differentiation matrix D such that for a function u, Du is
     * the approximation of the derivative.
     *
     * @return the differentiation operator
     */
    double[][] getDifferentiationMatrix();

    /**
     * Returns the rank of the truncated basis.
     *
     * @return the rank
     */
    int getRank();

    /**
     * Computes the integral of a function using Guassian quadrature.
     *
     * @param integrand The integrand to integrate, must be at least length
     * rank.
     * @return the integrated function.
     */
    double integrate(double[] integrand);

    /**
     * Sets the maximum rank of the truncated basis.
     *
     * @param rank the maximum rank of this truncated basis.
     */
    void setRank(int rank);

    /**
     * Returns a two element array with the left and right end points of the
     * domain over which the basis is orthogonal with respect to a given weight.
     *
     * @return the domain of the polynomial basis
     */
    public double[] getDomain();
}
