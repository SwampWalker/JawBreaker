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
     * @param rank the maximum rank of this truncated basis.
     */
    public void setRank(int rank) {
        if (rank <= 0) {
            throw new IllegalArgumentException("Rank must be greater than 0.");
        }
        this.rank = rank;
        this.abscissas = getAbscissas(rank);
    }
    
    /**
     * Returns the abscissas for the specified rank.
     * @param rank the desired number of abscissas
     * @return the abscissas for the given rank
     */
    public double[] getAbscissas(int rank) {
        double[] x = new double[rank];
        double inverseRankMinusOne = 1./(rank - 1);
        for (int i = 0; i < rank; i++) {
            x[rank-i-1] = Math.cos(Math.PI*i*inverseRankMinusOne);
        }
        return x;
    }
    
    /**
     * Returns the abscissas for this basis
     * @return the abscissas
     */
    public double[] getAbscissas() {
        return this.abscissas;
    }
}
