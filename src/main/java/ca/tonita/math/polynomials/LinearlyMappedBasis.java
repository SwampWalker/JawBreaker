package ca.tonita.math.polynomials;

/**
 *
 * @author atonita
 */
public class LinearlyMappedBasis implements PolynomialBasis {

    private PolynomialBasis basis;
    private double[] mappedAbscissas;
    private double[][] mappedDifferentiate;
    private double[] domain = {0, 1};
    private double slope;
    
    public LinearlyMappedBasis(PolynomialBasis basis) {
        this.basis = basis;
        slope = 1./(basis.getDomain()[1] - basis.getDomain()[0]);
    }

    @Override
    public double[] getAbscissas() {
        if (mappedAbscissas == null) {
            double[] x = basis.getAbscissas();
            mappedAbscissas = new double[basis.getRank()];
            for (int i = 0; i < getRank(); i++) {
                mappedAbscissas[i] = map(x[i]);
            }
        }
        return mappedAbscissas;
    }

    @Override
    public double[][] getDifferentiationMatrix() {
        if (mappedDifferentiate == null) {
            mappedDifferentiate = new double[getRank()][getRank()];
            double[][] differentiate = basis.getDifferentiationMatrix();
            for (int i = 0; i < getRank(); i++) {
                for (int j = 0; j < getRank(); j++) {
                    mappedDifferentiate[i][j] = differentiate[i][j]/dmap();
                }
            }
        }
        return mappedDifferentiate;
    }

    @Override
    public int getRank() {
        return basis.getRank();
    }

    @Override
    public double integrate(double[] integrand) {
        return basis.integrate(integrand)*dmap();
    }

    @Override
    public void setRank(int rank) {
        basis.setRank(rank);
        mappedAbscissas = null;
        mappedDifferentiate = null;
    }

    /**
     * Sets the domain.
     *
     * @param domain the domain to set
     */
    public void setDomain(double[] domain) {
        this.domain[0] = domain[0];
        this.domain[1] = domain[1];
    }

    /**
     * Sets the left endpoint.
     *
     * @param left the left end point of the domain.
     */
    public void setLeft(double left) {
        domain[0] = left;
    }

    /**
     * Sets the right endpoint.
     *
     * @param right the right end point of the domain
     */
    public void setRight(double right) {
        domain[1] = right;
    }

    @Override
    public double[] getDomain() {
        return domain;
    }

    private double map(double x) {
        return domain[0] + (domain[1] - domain[0])*(x - basis.getDomain()[0])*slope;
    }
    
    private double dmap() {
        return (domain[1] - domain[0])*slope;
    }
}
