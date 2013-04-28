package ca.tonita.math.numerical.spectral;

import ca.tonita.math.linearalgebra.LinearAlgebra;
import ca.tonita.math.polynomials.LinearlyMappedBasis;

/**
 *
 * @author atonita
 */
public class SpectralVector1D {

    private int nDomains;
    /**
     * Indexed in the following way: variables[iDomain][iVariable][iAbscissa]
     */
    private double[][][] variables;
    private double[][][] dvariables;
    private double[] parameters;
    private LinearlyMappedBasis[] bases;

    /**
     * Returns the number of domains.
     *
     * @return the number of domains.
     */
    public int getNDomains() {
        return nDomains;
    }

    public LinearlyMappedBasis[] getBases() {
        return bases;
    }

    /**
     * Creates the spectral vector with no storage.
     *
     * @param domains
     * @param nParameters
     */
    public SpectralVector1D(LinearlyMappedBasis[] bases, int nParameters) {
        this.bases = bases;
        this.parameters = new double[nParameters];
        int nDomains = bases.length;
        variables = new double[nDomains][][];
    }

    /**
     * Computes the derivatives of the variables.
     */
    public void differentiate() {
        for (int iDomain = 0; iDomain < bases.length; iDomain++) {
            double[][] diff = bases[iDomain].getDifferentiationMatrix();
            for (int iVariable = 0; iVariable < variables[iDomain].length; iVariable++) {
                LinearAlgebra.matrixVectorMultiply(diff, variables[iDomain][iVariable], dvariables[iDomain][iVariable]);
            }
        }
    }

    /**
     * Returns the vector of coordinate values for the indexed domain.
     *
     * @param iDomain The domain to get the vertices of.
     * @return the coordinate values of the vertices.
     */
    public double[] getX(int iDomain) {
        return bases[iDomain].getAbscissas();
    }
    
    /**
     * Returns the rank of basis in the given domain.
     * 
     * @param iDomain The index of the domain.
     * @return the rank of that domain
     */
    public int getRank(int iDomain) {
        return bases[iDomain].getRank();
    }

    /**
     * Returns the coordinate values for indexed vertex in the indexed domain.
     *
     * @param iDomain The domain to get the vertices of.
     * @param iX The index of the vertex.
     * @return the coordinate values of the vertex.
     */
    public double getX(int iDomain, int iX) {
        return bases[iDomain].getAbscissas()[iX];
    }

    /**
     * Returns the overall length of the data vector.
     *
     * @return the length of the data vector.
     */
    int getLength() {
        int length = parameters.length;
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            length += bases[iDomain].getRank() * variables[iDomain].length;
        }
        return length;
    }

    /**
     * Returns all the variables at the indexed point.
     *
     * @param iDomain The domain of the point.
     * @param iX The index of the point in the domain.
     * @return The variables at that point.
     */
    double[] getVariables(int iDomain, int iX) {
        int nVariables = variables[iDomain].length;
        double[] y = new double[nVariables];
        for (int iVariable = 0; iVariable < nVariables; iVariable++) {
            y[iVariable] = variables[iDomain][iVariable][iX];
        }
        return y;
    }

    /**
     * Returns all the variables at the indexed point.
     *
     * @param iDomain The domain of the point.
     * @param iX The index of the point in the domain.
     * @return The variables at that point.
     */
    double[] getDVariables(int iDomain, int iX) {
        int nVariables = dvariables[iDomain].length;
        double[] dy = new double[nVariables];
        for (int iVariable = 0; iVariable < nVariables; iVariable++) {
            dy[iVariable] = dvariables[iDomain][iVariable][iX];
        }
        return dy;
    }

    /**
     * Returns the parameters of the vector only.
     * @return the parameters
     */
    double[] getParameters() {
        return parameters;
    }

    /**
     * Returns the variable vector.
     * @param iDomain
     * @param iVariable
     * @return 
     */
    public double[] getVariable(int iDomain, int iVariable) {
        return variables[iDomain][iVariable];
    }
}
