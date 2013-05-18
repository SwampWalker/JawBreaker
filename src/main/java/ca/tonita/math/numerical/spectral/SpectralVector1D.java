package ca.tonita.math.numerical.spectral;

import ca.tonita.math.linearalgebra.LinearAlgebra;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;
import ca.tonita.math.polynomials.LinearlyMappedBasis;
import java.util.ArrayList;

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
    private ODEIndexer1D indexer;

    SpectralVector1D(NonLinearFirstOrderODESystem system, LinearlyMappedBasis[] bases) {
        this.bases = bases;
        nDomains = system.getNDomains();
        parameters = new double[system.getNParameters()];
        int[] rank = new int[nDomains];
        int[] nVariables = new int[nDomains];
        variables = new double[nDomains][][];
        dvariables = new double[nDomains][][];
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            rank[iDomain] = bases[iDomain].getRank();
            nVariables[iDomain] = system.getNVariables(iDomain);
            variables[iDomain] = new double[nVariables[iDomain]][rank[iDomain]];
            dvariables[iDomain] = new double[nVariables[iDomain]][rank[iDomain]];
        }
        indexer = new ODEIndexer1D(rank, nVariables, system.getNParameters(), nDomains);
    }

    /**
     * Returns the indexer for this vector.
     *
     * @return the ODEIndexer1D for this vector.
     */
    public ODEIndexer1D getIndexer() {
        return indexer;
    }

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
     * Returns the pointer to the variables.
     * @return the variables.
     */
    public double[][][] getVariables() {
        return variables;
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
     *
     * @return the parameters
     */
    public double[] getParameters() {
        return parameters;
    }

    /**
     * Returns the variable vector.
     *
     * @param iDomain
     * @param iVariable
     * @return
     */
    public double[] getVariable(int iDomain, int iVariable) {
        return variables[iDomain][iVariable];
    }

    /**
     * Sets a guess. The guess is indexed in the following fashion:
     * variables[iDomain][iVariable][iAbscissa].
     *
     * @param guess The guess to set.
     * @param parameters The guess at the parameter values.
     */
    void setGuess(double[][][] guess, double[] parameters) {
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            int nVariables = variables[iDomain].length;
            int nAbscissa = bases[iDomain].getRank();
            for (int iVariable = 0; iVariable < nVariables; iVariable++) {
                System.arraycopy(guess[iDomain][iVariable], 0, variables[iDomain][iVariable], 0, nAbscissa);
            }
        }
        if (this.parameters.length > 0) {
            System.arraycopy(parameters, 0, this.parameters, 0, this.parameters.length);
        }
        differentiate();
    }

    /**
     * Updates the data given an update vector. The vector is assumed to be the
     * solution to a Newton-Raphson step trying to solve f(V) = 0, the sign
     * convention used follows the derivation, f(V_0) - JdV = 0 -> f(V_0) = JdV.
     * So the update vector is subtracted from the current data vector.
     *
     * @param update The data to use as an update.
     */
    void update(double[] update) {
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            int nVariables = variables[iDomain].length;
            int nAbscissas = bases[iDomain].getRank();
            for (int iVariable = 0; iVariable < nVariables; iVariable++) {
                for (int iAbscissa = 0; iAbscissa < nAbscissas; iAbscissa++) {
                    int index = indexer.index(iDomain, iVariable, iAbscissa);
                    variables[iDomain][iVariable][iAbscissa] -= update[index];
                }
            }
        }
        for (int iParameter = 0; iParameter < parameters.length; iParameter++) {
            parameters[iParameter] -= update[indexer.index(iParameter)];
        }
        differentiate();
    }

    /**
     * Sets a guess. The guess is indexed in the following fashion:
     * guess[iDomain].get(iAbscissa)[iVariable].
     *
     * @param guess The guess to set.
     * @param parameters The guess at the parameter values.
     */
    void setGuess(ArrayList<double[]>[] guess, double[] parameters) {
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            int nVariables = guess[iDomain].get(0).length;
            int nAbscissa = bases[iDomain].getRank();
            for (int iVariable = 0; iVariable < nVariables; iVariable++) {
                for (int iAbscissa = 0; iAbscissa < nAbscissa; iAbscissa++) {
                    variables[iDomain][iVariable][iAbscissa] = guess[iDomain].get(iAbscissa)[iVariable];
                }
            }
        }
        if (this.parameters.length > 0) {
            System.arraycopy(parameters, 0, this.parameters, 0, this.parameters.length);
        }
        differentiate();
    }
    
    /**
     * Prints variables.
     */
    void printVariables() {
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            int nVariables = variables[iDomain].length;
            int nAbscissa = bases[iDomain].getRank();
            for (int iVariable = 0; iVariable < nVariables; iVariable++) {
                for (int iAbscissa = 0; iAbscissa < nAbscissa; iAbscissa++) {
                    System.out.print(variables[iDomain][iVariable][iAbscissa] + " ");
                }
                System.out.print("\n");
            }
        }
        for (int iDomain = 0; iDomain < nDomains; iDomain++) {
            int nVariables = variables[iDomain].length;
            int nAbscissa = bases[iDomain].getRank();
            for (int iVariable = 0; iVariable < nVariables; iVariable++) {
                for (int iAbscissa = 0; iAbscissa < nAbscissa; iAbscissa++) {
                    System.out.print(dvariables[iDomain][iVariable][iAbscissa] + " ");
                }
                System.out.print("\n");
            }
        }
    }
}
