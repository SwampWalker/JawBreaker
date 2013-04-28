package ca.tonita.math.numerical.spectral;

import ca.tonita.math.linearalgebra.LinearAlgebra;
import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;
import ca.tonita.math.polynomials.LinearlyMappedBasis;

/**
 *
 * @author atonita
 */
public class SpectralSolver1D {

    private NonLinearFirstOrderODESystem system;
    private SpectralVector1D vector;
    private LinearlyMappedBasis[] bases;
    private double[][] jacobian;
    private double[] rhs;

    public SpectralSolver1D(NonLinearFirstOrderODESystem system, LinearlyMappedBasis[] bases) {
        this.system = system;
        this.bases = bases;
        vector = new SpectralVector1D(system, bases);
    }

    public void solve() {
        for (int i = 0; i < 1; i++) {
            computeNewtonRaphsonTerms();
            LinearAlgebra.solve(jacobian, rhs);
            vector.update(rhs);
        }
    }

    /**
     * Sets a guess. The guess is indexed in the following fashion:
     * variables[iDomain][iVariable][iAbscissa].
     *
     * @param guess The guess to set.
     * @param parameters The guess at the parameter values.
     */
    public void setGuess(double[][][] guess, double[] parameters) {
        vector.setGuess(guess, parameters);
    }

    /**
     * Computes the jacobian and the right hand side.
     */
    private void computeNewtonRaphsonTerms() {
        int N = vector.getLength();
        jacobian = new double[N][N];
        rhs = new double[N];
        for (int i = 0; i < N; i++) {
            rhs[i] = 0;
            for (int j = 0; j < N; j++) {
                jacobian[i][j] = 0;
            }
        }
        double[] parameters = vector.getParameters();
        // Compute.
        ODEIndexer1D indexer = vector.getIndexer();
        for (int iDomain = 0; iDomain < vector.getNDomains(); iDomain++) {
            double[] x = vector.getX(iDomain);
            double[][] diff = vector.getBases()[iDomain].getDifferentiationMatrix();
            for (int iX = 0; iX < x.length; iX++) {
                int type = NonLinearFirstOrderODESystem.BULK;
                if (iX == 0) {
                    type = NonLinearFirstOrderODESystem.LEFTBOUNDARY;
                } else if (iX == x.length - 1) {
                    type = NonLinearFirstOrderODESystem.RIGHTBOUNDARY;
                }
                double[] y = vector.getVariables(iDomain, iX);
                double[] dy = vector.getDVariables(iDomain, iX);
                double[] d = diff[iX]; // The vector that dotted into a variable vector gives the derivative at the iX'th collocation. 
                NonLinearFirstOrderODEBean equations = system.equations(iX, x[iX], y, dy, parameters, type);
                int nVariables = indexer.getNVariables(iDomain);
                for (int iEquation = 0; iEquation < nVariables; iEquation++) {
                    int iRow = indexer.index(iDomain, iEquation, iX);
                    rhs[iRow] = equations.getResidue()[iEquation];

                    /* Indexing the jacobian is the complicated part. :)
                     * 
                     * The idea is that we have to iterate over the columns,
                     * which is an iteration over the variables in the same
                     * domain, and for each variable, an iteration over the
                     * abscissas.
                     */
                    for (int jVariable = 0; jVariable < nVariables; jVariable++) {
                        // Direct part first.
                        jacobian[iRow][indexer.index(iDomain, jVariable, iX)] += equations.getJacobian()[iEquation][jVariable];
                        // Derivative part.
                        double jac = equations.getJacobian()[iEquation][nVariables + jVariable];
                        for (int jX = 0; jX < x.length; jX++) {
                            int iColumn = indexer.index(iDomain, jVariable, jX);
                            jacobian[iRow][iColumn] += jac * d[jX];
                        }
                    }
                    // Parameter part.
                    for (int iParameter = 0; iParameter < parameters.length; iParameter++) {
                        int index = indexer.index(iParameter);
                        jacobian[iRow][index] += equations.getJacobian()[iEquation][nVariables * 2 + iParameter];
                    }
                }
            } // End loop over equations.
        } // End loop over domains.
        // The constraints.
        for (int iConstraint = 0; iConstraint < parameters.length; iConstraint++) {
            int iRow = indexer.index(iConstraint);
            rhs[iRow] = system.computeConstraint(iConstraint, indexer, vector, jacobian[iRow]);
        }
    }

    public double[][][] getSolution() {
        return vector.getVariables();
    }
}
