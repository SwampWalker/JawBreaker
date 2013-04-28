package ca.tonita.math.numerical.spectral;

import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;

/**
 *
 * @author atonita
 */
public class SpectralSolver1D {

    private NonLinearFirstOrderODESystem system;
    private SpectralVector1D vector;
    private ODEIndexer1D indexer;
    private double[][] jacobian;
    private double[] rhs;

    public SpectralSolver1D(NonLinearFirstOrderODESystem system, SpectralVector1D vector, ODEIndexer1D indexer) {
        this.system = system;
        this.vector = vector;
        this.indexer = indexer;
    }

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
        for (int iDomain = 0; iDomain <vector.getNDomains(); iDomain++) {
            double[] x = vector.getX(iDomain);
            double[][] diff = vector.getBases()[iDomain].getDifferentiationMatrix();
            for (int iX = 0; iX < x.length; iX++) {
                double[] y = vector.getVariables(iDomain, iX);
                double[] dy = vector.getDVariables(iDomain, iX);
                double[] d = diff[iX]; // The vector that dotted into a variable vector gives the derivative at the iX'th collocation. 
                NonLinearFirstOrderODEBean equations = system.equations(iX, x[iX], y, dy, parameters);
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
                            int iColumn = indexer.index(iDomain, iEquation, iX);
                            jacobian[iRow][iColumn] += jac * d[jX];
                        }
                    }
                    // Parameter part.
                    for (int iParameter = 0; iParameter < parameters.length; iParameter++) {
                        jacobian[iRow][indexer.index(nVariables)] += equations.getJacobian()[iEquation][nVariables * 2 + iParameter];
                    }
                }
            } // End loop over equations.
        } // End loop over domains.
        // The constraints.
        for (int iConstraint = 0; iConstraint < parameters.length; iConstraint++) {
            system.computeConstraint(iConstraint, indexer, vector, jacobian[indexer.index(iConstraint)]);
        }
    }
}
