package ca.tonita.math.numerical.spectral;

import ca.tonita.math.linearalgebra.LinearAlgebra;
import ca.tonita.math.numerical.NewtonRaphsonPostProcessor;
import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;
import ca.tonita.math.polynomials.LinearlyMappedBasis;
import java.util.ArrayList;

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

    /**
     * Solves the Nonlinear system by successive Newton-Raphson iterations. If
     * the system implements the NewtonRaphsonPostProcessor interface, this is
     * called after each step.
     *
     * @param tolerance The tolerance to attempt to achieve in the Euclidean
     * norm of the residual vector.
     * @param maxSteps The maximum number of NR steps to perform.
     */
    public void solve(double tolerance, int maxSteps) {
        double residue = 2 * tolerance;
        for (int i = 0; i < maxSteps && residue > tolerance * 0; i++) {
            residue = computeNewtonRaphsonTerms();
            System.out.println(i + " " + residue);
            LinearAlgebra.solve(jacobian, rhs);
            vector.update(rhs);
            if (system instanceof NewtonRaphsonPostProcessor) {
                ((NewtonRaphsonPostProcessor) system).postStepProcessing(vector, bases);
            }
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
     * Computes a finite difference version of the Jacobian. For testing
     * purposes.
     */
    public void finiteDifferenceJacobian() {
        int N = vector.getLength();
        double[][] rhspm = new double[2][N];
        double[] update = new double[N];
        double h = 1.0e-6;
        ODEIndexer1D indexer = vector.getIndexer();
        double[][] fdJacobian = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int iRHS = 0; iRHS < 2; iRHS++) {
                update[i] = iRHS == 0 ? h : -2 * h;
                vector.update(update);
                if (system instanceof NewtonRaphsonPostProcessor) {
                    ((NewtonRaphsonPostProcessor) system).postStepProcessing(vector, bases);
                }
                double[] parameters = vector.getParameters();
                for (int iDomain = 0; iDomain < vector.getNDomains(); iDomain++) {
                    double[] x = vector.getX(iDomain);
                    for (int iX = 0; iX < x.length; iX++) {
                        int type = NonLinearFirstOrderODESystem.BULK;
                        if (iX == 0) {
                            type = NonLinearFirstOrderODESystem.LEFTBOUNDARY;
                        } else if (iX == x.length - 1) {
                            type = NonLinearFirstOrderODESystem.RIGHTBOUNDARY;
                        }
                        double[] y = vector.getVariables(iDomain, iX);
                        double[] dy = vector.getDVariables(iDomain, iX);
                        NonLinearFirstOrderODEBean equations = system.equations(iX, x[iX], y, dy, parameters, type);
                        int nVariables = indexer.getNVariables(iDomain);
                        for (int iEquation = 0; iEquation < nVariables; iEquation++) {
                            int iRow = indexer.index(iDomain, iEquation, iX);
                            rhspm[iRHS][iRow] = equations.getResidue()[iEquation];
                        }
                    }
                }
                for (int iConstraint = 0; iConstraint < parameters.length; iConstraint++) {
                    int iRow = indexer.index(iConstraint);
                    double[] junk = new double[N];
                    rhspm[iRHS][iRow] = system.computeConstraint(iConstraint, indexer, vector, junk);
                }
            }
            update[i] = h;
            vector.update(update);
            update[i] = 0;
            for (int j = 0; j < N; j++) {
                fdJacobian[j][i] = -(rhspm[0][j] - rhspm[1][j]) / (2 * h);
            }
        }
        computeNewtonRaphsonTerms();
        LinearAlgebra.printMatrixDiff(jacobian, fdJacobian);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (Math.abs(fdJacobian[i][j] - jacobian[i][j]) > h) {
                    System.out.println("Problem at (i,j) = (" + i + ", " + j + "), Jacobian = " + jacobian[i][j] + " ~ " + fdJacobian[i][j]);
                }
            }
        }
    }

    /**
     * Computes the jacobian and the right hand side.
     */
    private double computeNewtonRaphsonTerms() {
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
        double residue = 0;
        for (int i = 0; i < rhs.length; i++) {
            if (Math.abs(rhs[i]) > residue) {
                residue = Math.abs(rhs[i]);
            }
        }
        return residue;
    }

    public double[][][] getSolution() {
        return vector.getVariables();
    }

    /**
     * Sets a guess. The guess is indexed in the following fashion:
     * guess[iDomain].get(iAbscissa)[iVariable].
     *
     * @param guess The guess to set.
     * @param parameters The guess at the parameter values.
     */
    public void setGuess(ArrayList<double[]>[] variables, double[] parameters) {
        vector.setGuess(variables, parameters);
    }
}
