package ca.tonita.math.numerical;

import ca.tonita.math.polynomials.PolynomialBasis;

/**
 *
 * @author atonita
 */
public class SpectralSolver1D {
    private PolynomialBasis basis;
    private NonLinearFirstOrderODESystem system;
    private double[][] variables;
    private double[][] dvariables;
    private double[] parameters;
    private double[][] jacobian;
    private double[] rhs;
    
    public SpectralSolver1D() {
    }
    
    public SpectralSolver1D(PolynomialBasis basis, NonLinearFirstOrderODESystem system) {
        this.basis = basis;
        this.system = system;
    }
    
    private void computeNewtonRaphsonTerms() {
        double[] x = basis.getAbscissas();
        double[][] diff = basis.getDifferentiationMatrix();
        int nAbscissas = x.length;
        int nVariables = variables[0].length;
        // initialise
        for (int iRow = 0; iRow < jacobian.length; iRow++) {
            for (int iCol = 0; iCol < jacobian[0].length; iCol++) {
                jacobian[iRow][iCol] = 0;
            }
        }
        
        // Compute
        for (int iX = 0; iX < nAbscissas; iX++) {
            int type = NonLinearFirstOrderODESystem.BULK;
            if (iX == 0) {
                type = NonLinearFirstOrderODESystem.LEFTBOUNDARY;
            } else if (iX == nAbscissas-1) {
                type = NonLinearFirstOrderODESystem.RIGHTBOUNDARY;
            }
            NonLinearFirstOrderODEBean bean = system.equations(x[iX], variables[iX], dvariables[iX], parameters, type);
            double[] residue = bean.getResidue();
            double[][] dresidue = bean.getJacobian();
            for (int iEquation = 0; iEquation < nVariables; iEquation++) {
                int iRow = index(iEquation, iX);
                rhs[iRow] = residue[iEquation];
                // Now jacobian.
                for (int iVariable = 0; iVariable < nVariables; iVariable++) {
                    // Direct part.
                    int iCol = index(iVariable, iX);
                    jacobian[iRow][iCol] = dresidue[iEquation][iVariable];
                    // Derivative part.
                    for (int i = 0; i < nVariables; i++) {
                        iCol = index(iVariable, i);
                        jacobian[iRow][iCol] += diff[iX][i]*dresidue[iEquation][nVariables + iVariable];
                    }
                }
                // Derivatives w.r.t. parameters.
                for (int iParameter = 0; iParameter < parameters.length; iParameter++) {
                    jacobian[iRow][index(iParameter)] = dresidue[iEquation][nVariables*2 + iParameter];
                }
            }
        }
        // constraints.
        for (int iParameter = 0; iParameter < parameters.length; iParameter++) {
            system.computeConstraints(x, variables, dvariables, parameters, jacobian[index(iParameter)], iParameter);
        }
    }

    public NonLinearFirstOrderODESystem getSystem() {
        return system;
    }

    public void setSystem(NonLinearFirstOrderODESystem system) {
        this.system = system;
    }

    public PolynomialBasis getBasis() {
        return basis;
    }

    public void setBasis(PolynomialBasis basis) {
        this.basis = basis;
    }
    
    /**
     * Returns the index into the total variable vector.
     * @param iVariable The index of the variable.
     * @param iAbscissa The index of the abscissa.
     * @return the index
     */
    private int index(int iVariable, int iAbscissa) {
        return variables[0].length*iVariable + iAbscissa;
    }
    
    private int index(int iParameter) {
        return variables[0].length*variables.length + iParameter;
    }
}
