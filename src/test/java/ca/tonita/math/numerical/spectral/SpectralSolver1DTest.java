package ca.tonita.math.numerical.spectral;

import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import ca.tonita.math.numerical.NonLinearFirstOrderODESystem;
import ca.tonita.math.numerical.ODEIndexer1D;
import ca.tonita.math.polynomials.ChebyshevExtrema;
import ca.tonita.math.polynomials.LinearlyMappedBasis;
import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class SpectralSolver1DTest extends TestCase implements NonLinearFirstOrderODESystem {
    private double tolerance = 1.0e-10;
    private int rank = 20;
    
    public SpectralSolver1DTest(String testName) {
        super(testName);
    }
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public NonLinearFirstOrderODEBean equations(int iX, double x, double[] y, double[] dy, double[] parameters, int type) {
        NonLinearFirstOrderODEBean bean = new NonLinearFirstOrderODEBean();
        double[] residue = new double[1];
        double[][] jacobian = new double[1][3];
        if (type == BULK || type == RIGHTBOUNDARY) {
        residue[0] = dy[0] - parameters[0]*y[0];
        jacobian[0][0] = -parameters[0];
        jacobian[0][1] = 1;
        jacobian[0][2] = -y[0];
        } else if (type == LEFTBOUNDARY) {
            residue[0] = y[0] - 1;
            jacobian[0][0] = 1;
        }
        bean.setResidue(residue);
        bean.setJacobian(jacobian);
        return bean;
    }

    public double computeConstraint(int iConstraint, ODEIndexer1D indexer, SpectralVector1D vector, double[] dconstraint) {
        for (int i = 0; i < dconstraint.length; i++) {
            dconstraint[i] = 0;
        }
        dconstraint[rank-1] = 1;
        return vector.getVariable(0, 0)[rank-1] - 2;
    }

    public int getNDomains() {
        return 1;
    }

    public int[] getNVariables() {
        return new int[]{1};
    }

    public int getNVariables(int iDomain) {
        return 1;
    }

    public int getNParameters() {
        return 1;
    }

    /**
     * Test of solve method, of class SpectralSolver1D.
     */
    public void testSolve() {
        System.out.println("solve");
        LinearlyMappedBasis[] bases = new LinearlyMappedBasis[]{new LinearlyMappedBasis(new ChebyshevExtrema())};
        bases[0].setRank(rank);
        bases[0].setDomain(new double[]{0, 1});
        SpectralSolver1D instance = new SpectralSolver1D(this, bases);
        double[] parameters = new double[]{Math.log(2) + 0.02*(Math.random() - 0.5)};
        instance.setGuess(getGuess(bases[0].getAbscissas()), parameters);
        instance.solve(tolerance, 5);
        double[][][] solution = instance.getSolution();
        double[] x = bases[0].getAbscissas();
        for (int i = 0; i < rank; i++) {
            assertEquals(Math.pow(2, x[i]), solution[0][0][i], 1.0E-6);
        }
    }
    
    private double[][][] getGuess(double[] x) {
        double[] guess = new double[x.length];
        for (int iX = 0; iX < x.length; iX++) {
            guess[iX] = Math.exp(x[iX]) + 0.02*(Math.random() - 0.5);
        }
        return new double[][][]{{guess}};
    }
}
