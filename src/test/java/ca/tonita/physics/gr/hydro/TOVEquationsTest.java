package ca.tonita.physics.gr.hydro;

import ca.tonita.jawbreaker.equationsOfState.Polytrope;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.math.numerical.NonLinearFirstOrderODEBean;
import junit.framework.TestCase;

/**
 *
 * @author atonita
 */
public class TOVEquationsTest extends TestCase {

    public TOVEquationsTest(String testName) {
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

    /**
     * Test of equations method, of class TOVEquations.
     */
    public void testEquations() {
        System.out.println("equations");
        int iR = 1;
        double r = Math.random();
        double[] y = new double[]{Math.random(), Math.random(), Math.random()};
        double[] dy = new double[]{Math.random(), Math.random(), Math.random()};
        double[] parameters = new double[]{Math.random()};
        Polytrope poly = new Polytrope(100, 2, 1);
        int nPoints = 50;
        double[] logn = new double[nPoints];
        double[] logp = new double[nPoints];
        double[] energyPerParticle = new double[nPoints];
        double[] A = new double[nPoints];
        double[] Z = new double[nPoints];
        double logNMin = 1.0E-3;
        double logNMax = 1;
        for (int i = 0; i < nPoints; i++) {
            logn[i] = logNMin + i * (logNMax - logNMin) / (nPoints - 1);
            double n = Math.pow(10, logn[i]);
            logp[i] = Math.log10(poly.pressure(n));
            energyPerParticle[i] = poly.energyPerParticle(n);
            A[i] = 56;
            Z[i] = 26;
        }
        TabulatedHermite eos = new TabulatedHermite(logn, logp, energyPerParticle, 1, A, Z);
        TOVEquations instance = new TOVEquations(eos);
        for (int iType = 0; iType < 3; iType++) {
            int type = iType;
            NonLinearFirstOrderODEBean result = instance.equations(iR, r, y, dy, parameters, type);
            double[][] jacobian = new double[3][7];
            double[] yp = new double[3];
            double[] ym = new double[3];
            double[] dyp = new double[3];
            double[] dym = new double[3];
            double h = 1.0E-6;
            for (int iVariable = 0; iVariable < 3; iVariable++) {
                System.arraycopy(y, 0, yp, 0, 3);
                System.arraycopy(y, 0, ym, 0, 3);
                System.arraycopy(dy, 0, dyp, 0, 3);
                System.arraycopy(dy, 0, dym, 0, 3);
                yp[iVariable] = yp[iVariable] + h;
                ym[iVariable] = ym[iVariable] - h;
                double[] residueP = instance.equations(iR, r, yp, dy, parameters, type).getResidue();
                double[] residueM = instance.equations(iR, r, ym, dy, parameters, type).getResidue();
                for (int iEquation = 0; iEquation < 3; iEquation++) {
                    jacobian[iEquation][iVariable] = (residueP[iEquation] - residueM[iEquation]) / (2 * h);
                }
                dyp[iVariable] = dyp[iVariable] + h;
                dym[iVariable] = dym[iVariable] - h;
                residueP = instance.equations(iR, r, y, dyp, parameters, type).getResidue();
                residueM = instance.equations(iR, r, y, dym, parameters, type).getResidue();
                for (int iEquation = 0; iEquation < 3; iEquation++) {
                    jacobian[iEquation][3 + iVariable] = (residueP[iEquation] - residueM[iEquation]) / (2 * h);
                }
            }
            // Parameter derivatives.
            double rp = r / parameters[0] * (parameters[0] + h);
            double rm = r / parameters[0] * (parameters[0] - h);
            System.arraycopy(dy, 0, dyp, 0, 3);
            System.arraycopy(dy, 0, dym, 0, 3);
            for (int iVariable = 0; iVariable < 3; iVariable++) {
                dyp[iVariable] *= parameters[0] / (parameters[0] + h);
                dym[iVariable] *= parameters[0] / (parameters[0] - h);
            }
            double[] residueP = instance.equations(iR, rp, y, dyp, parameters, type).getResidue();
            double[] residueM = instance.equations(iR, rm, y, dym, parameters, type).getResidue();
            for (int iEquation = 0; iEquation < 3; iEquation++) {
                jacobian[iEquation][6] = (residueP[iEquation] - residueM[iEquation]) / (2 * h);
            }

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    assertEquals(result.getJacobian()[i][j], jacobian[i][j], 0.001);
                }
            }
        }
    }
}
