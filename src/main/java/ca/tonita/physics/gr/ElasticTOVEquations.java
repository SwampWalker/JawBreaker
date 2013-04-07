package ca.tonita.physics.gr;

/**
 *
 * @author atonita
 */
public class ElasticTOVEquations {
    // Indexes

    protected static int iprr = 0;
    protected static int im = 1;
    protected static int ilambda = 2;
    protected static int ixi = 3;
    protected static int nVariables = 4;
    private static int RESIDUE = 0;
    private static int JACOBIAN = 1;
    private SphericalElasticBody body;

    private Object compute(double[] y, double[] dy, double r, int quantities) {
        double[] residue = new double[nVariables];
        double[] bodyVariables = body.getVariables(y[ixi], r, y[im]);
        double rho = bodyVariables[body.in] * bodyVariables[body.iepsilon];
        residue[im] = dy[im] - 4 * Math.PI * rho * r * r;
        double dlambda = (y[im] + 4 * Math.PI * Math.pow(r, 3) * y[iprr])
                / (r * (r - 2 * y[im]));
        residue[ilambda] = dy[ilambda] - dlambda;
        residue[iprr] = dy[iprr] + dlambda * (y[iprr] + rho) + 2 / r * (y[iprr] - bodyVariables[body.iptt]);
        residue[ixi] = y[iprr] - bodyVariables[body.iprr];
        if (quantities == RESIDUE) {
            return residue;
        }
        throw new UnsupportedOperationException("Not yet implemented");
    }
    
    public double[] residue(double[] y, double[] dy, double r) {
        return (double[])compute(y, dy, r, RESIDUE);
    }
}
