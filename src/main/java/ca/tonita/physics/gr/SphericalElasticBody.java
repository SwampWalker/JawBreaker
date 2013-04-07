package ca.tonita.physics.gr;

/**
 *
 * @author atonita
 */
public class SphericalElasticBody {
    protected static int in = 0;
    protected static int iepsilon = 1;
    protected static int iprr = 3;
    protected static int iptt = 4;
    /**
     * Returns an array holding the number density, energy per particle,
     * radial pressure and tangential pressure of the strained body.
     * @param xi the corresponding areal radial coordinate in the body manifold
     * @param r the areal coordinate radius in the spherically symmetry spacetime in which the body is embedded
     * @param m the mass potential of the spacetime in which the body is embedded.
     * @return 
     */
    public double[] getVariables(double xi, double r, double m) {
        throw new UnsupportedOperationException("Not yet implemented.");
    }
}
