package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.EOSHandler;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.TOVFamily;
import ca.tonita.physics.gr.hydro.TOVData;
import ca.tonita.math.numerical.EvolutionTerminator;
import ca.tonita.math.numerical.RK4;
import ca.tonita.physics.gr.hydro.TOVBuilder;
import ca.tonita.physics.gr.hydro.TOVEquations;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Builds neutron stars with elastic crusts.
 *
 * @author atonita
 */
public final class ElasticTOVBuilder {

    private static final int UNKNOWN = 0;
    private static final int STABLE = 1;
    private static final int UNSTABLE = 2;

    /**
     * Solves the elastic TOV equations using RK4 method.
     *
     * @param body the background TOV to use as the body manifold
     * @param eos the equation of state to use for the hydro core
     * @param stepSize the step size to use for the RK4
     * @param outputEvery how often to output data
     * @param minPressure the value of pressure to terminate evolution at (the
     * effective surface pressure)
     * @param eosFamily a TOVFamily for the equation of state to help with the
     * shooting
     * @return the TOVData solution
     */
    public static TOVData createRK4(SphericalBodyManifoldRK4 body, TabulatedHermite eos, double stepSize, int outputEvery, double minPressure, TOVFamily eosFamily) {
        TOVEquations eqns = new TOVEquations(eos);
        // Create guess.
        ArrayList<double[]> points = new ArrayList<double[]>();
        ArrayList<Double> radii = new ArrayList<Double>();
        double centralPressure = ((TOVData) body.getBackground()).getPressure(0);
        double rCore = body.getBackground().getCoreRadius();
        double mCore = body.getBackground().getCoreMass();
        double mRestCore = body.getBackground().getCoreRestMass();
        // Bound the rest mass.
        int rightBoundary = eosFamily.getIMaximum();
        int leftBoundary = 0;
        for (int i = 0; i < rightBoundary && leftBoundary == 0; i++) {
            System.out.println(eosFamily.get(i).getConservedMass() + " " + mRestCore);
            if (eosFamily.get(i).getConservedMass() >= mRestCore) {
                leftBoundary = i - 1;
            } else if (i + 1 == rightBoundary) {
                leftBoundary = rightBoundary;
            }
        }
        if (leftBoundary < 0) {
            // TODO: create reference at something like 2*minPressure and shoot between minimum and this.
            throw new UnsupportedOperationException("Cannot find solution if TOVFamily's minimum rest mass is greater than the rest mass of the background solution.");
        } else if (leftBoundary == rightBoundary) {
            throw new UnsupportedOperationException("Black hole core not implemented.");
        }
        System.out.println(leftBoundary + " " + rightBoundary);

        // Start checking xi'.
        for (int i = leftBoundary + 1; i < rightBoundary; i++) {
            
        }
        throw new UnsupportedOperationException("Not finished yet.");
    }

    private static void createGuessCore(ArrayList<double[]> points, double centralPressure, ArrayList<Double> radii, double stepSize, TOVEquations eqns, int outputEvery, EvolutionTerminator terminator) {
        points.add(new double[]{centralPressure, 0, 0, 0});
        radii.add(0.);
        RK4.evolve(points, radii, eqns, stepSize, outputEvery, terminator);
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {
        TabulatedHermite eos2 = EOSHandler.readTableFromFile(new File("F:\\Projects\\JawBreaker\\data\\eos2.shear.betaEquilibrium.2.t00"));
        TabulatedHermite eos3 = EOSHandler.readTableFromFile(new File("F:\\Projects\\JawBreaker\\data\\eos3.shear.betaEquilibrium.2.t00"));
        TOVData background = new TOVData();
        double centralPressure = 1.5e-4;
        double stepSize = 1.0e-3;
        int outputEvery = 1;
        double minPressure = 1.0e-8;
        TOVBuilder.evolve(background, eos2, centralPressure, stepSize, outputEvery, minPressure);
        SphericalBodyManifoldRK4 body = new SphericalBodyManifoldRK4(background, eos2);
        TOVFamily family = new TOVFamily(eos3);
        int nTOVs = 100;
        double minFamilyPressure = 1.0e-4;
        double maxFamilyPressure = 0.001;
        TOVBuilder.fillFamily(family, nTOVs, minFamilyPressure, maxFamilyPressure, TOVBuilder.QUADRATICLY_SPACED, stepSize, outputEvery, minPressure);
        TOVBuilder.findMaxima(family, stepSize, outputEvery, minPressure, 1.0e-5);
        ElasticTOVBuilder.createRK4(body, eos3, stepSize, outputEvery, minPressure, family);
    }
}
