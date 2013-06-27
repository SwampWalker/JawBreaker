package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.EOSHandler;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.gr.hydro.TOVData;
import ca.tonita.math.numerical.EvolutionTerminator;
import ca.tonita.math.numerical.RK4;
import ca.tonita.physics.gr.hydro.TOVBuilder;
import ca.tonita.physics.gr.hydro.TOVEquations;
import ca.tonita.physics.gr.hydro.TOVIndex;
import ca.tonita.physics.gr.hydro.TOVTerminator;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Builds neutron stars with elastic crusts.
 *
 * @author atonita
 */
public class ElasticTOVBuilder {

    /**
     * Solves the elastic TOV equations using RK4 method.
     *
     * @param body the background TOV to use as the body manifold
     * @param eos the equation of state to use for the hydro core
     * @param stepSize the step size to use for the RK4
     * @param outputEvery how often to output data
     * @param minPressure the value of pressure to terminate evolution at (the
     * effective surface pressure)
     * @return the TOVData solution
     */
    public static TOVData createRK4(SphericalBodyManifoldRK4 body, TabulatedHermite eos, double stepSize, int outputEvery, double minPressure) {
        TOVEquations eqns = new TOVEquations(eos);
        // Create guess.
        ArrayList<double[]> points = new ArrayList<double[]>();
        ArrayList<Double> radii = new ArrayList<Double>();
        double centralPressure = ((TOVData)body.getBackground()).getPressure(0);
        double rCore = body.getBackground().getCoreRadius();
        double mCore = body.getBackground().getCoreMass();
        int maxSteps = Integer.MAX_VALUE;
        ElasticTOVCoreTerminator terminator = new ElasticTOVCoreTerminator(rCore*2, maxSteps, body.getBackground().getCoreRestMass(), minPressure);
        createGuessCore(points, centralPressure, radii, stepSize, eqns, outputEvery, terminator);
        System.out.println(body.getBackground().getCoreRadius() + " " + radii.get(radii.size() - 1) + " " + points.get(points.size() - 1)[TOVIndex.RESTMASS]);
        
        // Elastic equations.
        double r = radii.get(radii.size() - 1);
        double m = points.get(radii.size() - 1)[TOVIndex.MASS];
        double xiPrime = Math.sqrt((rCore*rCore - 2*mCore*rCore)/(r*r - 2*m*r));
        ElasticTOVEquations eqns2 = new ElasticTOVEquations(body);
        SphericalElasticBean quantities = body.getQuantities(rCore);
        double prr = eqns2.pressureRadial(r, quantities, rCore, xiPrime, m);
        System.out.println(prr + " " + points.get(points.size()-1)[TOVIndex.PRESSURE]);
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
        double centralPressure = 1e-4;
        double stepSize = 1.0e-3;
        int outputEvery = 1;
        double minPressure = 1.0e-8;
        TOVBuilder.evolve(background, eos2, centralPressure, stepSize, outputEvery, minPressure);
        SphericalBodyManifoldRK4 body = new SphericalBodyManifoldRK4(background, eos2);
        ElasticTOVBuilder.createRK4(body, eos3, stepSize, outputEvery, minPressure);
    }
}
