package ca.tonita.physics.gr.elastic;

import ca.tonita.jawbreaker.equationsOfState.EOSHandler;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.TOVFamily;
import ca.tonita.physics.gr.hydro.TOVData;
import ca.tonita.math.numerical.EvolutionTerminator;
import ca.tonita.math.numerical.RK4;
import ca.tonita.physics.gr.hydro.TOVBuilder;
import ca.tonita.physics.gr.hydro.TOVEquations;
import ca.tonita.physics.gr.hydro.TOVIndex;
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
        // Bound the rest mass.
        double mRestCore = body.getBackground().getCoreRestMass();
        int rightBoundary = eosFamily.getIMaximum();
        int leftBoundary = 0;
        for (int i = 0; i < rightBoundary && leftBoundary == 0; i++) {
            if (eosFamily.get(i).getRestMass() >= mRestCore) {
                leftBoundary = i - 1;
            } else if (i + 1 == rightBoundary) {
                leftBoundary = rightBoundary;
            }
        }
        if (leftBoundary < 0) {
            // TODO: create reference at something like 2*minPressure and shoot between minimum and this.
            throw new UnsupportedOperationException("Cannot find solution if TOVFamily's minimum rest mass is greater than the rest mass of the background solution.");
        } else if (leftBoundary == rightBoundary) {
            // TODO: black hole + crust.
            throw new UnsupportedOperationException("Black hole core not implemented.");
        }
        for (int i = leftBoundary; i <= rightBoundary; i++) {
            System.out.println(mRestCore + " " + eosFamily.get(i).getRestMass());
        }
        System.out.println(leftBoundary + " " + rightBoundary);

        // Find the point in the solutions where the elastic pressure becomes greater than the hydro pressure.
        ElasticTOVEquations elasticEqns = new ElasticTOVEquations(body);
        // TODO Capture numeric parameters in a single bean and use it here instead of minPressure.
        double r = eosFamily.getMaximum().getRadiusByRestMass(mRestCore, minPressure);
        double[] variables = eosFamily.getMaximum().getVariables(r);
        double xiPrime = elasticEqns.xiPrimeCoreIsotropic(r, variables[TOVIndex.MASS]);
        SphericalElasticBean coreQuantities = body.getCoreQuantities();
        double prr = elasticEqns.pressureRadial(r, coreQuantities, body.getBackground().getCoreRadius(), xiPrime, variables[TOVIndex.MASS]);
        int sign = 1;
        if (prr > variables[TOVIndex.PRESSURE]) {
            // TODO: TOV + Vacuum + crust solution.
            throw new UnsupportedOperationException("Triple domain solution not supported.");
        }
        System.out.println(prr + " " + variables[TOVIndex.PRESSURE]);
        int iLeft = rightBoundary;
        while (prr < variables[TOVIndex.PRESSURE] && iLeft >= leftBoundary) {
            iLeft--;
            TOVData star = eosFamily.get(iLeft);
            r = star.getRadiusByRestMass(mRestCore, minPressure);
            variables = star.getVariables(r);
            xiPrime = elasticEqns.xiPrimeCoreIsotropic(r, variables[TOVIndex.MASS]);
            coreQuantities = body.getCoreQuantities();
            prr = elasticEqns.pressureRadial(r, coreQuantities, body.getBackground().getCoreRadius(), xiPrime, variables[TOVIndex.MASS]);
            System.out.println(r + " " + prr + " " + variables[TOVIndex.PRESSURE]);
        }
        if (iLeft < leftBoundary) {
            // We should not get here...
            throw new RuntimeException("This case should have been caught above: triple domain solution detected.");
        }

        // We now need to bisect between family.get(iLeft) and family.get(iLeft+1)
        //      until prr = p
        TOVData left = eosFamily.get(iLeft);
        TOVData right = eosFamily.get(iLeft + 1);
        TOVData middle; // For storing the final solution.
        int iterations = 0;
        while (Math.abs(prr - variables[TOVIndex.PRESSURE]) > 1.0e-8*prr && iterations < 30) { // TODO: parameterise these numeric parameters.
            iterations++;
            middle = new TOVData(eos);
            TOVBuilder.evolve(middle, eos, 0.5 * (left.getPressure(0) + right.getPressure(0)), stepSize, outputEvery, minPressure);
            r = middle.getRadiusByRestMass(mRestCore, minPressure);
            variables = middle.getVariables(r);
            xiPrime = elasticEqns.xiPrimeCoreIsotropic(r, variables[TOVIndex.MASS]);
            coreQuantities = body.getCoreQuantities();
            prr = elasticEqns.pressureRadial(r, coreQuantities, body.getBackground().getCoreRadius(), xiPrime, variables[TOVIndex.MASS]);
            System.out.println(iterations + " " + r + " " + prr + " " + variables[TOVIndex.PRESSURE] + " " + xiPrime);
            if (prr >= variables[TOVIndex.PRESSURE]) {
                left = middle;
            } else {
                right = middle;
            }
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
        TOVData background = new TOVData(eos2);
        double centralPressure = 1.1e-4;
        double stepSize = 1.0e-3;
        int outputEvery = 1;
        double minPressure = 1.0e-8;
        TOVBuilder.evolve(background, eos2, centralPressure, stepSize, outputEvery, minPressure);
        SphericalBodyManifoldRK4 body = new SphericalBodyManifoldRK4(background, eos2);
        TOVFamily family = new TOVFamily(eos3);
        int nTOVs = 100;
        double minFamilyPressure = 1.0e-5;
        double maxFamilyPressure = 0.001;
        TOVBuilder.fillFamily(family, nTOVs, minFamilyPressure, maxFamilyPressure, TOVBuilder.QUADRATICLY_SPACED, stepSize, outputEvery, minPressure, true);
        ElasticTOVBuilder.createRK4(body, eos3, stepSize, outputEvery, minPressure, family);
    }
}
