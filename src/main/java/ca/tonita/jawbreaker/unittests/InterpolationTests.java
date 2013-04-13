package ca.tonita.jawbreaker.unittests;

import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.shenzerotemperature.drivers.interpolators.Polynomials;

/**
 *
 * @author atonita
 */
public class InterpolationTests {
    
    public static void main(String[] args) {
        double[] x = {1, 2, 4, 8};
        double[] f = {f(x[0]), f(x[1]), f(x[2]), f(x[3])};
        
        // Test generic interpolation.
        double[] a = Polynomials.interpolatingCoefficients(f, x);
        double[] b = Polynomials.interpolatingCoefficients(new double[]{f(x[0]), f(x[1])}, new double[]{df(x[0]), df(x[1])}, new double[]{x[0], x[1]});
        System.out.println("Coefficients:");
        System.out.println("4 -1 2 1");
        System.out.println(a[0] + " " + a[1] + " " + a[2] + " " + a[3]);
        System.out.println(b[0] + " " + b[1] + " " + b[2] + " " + b[3]);
        System.out.println("Value and derivative:");
        System.out.println(f(16) + " = " + Polynomials.interpolate(a, 16));
        double h = 0.0001;
        System.out.println((f(16+h) - f(16-h))/(2*h) + " = " + Polynomials.differentiate(a, 16));
        
        // Test tabulated equation of state's ability to recover energy density.
        double logn0 = Math.log10(0.0005);
        double lognStep = Math.log10(1.1); 
        int nPoints = 15;
        double[] logn = new double[nPoints];
        double[] logp = new double[nPoints];
        double[] spE = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            logn[i] = logn0 + i*lognStep;
            logp[i] = Math.log10(polyP(Math.pow(10, logn[i])));
            spE[i] = polySpE(Math.pow(10, logn[i]));
        }
        TabulatedHermite eos = new TabulatedHermite(logn, logp, spE, 1);
        System.out.println("\nEnergy Density:");
        double n = Math.pow(10, logn[1]);
        double p = polyP(n);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        n = Math.pow(10, logn[1]) + 0.01*Math.pow(10, logn[1]);
        p = polyP(n);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        n = Math.pow(10, logn[1]) - 0.01*Math.pow(10, logn[1]);
        p = polyP(n);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        n = (Math.pow(10, logn[1]) + Math.pow(10, logn[2]))*0.5;
        p = polyP(n);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        
        // Test extrapolations.
        n = Math.pow(10, logn[0]);
        p = Math.pow(10, logp[0]);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        n = n - 0.01*n;
        p = polyP(n);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        n = Math.pow(10, logn[nPoints-1]);
        p = Math.pow(10, logp[nPoints-1]);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
        n = n + 0.01*n;
        p = polyP(n);
        System.out.println("n=" + n +", rho = " + polyRho(p) + " ~ " + eos.energyDensity(p) + ", drho = " + dpolyRho(p) + " ~ " + eos.denergyDensity(p));
    }
    
    private static double f(double x) {
        return 4 + x*(-1 + x*(2 + x)); // x^3 + 2x^2 - x + 4
    }
    
    private static double df(double x) {
        return -1 + x*(4 + 3*x);  // 3x^2 + 4x - 1
    }
    
    private static double polyP(double n) {
        return 100*n*n;
    }
    
    private static double polySpE(double n) {
        return 100*n;
    }
    
    private static double polyRho(double p) {
        return p + 0.1*Math.sqrt(p);
    }
    
    private static double dpolyRho(double p) {
        return 1 + 0.05/Math.sqrt(p);
    }
}
