package ca.tonita.jawbreaker.equationsOfState;

public class AparicioParameters {
	public static double [] a = {6.7774,		3.7601,			7.5669};
	public static double [] b = {1.1418,		9.3719E-2,		1.1695};
	public static double [] c = {2.9826,		2.1063E-2,		7.5416E-1};
	public static double [] d = {0.0000,		3.1084E+1,		6.6558};
	public static double [] e = {0.0000,		1.0056,			-1.2819E-1};
	public static double D = 3.3609;
	public static double sigma = 9.1186E-2;
	
	/**
	 * This is the second split, counting up from 0.
	 * @param eta the degeneracy parameter, eta
	 * @return the location of the second domain's endpoint, the third domain's starting point
	 */
	public double split1(double eta) {
		double xi = xi(eta);
		double split = (a[0] + b[0]*xi + c[0]*xi*xi) / (1.0 + c[0]*xi);
		return split;
	}	

	/**
	 * This is the first split, counting up from 0.
	 * @param eta the degeneracy parameter, eta
	 * @return the location of the first domain's endpoint, the second domain's starting point
	 */
	public double split0(double eta) {
		double xi = xi(eta);
		double split = (a[1] + b[1]*xi + c[1]*d[1]*xi*xi) / (1.0 + e[1]*xi + c[1]*xi*xi);
		double split0 = split1(eta)-split;
		return split0;
	}
	
	/**
	 * This is the third split, counting up from 0.
	 * @param eta the degeneracy parameter, eta
	 * @return the location of the third domain's endpoint, the fourth domain's starting point
	 */
	public double split2(double eta) {
		double xi = xi(eta);
		double split = (a[2] + b[2]*xi + c[2]*d[2]*xi*xi) / (1.0 + e[2]*xi + c[2]*xi*xi);
		return split1(eta)+split;
	}
	
	/**
	 * Returns a normalized quantity which maps the original domain of sigma (-infinity,infinity) onto (0,infinity).
	 * @param eta the degeneracy parameter
	 * @return the normalize degeneracy
	 */
	public double xi(double eta) {
		double exponent = sigma*(eta-D);
		double xi = 0.0;
		if (exponent < 20) {
			xi = Math.log(1. + Math.exp(sigma*(eta-D)))/sigma;
		} else {
			xi = exponent/sigma;
		}
		return xi;
	}
}
