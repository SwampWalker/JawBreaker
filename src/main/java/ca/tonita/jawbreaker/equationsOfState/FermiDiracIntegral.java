package ca.tonita.jawbreaker.equationsOfState;

public class FermiDiracIntegral {
	public enum IntegrandType {
		INTEGRAND, ETADERIVATIVE, BETADERIVATIVE
	}

	private AparicioParameters aparicio = new AparicioParameters();
	private LegendreQuadrature20 legendre = new LegendreQuadrature20();
	private LaguerreQuadrature20 laguerre = new LaguerreQuadrature20();
	
	
	/**
	 * The fermi dirac integral.
	 * @param k the order of the integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return F_k(eta,beta)
	 */
	public double f(double k, double eta, double beta) {
		double integral = (	integral0(k,eta,beta,IntegrandType.INTEGRAND) + integral1(k,eta,beta,IntegrandType.INTEGRAND) + 
				 			integral2(k,eta,beta,IntegrandType.INTEGRAND) + integral3(k,eta,beta,IntegrandType.INTEGRAND));
		return integral;
	}
	
	/**
	 * The derivative of the fermi dirac integral with respect to the degeneracy parameter.
	 * @param k the order of the integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return dF/deta_k(eta,beta)
	 */
	public double dfdeta(double k, double eta, double beta) {
		return 	(integral0(k,eta,beta,IntegrandType.ETADERIVATIVE) + integral1(k,eta,beta,IntegrandType.ETADERIVATIVE) + 
				 integral2(k,eta,beta,IntegrandType.ETADERIVATIVE) + integral3(k,eta,beta,IntegrandType.ETADERIVATIVE));
	}
	
	/**
	 * The derivative of the fermi dirac integral with respect to the relativity parameter.
	 * @param k the order of the integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return dF/dbeta_k(eta,beta)
	 */
	public double dfdbeta(double k, double eta, double beta) {
		return 	(integral0(k,eta,beta,IntegrandType.BETADERIVATIVE) + integral1(k,eta,beta,IntegrandType.BETADERIVATIVE) + 
				 integral2(k,eta,beta,IntegrandType.BETADERIVATIVE) + integral3(k,eta,beta,IntegrandType.BETADERIVATIVE));
	}

	/**
	 * The integral of the first integral.
	 * @param k the order of the fermi dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return the second integral
	 */
	public double integral0(double k, double eta, double beta, IntegrandType integrandType) {
		double s0 = aparicio.split0(eta);
		double a = 0.0;
		double b = Math.sqrt(s0);
		double jacobian = (b-a)/2.0;
		double integral = 0.0;
		for (int i = 0; i < 20; i++) {
			double x = x(a,b,i);
			switch (integrandType) {
			case INTEGRAND:
				integral += legendre.w[i]*2.0*integrand(x*x,k+0.5,eta,beta);
				break;
			case ETADERIVATIVE:
				integral += legendre.w[i]*2.0*derivativeIntegrandwrtEta(x*x,k+0.5,eta,beta);
				break;
			case BETADERIVATIVE:
				integral += legendre.w[i]*2.0*derivativeIntegrandwrtBeta(x*x,k+0.5,eta,beta);
				break;
			}
		}
		integral *= jacobian;
		return integral;
	}
	
	/**
	 * The integral of the second integral.
	 * @param k the order of the fermi dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return the second integral
	 */
	public double integral1(double k, double eta, double beta, IntegrandType integrandType) {
		double s0 = aparicio.split0(eta); double a = s0;
		double s1 = aparicio.split1(eta); double b = s1;
		double jacobian = (b-a)/2.0;
		double integral = 0.0;
		for (int i = 0; i < 20; i++) {
			double x = x(a,b,i);
			switch (integrandType) {
			case INTEGRAND:
				integral += legendre.w[i]*integrand(x,k,eta,beta);
				break;
			case ETADERIVATIVE:
				integral += legendre.w[i]*derivativeIntegrandwrtEta(x,k,eta,beta);
				break;
			case BETADERIVATIVE:
				integral += legendre.w[i]*derivativeIntegrandwrtBeta(x,k,eta,beta);
				break;
			}
		}
		integral *= jacobian;
		return integral;
	}
	
	/**
	 * The integral of the third integral.
	 * @param k the order of the fermi dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return the second integral
	 */
	public double integral2(double k, double eta, double beta, IntegrandType integrandType) {
		double a = aparicio.split1(eta);
		double b = aparicio.split2(eta);
		double jacobian = (b-a)/2.0;
		double integral = 0.0;
		for (int i = 0; i < 20; i++) {
			double x = x(a,b,i);
			switch (integrandType) {
			case INTEGRAND:
				integral += legendre.w[i]*integrand(x,k,eta,beta);
				break;
			case ETADERIVATIVE:
				integral += legendre.w[i]*derivativeIntegrandwrtEta(x,k,eta,beta);
				break;
			case BETADERIVATIVE:
				integral += legendre.w[i]*derivativeIntegrandwrtBeta(x,k,eta,beta);
				break;
			}
		}
		integral *= jacobian;
		return integral;
	}
	
	/**
	 * The integral of the fourth integral.
	 * @param k the order of the fermi dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta relativity parameter
	 * @return the second integral
	 */
	public double integral3(double k, double eta, double beta, IntegrandType integrandType) {
		double a = aparicio.split2(eta);
		double integral = 0.0;
		for (int i = 0; i < 20; i++) {
			double x = laguerre.x[i] + a;
			switch (integrandType) {
			case INTEGRAND:
				integral += laguerre.w[i]*integrand(x,k,eta,beta);
				break;
			case ETADERIVATIVE:
				integral += laguerre.w[i]*derivativeIntegrandwrtEta(x,k,eta,beta);
				break;
			case BETADERIVATIVE:
				integral += laguerre.w[i]*derivativeIntegrandwrtBeta(x,k,eta,beta);
				break;
			}
		}
		return integral;
	}
	
	/**
	 * This helper function takes roots of the Legendre polynomials and maps them onto the interval [a,b].
	 * @param a the left endpoint of the interval
	 * @param b the right endpoint of the interval
	 * @param i the index of the root
	 * @return the remapped root
	 */
	private double x(double a, double b, int i) {
		double x = (b-a)/2.0*(legendre.x[i]+1.0)+a;
		return x;
	}
	
	/**
	 * This is the integrand of the Fermi-Dirac integral.
	 * @param x the normalized energy epsilon/kT
	 * @param k the order of the Fermi-Dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta the relativity parameter
	 * @return the integrand of the Fermi-Dirac integral
	 */
	private double integrand(double x, double k, double eta, double beta) {
		double integrand = Math.pow(x, k)*Math.sqrt(1.0 + x*beta/2.0)/(Math.exp(x-eta) + 1.0);
		return integrand;
	}
	

	/**
	 * This is the derivative of the integrand of the Fermi-Dirac integral with respect to the relativity parameter.
	 * @param x the normalized energy epsilon/kT
	 * @param k the order of the Fermi-Dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta the relativity parameter
	 * @return the integrand of the Fermi-Dirac integral
	 */
	private double derivativeIntegrandwrtBeta(double x, double k, double eta, double beta) {
		double integrand = integrand(x,k+1,eta,beta)/(4.0+2.0*beta*x);
		return integrand;
	}
	

	/**
	 * This is the derivative of the integrand of the Fermi-Dirac integral with respect to the degeneracy parameter.
	 * @param x the normalized energy epsilon/kT
	 * @param k the order of the Fermi-Dirac integral
	 * @param eta the degeneracy parameter
	 * @param beta the relativity parameter
	 * @return the integrand of the Fermi-Dirac integral
	 */
	private double derivativeIntegrandwrtEta(double x, double k, double eta, double beta) {
		double integrand = integrand(x,k,eta,beta)/(1.0+Math.exp(eta-x));
		return integrand;
	}
}
