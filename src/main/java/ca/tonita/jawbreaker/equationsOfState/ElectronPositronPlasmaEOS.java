package ca.tonita.jawbreaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.PhysicalQuantity;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;

/**
 * One may neglect positrons when eta >> -1/beta.  Which is basically
 * never inside of a neutron-star. Therefore, we need to account for 
 * the full degeneracy and the existence of positrons in general.<br>
 * <br>
 * The relations in this class can be found in Cox and Giuli's: 
 * Principles of Stellar Structure in section 24.9.
 * @author atonita
 *
 */
public class ElectronPositronPlasmaEOS {
	private UnitSystem myUnits = null;
	private UnitSystem units = null;
	private FermiDiracIntegral fd = null;
	private double k_b	= 1.0;	//1.3806504E-23;
	private double m_e	= 1.0;	//9.10938215E-31;
	private double c		= 1.0;	//2.99792458E8;
	private double h		= 1.0;	//6.62606896E-34;
	private int iterationMax = 75;
	private double desiredPrecision = 1.0E-8;
	
	/**
	 * Blank constructor. Use it.
	 */
	public ElectronPositronPlasmaEOS() {
		fd = new FermiDiracIntegral();
		myUnits =  new UnitSystem(SIConstants.c, SIConstants.melectron, SIConstants.h, SIConstants.k, SIConstants.COULOMB);
	}
	
	/**
	 * Gives the number density of the positrons as a function of the chemical potential, mu, and the temperature, T.
	 * This number density is defined as the number of electrons minus the number of positrons.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the number density of leptons
	 */
	public double positronNumberDensity(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double positronEta = -eta-2.0/beta;
		double numberDensity = fd.f(0.5, positronEta, beta) + beta*fd.f(1.5, positronEta, beta);		// Positron contribution
		numberDensity *= Math.pow(beta, 1.5);
		double frontFactor = 8.0*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.sqrt(2.0);	// Statistical weight of Cox & Giuli.
		numberDensity *= frontFactor;
		return numberDensity;
	}
	
	/**
	 * Returns the chemical potential as a function of lepton density and temperature.
	 * @param leptonNumberDensity the number of leptons per unit volume (electrons minus positrons) ie. the charge number density
	 * @param temperature the temperature of the plasma, T
	 * @return the chemical potential mu(n_l,T)
	 */
	public double chemicalPotential(double leptonNumberDensity, double temperature) {
		//double chemicalPotential = h*h/(2.0*m_e)*Math.pow(3.0*leptonNumberDensity/(8.0*Math.PI), 2.0/3.0); // Initial guess as the degenerate chemical potential.
		// The lepton density, chemical potential and temperature are in the input units now. As the functions require.
		double chemicalPotential = this.guessMu(leptonNumberDensity, temperature);
		return chemicalPotential(leptonNumberDensity,temperature,chemicalPotential);
	}
	
	/**
	 * Returns the chemical potential as a function of lepton density and temperature.<br>
	 * Checked 09/08/2010.
	 * @param leptonNumberDensity the number of leptons per unit volume (electrons minus positrons) ie. the charge number density
	 * @param temperature the temperature of the plasma, T
	 * @param guess the initial guess for the chemical potential
	 * @return the chemical potential mu(n_l,T)
	 */
	public double chemicalPotential(double leptonNumberDensity, double temperature, double guess) {
		double error = Double.MAX_VALUE;
		int iterations = 0;
		double chemicalPotential = guess;
		chemicalPotential = guessMu(leptonNumberDensity,temperature);
		while (Math.abs(error) > desiredPrecision && iterations < iterationMax) {
			double derivative = dNumberDensitydChemicalPotential(chemicalPotential,temperature);
			error = (leptonNumberDensity - numberDensity(chemicalPotential,temperature))/derivative;
			chemicalPotential += error;
			iterations++;
		}
		if (iterations == iterationMax) {
			System.err.println("Exceeded maximum number of iterations ("+iterationMax+").\n" +
					"leptonNumberDensity = "+leptonNumberDensity+" temperature = "+temperature+
					"\nLast error was "+error+" with chemical potential "+chemicalPotential+".");
			double beta = k_b*temperature/(m_e*c*c);
			System.err.println("beta = "+beta);
			System.err.flush();
		} 
		return chemicalPotential;
	}
	
	/**
	 * Returns the pressure of the plasma as a function of the chemical potential and temperature.<br>
	 * <br>
	 * Uses equations 24.363 and 24.367a of C&G.<br>
	 * Checked 09/08/2010.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the pressure of the plasma P(mu,T)
	 */
	public double pressure(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double pressure = fd.f(1.5, eta, beta) + 0.5*beta*fd.f(2.5, eta, beta);	// Electron contribution.
		double positronEta = -eta-2.0/beta;
		pressure += fd.f(1.5, positronEta, beta) + 0.5*beta*fd.f(2.5, positronEta, beta);
		pressure *= 16.0*Math.PI*Math.sqrt(2.)/(3.0*h*h*h)*Math.pow(m_e, 4.0)*Math.pow(c,5.0)*Math.pow(beta, 2.5);
		return pressure; // checked 27.10.2009
	}
	
	/**
	 * Returns the pressure of the plasma as a function of the chemical potential, the
	 * temperature is assumed to be zero.<br>
	 * <br> 
	 * Fermi Momentum comes from eq. 24.152 of Cox & Giuli.<br>
	 * Constant x comes from eq. 24.156.<br>
	 * The function f(x) is eq. 24.157.<br>
	 * The pressure is eq 24.158.
	 * @param chargeDensity the charge density (electron number density)
	 * @return the pressure of the plasma P(n_e)
	 */
	public double pressure(double chargeDensity) {
		double fermiMomentum = Math.pow(chargeDensity*3.*h*h*h/(8.*Math.PI), 1./3.);
		double x = fermiMomentum/(m_e*c);
		double fOfX = x*Math.sqrt(x*x + 1)*(2*x*x - 3) + 3.*Math.log(x + Math.sqrt(1 + x*x));
		double pressure = Math.PI*Math.pow(m_e*c/h, 3)*m_e*c*c*fOfX;
		return pressure;
	}

	/**
	 * Returns the pressure of the plasma as a function of the chemical potential, the
	 * temperature is assumed to be zero.<br>
	 * <br>
	 * Momentum comes from eq. 24.152 of Cox & Giuli.<br>
	 * Constant x comes from eq. 24.156.<br>
	 * The function g(x) is eq. 24.162.<br>
	 * The internal energy is formula 24.161
	 * @param chargeDensity the charge density (electron number density)
	 * @return the pressure of the plasma P(n_e)
	 */
	public double internalEnergyDensity(double chargeDensity) {
		double fermiMomentum = Math.pow(chargeDensity*3.*h*h*h/(8.*Math.PI), 1./3.);
		double x = fermiMomentum/(m_e*c);
		double fOfX = x*Math.sqrt(x*x + 1)*(2*x*x - 3) + 3.*Math.log(x + Math.sqrt(1 + x*x));
		double gOfX = 8*x*x*x*(Math.sqrt(1 + x*x) - 1) - fOfX;
		return Math.PI*Math.pow(m_e*c/h, 3)*m_e*c*c*gOfX;
	}

	/**
	 * Returns the internal energy density of the plasma as a function of the chemical potential and temperature.
	 * That is, the internal energy per unit volume.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the internal energy density of the plasma U(mu,T)/V
	 */
	public double internalEnergyDensity(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double internalEnergyDensity = fd.f(1.5, eta, beta) + beta*fd.f(2.5, eta, beta);
		double positronEta = -eta-2.0/beta;
		internalEnergyDensity += fd.f(1.5, positronEta, beta) + beta*fd.f(2.5, positronEta, beta);
		internalEnergyDensity *= 8.0*Math.PI*Math.sqrt(2.0)/(h*h*h)*Math.pow(m_e, 4)*Math.pow(c, 5.0)*Math.pow(beta, 2.5);
		// Need the contribution from the electron-positron pair rest mass.
		double n_pos = fd.f(0.5, positronEta, beta) + beta*fd.f(1.5, positronEta, beta);
		n_pos *= 8.0*Math.sqrt(2.0)*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.pow(beta, 1.5);
		internalEnergyDensity += 2.0*m_e*c*c*n_pos;
		return internalEnergyDensity; // checked 27.10.2009
	}
	
	/**
	 * Who likes using weird units? I do. Use this routine to tell me what units you want to give and receive physical quantities in.
	 * @param preferredUnits your preferred units
	 */
	public void setUnits(UnitSystem preferredUnits) {
		this.units = preferredUnits;
		PhysicalQuantity x = new PhysicalQuantity(1.0,Dimension.VELOCITY,myUnits); // c
		x = units.convert(x);
		this.c = x.getValue();
		x = new PhysicalQuantity(1.0,Dimension.MASS,myUnits); // m_e
		x = units.convert(x);
		this.m_e = x.getValue();
		x = new PhysicalQuantity(1.0,SIConstants.h.getDimension(),myUnits); // h
		x = units.convert(x);
		this.h = x.getValue();
		x = new PhysicalQuantity(1.0,SIConstants.k.getDimension(),myUnits); // k_boltzmann
		x = units.convert(x);
		this.k_b = x.getValue();
	}
	
	/**
	 * Shamefully stolen from Timmes equation of state, this function
	 * returns a pretty good guess of the chemical potential
	 * @param numberDensity the charge number density
	 * @param temperature the temperature
	 * @return a good guess of the chemical potential
	 */
	protected double guessMu(double numberDensity, double temperature) {
		// First we determine if the (numberDensity,temperature) point lies to the left of the turning point curve.
		// This would imply that positrons are important and we don't have a good estimate of eta. See Cox & Giuli sec 24.9.
		double beta = k_b*temperature/(m_e*c*c);
		double turningPointDensity;
		if (beta > 1.0) {
			turningPointDensity = beta*beta*4.0*8.0*Math.PI*Math.pow(m_e*c/h, 3.0)/3.0;	// C&G 24.345
		} else {
			turningPointDensity = beta*(1.0+0.75*beta)*Math.exp(-1.0/beta)*2.0*Math.sqrt(1.5*Math.PI)*8.0*Math.PI*Math.pow(m_e*c/h, 3.0)/3.0;
		}
		double eta;
		if (turningPointDensity > numberDensity) {
			// Positrons are important, most likely mu is between 0 and -1, so...
			eta = -0.5;
		} else {
			// We estimate eta through equation 24.309, so we compute the dimensionless number density
			double dimensionlessNumberDensitySquared = Math.pow(h/(m_e*c)*Math.pow(3.0*numberDensity/(8.0*Math.PI), 1.0/3.0),2.0);
			double firstGuess;
			if (dimensionlessNumberDensitySquared > 1.0e-6) {
				firstGuess = (Math.sqrt(1.0+dimensionlessNumberDensitySquared)-1.0)/beta;
			} else {
				// Root will give crap so we use a binomial expansion
				firstGuess = dimensionlessNumberDensitySquared*(1.0 - dimensionlessNumberDensitySquared/4.0)*0.5/beta;
			}
			// Now we run the divine approximation in reverse if the number density is large. Otherwise we use eta as is.
			//		Supposedly this is C&G 43, I copied it from Timmes.
			double divineNumberDensity = Math.log10(Math.pow(numberDensity, 0.6)/temperature);
			if (divineNumberDensity < 9.5) {
				double z = ((1.0+64.0/(9.0*Math.PI)*beta)*Math.sqrt(1.0+64.0/(9.0*Math.PI)*beta*0.5)-1.0)/(64.0/(9.0*Math.PI));
				eta = -Math.log(4.0*Math.PI*Math.pow(2.0*m_e*k_b*temperature/(h*h), 1.5)*(0.5+0.74*z)/numberDensity);
				if (eta > 8.5) {
					eta = eta*(9.5-divineNumberDensity)+firstGuess*(1.0-9.5+divineNumberDensity);
				}
			} else {
				eta = firstGuess;
			}
		}
		//return eta*k_b*temperature;
		return eta;
	}
	
	/**
	 * In case you want to know what units this routine wants to use.
	 * @return the unit system used by this subroutine.
	 */
	public UnitSystem getUnits() {
		return myUnits;
	}

	/**
	 * Gives the charge number density of the plasma as a function of the chemical potential, mu, and the temperature, T.
	 * This number density is defined as the number of electrons minus the number of positrons.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the number density of leptons
	 */
	public double numberDensity(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double numberDensity = fd.f(0.5, eta, beta) + beta*fd.f(1.5, eta, beta); // Electron contribution
		double positronEta = -eta-2.0/beta;
		numberDensity -= fd.f(0.5, positronEta, beta) + beta*fd.f(1.5, positronEta, beta);		// Positron contribution
		numberDensity *= Math.pow(beta, 1.5);
		double frontFactor = 8.0*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.sqrt(2.0);	// Statistical weight of Cox & Giuli.
		numberDensity *= frontFactor;
		return numberDensity;
	}
	
	/**
	 * The derivative of the number density with respect to the chemical potential.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return del/del mu (n_l(mu,T))
	 */
	public double dNumberDensitydChemicalPotential(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double derivativeNumberDensity = fd.dfdeta(0.5, eta, beta) + beta*fd.dfdeta(1.5, eta, beta);						// Electron contribution
		double positronEta = -eta-2.0/beta;
		derivativeNumberDensity += fd.dfdeta(0.5, positronEta, beta) + beta*fd.dfdeta(1.5, positronEta, beta);		// Positron contribution
		derivativeNumberDensity *= 8.0*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.sqrt(2.0)*Math.pow(beta, 1.5);	// Statistical weight of Cox & Giuli.
		derivativeNumberDensity *= 1.0/(k_b*temperature);	// Derivative of eta wrt chemicalPotential via the chain rule.
		return derivativeNumberDensity;
	}

	/**
	 * Gives the number density of electron-positron pairs as a function of the
	 * chemical potential, mu, and the temperature, T.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the number density of leptons
	 */
	public double pairNumberDensity(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double positronEta = -eta-2.0/beta;
		double numberDensity = fd.f(0.5, positronEta, beta) + beta*fd.f(1.5, positronEta, beta);		// Positron contribution
		numberDensity *= Math.pow(beta, 1.5);
		double frontFactor = 8.0*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.sqrt(2.0);	// Statistical weight of Cox & Giuli.
		numberDensity *= frontFactor;
		return numberDensity;
	}
	
	/**
	 * Returns the kinetic chemical potential parameter, which is one of the parameters of the fermi dirac integrals.
	 * @param chemicalPotential the chemical potential of the electrons
	 * @param temperature the temperature of the plasma
	 * @return the kinetic chemical potential mu/(k*T)
	 */
	public double eta(double chemicalPotential, double temperature) {
		return chemicalPotential/(k_b*temperature);
	}
	
	/**
	 * Returns the rest mass normalized thermal energy parameter, 
	 * which is one of the parameters of the fermi dirac integrals.
	 * @param chemicalPotential the chemical potential of the electrons
	 * @param temperature the temperature of the plasma
	 * @return the normalized thermal energy k*T/(m_e*c^2)
	 */
	public double beta(double chemicalPotential, double temperature) {
		return k_b*temperature/(m_e*c*c);
	}
	
	/**
	 * The entropy density.
	 * @param chemicalPotential the chemical potential of the lepton gas
	 * @param temperature the temperature of the lepton gas
	 * @return the entropy density S(mu,T)/V
	 */
	public double entropyDensity(double chemicalPotential, double temperature) {
		double u = internalEnergyDensity(chemicalPotential, temperature);
		double p = pressure(chemicalPotential, temperature);
		double n = numberDensity(chemicalPotential, temperature);
		double eta = eta(chemicalPotential, temperature);
		double entropyDensity = (u+p)/temperature - eta*k_b*n;
		return entropyDensity; //checked 27.10.2009
	}
	
	/**
	 * Returns the Gibb's free energy density as a function of the chemical
	 * potential and temperature. The free energy density is computed using
	 * the definition of the thermodynamic potential which follows from the
	 * Legendre transformation:<br>
	 * <br>
	 * <code> G = U - TS + PV</code>,<br>
	 * <br>
	 * yielding,<br>
	 * <br>
	 * <code>g = u - Ts + P</code>
	 * @param chemicalPotential the chemical potential of the lepton gas
	 * @param temperature the temperature of the lepton gas
	 * @return the Gibb's free energy density f = u - Ts + P
	 */
	public double gibbsFreeEnergyDensity(double chemicalPotential, double temperature) {
		double internalEnergyDensity = this.internalEnergyDensity(chemicalPotential, temperature);
		double entropyDensity = this.entropyDensity(chemicalPotential, temperature);
		double pressure = this.pressure(chemicalPotential, temperature);
		double freeEnergy = internalEnergyDensity - temperature*entropyDensity + pressure;
		return freeEnergy; // Checked 27.10.2009
	}

	/**
	 * Gives the number density for a zero temperature fermi gas of non-interacting electrons.
	 * @param chemicalPotential the chemcial potential of the electrons
	 * @return the number density of electrons
	 */
	public double numberDensity(double chemicalPotential) {
		return Math.pow(2.0*m_e*chemicalPotential/(h*h), 1.5)*8.0*Math.PI/3.0;
	}
	
	/**
	 * Returns the chemical potential for a zero temperature gas of electrons.<br>
	 * Uses eq. 24.152 to compute the fermi momentum, then the fermi energy, and thus the
	 * chemical potential is p_F^2/(2m)
	 * 
	 * @param electronNumberDensity the number density of electrons
	 * @return the chemcial potential at zero temperature
	 */
	public double chemicalPotential(double electronNumberDensity) {
		return h*h*0.5/m_e*Math.pow(3.0*electronNumberDensity/(8.0*Math.PI), 2.0/3.0);
	}

	/**
	 * Returns the Helmholtz's free energy density as a function of the chemical
	 * potential and temperature. The free energy density is computed using
	 * the definition of the thermodynamic potential which follows from the
	 * Legendre transformation:<br>
	 * <br>
	 * <code> F = U - TS</code>,<br>
	 * <br>
	 * yielding,<br>
	 * <br>
	 * <code>f = u - Ts</code>
	 * @param chemicalPotential the chemical potential of the lepton gas
	 * @param temperature the temperature of the lepton gas
	 * @return the Helmholtz's free energy density f = u - Ts
	 */
	public double helmholtzFreeEnergyDensity(double chemicalPotential, double temperature) {
		double internalEnergyDensity = this.internalEnergyDensity(chemicalPotential, temperature);
		double entropyDensity = this.entropyDensity(chemicalPotential, temperature);
		double freeEnergy = internalEnergyDensity - temperature*entropyDensity;
		return freeEnergy; // Checked 27.10.2009
	}
}
