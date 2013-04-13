package ca.tonita.jawbreaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.PhysicalQuantity;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;

/**
 * One may neglect positrons when eta >> -1/beta.
 * @author atonita
 *
 */
public class ElectronPlasmaEOS {
	private UnitSystem myUnits = null;
	private UnitSystem inputUnits = null;
	private FermiDiracIntegral fd = null;
	private double k_b	= 1.0;	//1.3806504E-23;
	private double m_e	= 1.0;	//9.10938215E-31;
	private double c		= 1.0;	//2.99792458E8;
	private double h		= 1.0;	//6.62606896E-34;
	private double effectiveZeroTemperature = 1.0E-16; // About 6.0E-7 Kelvin. Pretty damn cold.
	
	/**
	 * Blank constructor. Use it.
	 */
	public ElectronPlasmaEOS() {
		fd = new FermiDiracIntegral();
		myUnits =  new UnitSystem(SIConstants.c, SIConstants.melectron, SIConstants.h, SIConstants.k, SIConstants.COULOMB);
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
	 * The derivative of the number density with respect to the chemical potential.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return del/del mu (n_l(mu,T))
	 */
	public double dNumberDensitydChemicalPotential(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double derivativeNumberDensity = fd.dfdeta(0.5, eta, beta) + fd.dfdeta(1.5, eta, beta);						// Electron contribution
		derivativeNumberDensity *= 8.0*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.sqrt(2.0)*Math.pow(beta, 1.5);	// Statistical weight of Cox & Giuli.
		derivativeNumberDensity *= 1.0/(k_b*temperature);	// Derivative of eta wrt chemicalPotential via the chain rule.
		return derivativeNumberDensity;
	}
	
	/**
	 * Returns the chemical potential as a function of lepton density and temperature.
	 * @param leptonNumberDensity the number of leptons per unit volume (electrons minus positrons) ie. the charge number density
	 * @param temperature the temperature of the plasma, T
	 * @return the chemical potential mu(n_l,T)
	 */
	public double chemicalPotential(double leptonNumberDensity, double temperature) {
		if (temperature < this.effectiveZeroTemperature) return chemicalPotential(leptonNumberDensity);
		//double chemicalPotential = h*h/(2.0*m_e)*Math.pow(3.0*leptonNumberDensity/(8.0*Math.PI), 2.0/3.0); // Initial guess as the degenerate chemical potential.
		// The lepton density, chemical potential and temperature are in the input units now. As the functions require.
		double chemicalPotential = this.guessMu(leptonNumberDensity, temperature);
		double error = Double.MAX_VALUE;
		int iterations = 0;
		while (Math.abs(error/chemicalPotential) > 1.E-8 && iterations < 5000) {
			error = (leptonNumberDensity - numberDensity(chemicalPotential,temperature))/dNumberDensitydChemicalPotential(chemicalPotential,temperature);
			chemicalPotential += error;
			iterations++;
		}
		return chemicalPotential;
	}
	
	/**
	 * Returns the chemical potential for a zero temperature gas of electrons.
	 * @param electronNumberDensity the number density of electrons
	 * @return the chemcial potential at zero temperature
	 */
	public double chemicalPotential(double electronNumberDensity) {
		return h*h/(2.0*m_e)*Math.pow(3.0*electronNumberDensity/(8.0*Math.PI), 2.0/3.0);
	}
	
	/**
	 * Returns the pressure of the plasma as a function of the chemical potential and temperature.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the pressure of the plasma P(mu,T)
	 */
	public double pressure(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double pressure = fd.f(1.5, eta, beta) + 0.5*beta*fd.f(2.5, eta, beta);	// Electron contribution.
		pressure *= 16.0*Math.PI*Math.sqrt(2.)/(3.0*h*h*h)*Math.pow(m_e, 4.0)*Math.pow(c,5.0)*Math.pow(beta, 2.5);
		return pressure;
	}
	
	/**
	 * Returns the internal energy density of the plasma as a function of the chemical potential and temperature.
	 * That is, the internal energy per unit volume.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the pressure of the plasma P(mu,T)
	 */
	public double internalEnergyDensity(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double internalEnergyDensity = fd.f(1.5, eta, beta) + beta*fd.f(2.5, eta, beta);
		internalEnergyDensity *= 8.0*Math.PI*Math.sqrt(2.0)/(h*h*h)*Math.pow(m_e, 4)*Math.pow(c, 5.0)*Math.pow(beta, 2.5);
		return internalEnergyDensity;
	}
	
	/**
	 * Who likes using weird units? I do. Use this routine to tell me what units you want to give and receive physical quantities in.
	 * @param preferredUnits your preferred units
	 */
	public void setInputUnits(UnitSystem preferredUnits) {
		this.inputUnits = preferredUnits;
		PhysicalQuantity x = new PhysicalQuantity(1.0,Dimension.VELOCITY,myUnits); // c
		x = inputUnits.convert(x);
		this.c = x.getValue();
		x = new PhysicalQuantity(1.0,Dimension.MASS,myUnits); // m_e
		x = inputUnits.convert(x);
		this.m_e = x.getValue();
		x = new PhysicalQuantity(1.0,SIConstants.h.getDimension(),myUnits); // h
		x = inputUnits.convert(x);
		this.h = x.getValue();
		x = new PhysicalQuantity(1.0,SIConstants.k.getDimension(),myUnits); // k_boltzmann
		x = inputUnits.convert(x);
		this.k_b = x.getValue();
		x = new PhysicalQuantity(this.effectiveZeroTemperature,Dimension.TEMPERATURE,myUnits);
		x = inputUnits.convert(x);
		this.effectiveZeroTemperature = x.getValue(); 
	}
	
	public double guessMu(double numberDensity, double temperature) {
		// First we determine if the (numberDensity,temperature) point lies to the left of the turning point curve.
		// This would imply that positrons are important and we don't have a good estimate of eta. See Cox & Giuli sec 24.9.
		double beta = k_b*temperature/(m_e*c*c);
		double turningPointDensity = 0.0;
		if (beta > 1.0) {
			turningPointDensity = beta*beta*4.0*8.0*Math.PI*Math.pow(m_e*c/h, 3.0)/3.0;	// C&G 24.345
		} else {
			turningPointDensity = beta*(1.0+0.75*beta)*Math.exp(-1.0/beta)*2.0*Math.sqrt(1.5*Math.PI)*8.0*Math.PI*Math.pow(m_e*c/h, 3.0)/3.0;
		}
		double eta = 0.0;
		if (turningPointDensity > numberDensity) {
			// Positrons are important, most likely mu is between 0 and -1, so...
			eta = -0.5;
		} else {
			// We estimate eta through equation 24.309, so we compute the dimensionless number density
			double dimensionlessNumberDensitySquared = Math.pow(h/(m_e*c)*Math.pow(3.0*numberDensity/(8.0*Math.PI), 1.0/3.0),2.0);
			double firstGuess = 0.0;
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
		return eta*k_b*temperature;
	}
	
	/**
	 * In case you want to know what units this routine wants to use.
	 * @return the unit system used by this subroutine.
	 */
	public UnitSystem getUnits() {
		return myUnits;
	}

	/**
	 * Gives the number density of the plasma as a function of the chemical potential, mu, and the temperature, T.
	 * This number density is defined as the number of electrons minus the number of positrons.
	 * @param chemicalPotential the chemical potential mu
	 * @param temperature the temperature T
	 * @return the number density of leptons
	 */
	public double numberDensity(double chemicalPotential, double temperature) {
		double eta = chemicalPotential/(k_b*temperature);
		double beta = k_b*temperature/(m_e*c*c);
		double numberDensity = fd.f(0.5, eta, beta) + beta*fd.f(1.5, eta, beta); // Electron contribution
		numberDensity *= Math.pow(beta, 1.5);
		double frontFactor = 8.0*Math.PI*Math.pow(m_e*c/h, 3.0)*Math.sqrt(2.0);	// Statistical weight of Cox & Giuli.
		numberDensity *= frontFactor;
		return numberDensity;
	}
}
