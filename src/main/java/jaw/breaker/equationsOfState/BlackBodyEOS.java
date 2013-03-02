package jaw.breaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.PhysicalQuantity;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;

public class BlackBodyEOS {
	private UnitSystem inputUnits = null;
	private UnitSystem myUnits = null;
	// These are here merely for completeness.
	private double k_b	= 1.0;	//1.3806504E-23;
	//private double m_e	= 1.0;	//9.10938215E-31;
	private double c	= 1.0;	//2.99792458E8;
	private double h	= 1.0;	//6.62606896E-34;
	
	public BlackBodyEOS() {
		myUnits =  new UnitSystem(SIConstants.c, SIConstants.melectron, SIConstants.h, SIConstants.k, SIConstants.COULOMB);
	}
	
	public double pressure(double temperature) {
		return a()*Math.pow(temperature, 4.0)/3.0; // checked 27.10.2009
	}
	
	/**
	 * Returns the total internal energy per volume.
	 * @param temperature the temperature of the photon gas
	 * @return the internal energy per unit volume
	 */
	public double internalEnergyDensity(double temperature) {
		return 3.0*pressure(temperature); // checked 27.10.2009
	}
	
	/**
	 * Who likes using weird units? I do. Use this routine to tell me what units you want to give and receive physical quantities in.
	 * @param preferredUnits your preferred units
	 */
	public void setInputUnits(UnitSystem preferredUnits) {
		this.inputUnits = preferredUnits;
		PhysicalQuantity k = new PhysicalQuantity(1.0,SIConstants.k.getDimension(),myUnits);
		k_b = inputUnits.convert(k).getValue();
		PhysicalQuantity c = new PhysicalQuantity(1.0,Dimension.VELOCITY,myUnits);
		this.c = inputUnits.convert(c).getValue();
		PhysicalQuantity h = new PhysicalQuantity(1.0,SIConstants.h.getDimension(),myUnits);
		this.h = inputUnits.convert(h).getValue();
	}
	
	/**
	 * Returns the 'a' parameter, which is the constant multiplying the T^4/3
	 * in the pressure.
	 * @return 4.0*sigma/c
	 */
	private double a() {
		double sigma = 2.0*Math.pow(Math.PI,5)*Math.pow(k_b, 4.0)/(15.0*c*c*h*h*h);
		return 4.0*sigma/c; // checked 27.10.2009
	}
	
	/**
	 * Returns the entropy density of the gas.
	 * @param temperature the temperature of the photon gas
	 * @return the entropy density of the photon gas
	 */
	public double entropyDensity(double temperature) {
		double s = 4.0/3.0*a()*Math.pow(temperature, 3);
		return s; // checked 27.10.2009
	}
	
	/**
	 * Returns the Gibb's free energy density as a function of the temperature.
	 * The free energy density is computed using
	 * the definition of the thermodynamic potential which follows from the
	 * Legendre transformation:<br>
	 * <br>
	 * <code> G = U - TS + PV</code>,<br>
	 * <br>
	 * yielding,<br>
	 * <br>
	 * <code>g = u - Ts + P</code>
	 * @param temperature the temperature of the boson gas
	 * @return the Gibb's free energy density g = u - Ts + P
	 */
	public double gibbsFreeEnergyDensity(double temperature) {
		double internalEnergyDensity = this.internalEnergyDensity(temperature);
		double entropyDensity = this.entropyDensity(temperature);
		double pressure = this.pressure(temperature);
		double freeEnergy = internalEnergyDensity - temperature*entropyDensity + pressure;
		return freeEnergy; // Checked 27.10.2009
	}

	/**
	 * Returns the Helmholtz's free energy density as a function of the temperature.
	 * The free energy density is computed using
	 * the definition of the thermodynamic potential which follows from the
	 * Legendre transformation:<br>
	 * <br>
	 * <code> F = U - TS</code>,<br>
	 * <br>
	 * yielding,<br>
	 * <br>
	 * <code>f = u - Ts</code>
	 * @param temperature the temperature of the boson gas
	 * @return the Helmholtz's free energy density f = u - Ts
	 */
	public double helmholtzFreeEnergyDensity(double temperature) {
		double internalEnergyDensity = this.internalEnergyDensity(temperature);
		double entropyDensity = this.entropyDensity(temperature);
		double freeEnergy = internalEnergyDensity - temperature*entropyDensity;
		return freeEnergy;
	}
}
