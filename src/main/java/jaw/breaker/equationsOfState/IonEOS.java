package jaw.breaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.PhysicalQuantity;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;

public class IonEOS {
	private double abar;
	private UnitSystem inputUnits = null;
	private UnitSystem myUnits = null;
	// These are here merely for completeness.
	private double k_b	= 1.0;	//1.3806504E-23;
	//private double m_e	= 1.0;	//9.10938215E-31;
	//private double c	= 1.0;	//2.99792458E8;
	//private double h	= 1.0;	//6.62606896E-34;
	
	public IonEOS(double abar) {
		this.abar = abar;
		myUnits =  new UnitSystem(SIConstants.c, SIConstants.melectron, SIConstants.h, SIConstants.k, SIConstants.COULOMB);
	}
	
	public double pressure(double numberDensity, double temperature) {
		return numberDensity*k_b*temperature;
	}
	
	public double internalEnergyDensity(double numberDensity, double temperature) {
		return 1.5*pressure(numberDensity,temperature);
	}
	
	public double numberDensity(double massDensity) {
		return SIConstants.AVOGADROSNUMBER.getValue()*massDensity/abar;
	}
	
	/**
	 * Who likes using weird units? I do. Use this routine to tell me what units you want to give and receive physical quantities in.
	 * @param preferredUnits your preferred units
	 */
	public void setInputUnits(UnitSystem preferredUnits) {
		this.inputUnits = preferredUnits;
		PhysicalQuantity k = new PhysicalQuantity(1.0,SIConstants.k.getDimension(),myUnits);
		k_b = inputUnits.convert(k).getValue();
	}
}
