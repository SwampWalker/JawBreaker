package jaw.breaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.PhysicalQuantity;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;

/**
 * Computes the Shear modulus.
 * @author atonita
 *
 */
public class StrohmayerShearModulus {
	
	/**
	 * Coulomb's constant.
	 */
	private double k; 
	
	/**
	 * The fundamental charge.
	 */
	private double e;
	
	/**
	 * The unit system for inputs and outputs.
	 */
	private UnitSystem units;
	
	/**
	 * Returns the Unit System used for input and output.
	 * @return units The unit system of this class.
	 */
	public UnitSystem getUnits() {
		return units;
	}

	/**
	 * Constructor.
	 */
	public StrohmayerShearModulus() {
		k = SIConstants.k_c.getValue();
		e = SIConstants.e.getValue();
		units = CommonUnits.SI;
	}
	
	/**
	 * Who likes using weird units? I do. Use this routine to tell me what units you want to give and receive physical quantities in.
	 * @param preferredUnits your preferred units
	 */
	public void setUnits(UnitSystem preferredUnits) {
		this.units = preferredUnits;
		PhysicalQuantity x = new PhysicalQuantity(SIConstants.e);
		x = units.convert(x);
		this.e = x.getValue();
		x = new PhysicalQuantity(SIConstants.k_c);
		x = units.convert(x);
		this.k = x.getValue();
	}
	
	/**
	 * This is the T=0 formula for \mu_eff / n from Strohmayer et al. APJ 1991, 375. There
	 * the authors used electrostatic units and the ion density. Here I have converted the
	 * formula to generic units. Note that it has the units of energy (Shear has the units of
	 * pressure, and multiplied by a volume gives energy).
	 * @param density The number density of baryons for the matter.
	 * @param Z The number of protons per ion composing the BCC lattice.
	 * @param A The number of baryons per ion composing the BCC lattice.
	 * @return The effective shear modulus over the number density.
	 */
	public double shearModulusOverNumberDensity(double density, double Z, double A) {
		double a = Math.pow(3/(4*Math.PI*density), 1./3.);
		if (A == 0) return 0; // Z should be zero also... so Z^2/Z -> 0?
		return 0.1194/A*k*(Z*e)*(Z*e)/a;
	}

}
