package ca.tonita.jawbreaker.shenzerotemperature.drivers;

import java.io.File;

public interface InputOutput {
	public static final int SUCCESS = 0;
	public static final int UNNAMEDERROR = -1;
	public static final int INCORRECTINPUT = -2;
	public int readAndWrite(File input, File output);
	
	public static int density = 0;
	public static int numberdensity = 1;
	public static int protonFraction = 2;
	public static int freeEnergy = 3;
	public static int internalEnergy = 4;
	public static int entropy = 5;
	public static int massNumber = 6;
	public static int chargeNumber = 7;
	public static int effectiveMass = 8;
	public static int freeNeutronFraction = 9;
	public static int freeProtonFraction = 10;
	public static int alphaFraction = 11;
	public static int ionFraction = 12;
	public static int pressure = 13;
	public static int neutronPotential = 14;
	public static int protonPotential = 15;
	public static int hyperonMass = 16;
	public static int hyperonFraction = 17;
	public static int electronPotential = 18;
	public static int nFields = 19;
	
	public static int shear = 19;
}
