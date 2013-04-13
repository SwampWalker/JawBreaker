package ca.tonita.jawbreaker.shenzerotemperature.drivers.interpolators;

public abstract class Interpolator {
	/**
	 * Finds the zero of the supplied function.
	 * @param function The data to find the zero of.
	 * @return The double valued index of the zero.
	 */
	public abstract double zeroPoint(double [] function);
	
	/**
	 * Interpolates the function to i, i.e. computes function(i).
	 * For integer i, function(i) = function[i].
	 * @param i The double valued index.
	 * @param function The function to interpolate.
	 * @return The interpolated value.
	 */
	public abstract double interpolate(double i, double[] function);
	
	public abstract String getDescriptor();

	protected double linearZeroPoint(double[] function) {
		double sign = 1;
		if (function[0] == 0) return 0;
		else
			sign = function[0]/Math.abs(function[0]);
		
		// Find the zero crossing.
		int i = 1;
		while (sign*function[i] > 0 && i < function.length-1) i++;
		
		// interpolate
		if (sign*function[i] < 0)
			return i - 1 - function[i-1]/(function[i] - function[i-1]);
		else
			return -1;
	}

	/**
	 * Computes the integer floor of the double.
	 * @param i
	 * @return
	 */
	protected int floor(double i) {
		int floor = 0;
		while (floor + 1 < i) floor++;
		return floor;
	}
}
