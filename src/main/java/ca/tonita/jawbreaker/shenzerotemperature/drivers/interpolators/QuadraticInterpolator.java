package ca.tonita.jawbreaker.shenzerotemperature.drivers.interpolators;

public class QuadraticInterpolator extends Interpolator {

	@Override
	public double zeroPoint(double[] function) {
		double i = linearZeroPoint(function);
		if (i < 0) return i;
		int index = floor(i);
		
		if ((index + 0.5 > i || index == function.length - 2) && index != 0) {
			return leftQuadraticZero(index, function);
		} else {
			return rightQuadraticZero(index, function);
		}
	}

	/**
	 * Computes a quadratic interpolant using the points [index - 1, index, index + 1].
	 * @param index The index to the left of the zero.
	 * @param function The function to find the zero of.
	 * @return The location of the zero according to quadratic interpolation.
	 */
	private double leftQuadraticZero(int index, double[] function) {
		double c = function[index];
		double a = 0.5*(function[index-1] + function[index+1]) - c;
		double b = function[index+1] - a - c;
		double radicand = Math.sqrt(b*b - 4*a*c);
		double rp = (-b + radicand)/(2*a);
		double rm = (-b - radicand)/(2*a);
		
		if ( rp > 0 && rp < 1) {
			return index + rp;
		} else {
			return index + rm;
		}
	}

	/**
	 * Computes a quadratic interpolant using the points [index, index + 1, index + 2].
	 * @param index The index to the left of the zero.
	 * @param function The function to find the zero of.
	 * @return The location of the zero according to quadratic interpolation.
	 */
	private double rightQuadraticZero(int index, double[] function) {
		double c = function[index];
		double a = 0.5*(function[index+2] - 2*function[index+1] + c);
		double b = function[index+1] - a - c;
		double radicand = Math.sqrt(b*b - 4*a*c);
		double rp = (-b + radicand)/(2*a);
		double rm = (-b - radicand)/(2*a);
		
		if ( rp > 0 && rp < 1) {
			return index + rp;
		} else {
			return index + rm;
		}
	}

	@Override
	public double interpolate(double i, double[] function) {
		int index = floor(i);
		
		if ((index + 0.5 > i || index == function.length - 2) && index != 0) {
			return leftInterpolate(i, index, function);
		} else {
			return rightInterpolate(i, index, function);
		}
	}

	private double leftInterpolate(double i, int index, double[] function) {
		double c = function[index];
		double a = 0.5*(function[index-1] + function[index+1]) - c;
		double b = function[index+1] - a - c;
		double chi = i - index;
		double value = c + chi*b + chi*chi*a;
		return value;
	}

	private double rightInterpolate(double i, int index, double[] function) {
		double c = function[index];
		double a = 0.5*(function[index+2] - 2*function[index+1] + c);
		double b = function[index+1] - a - c;
		double chi = i - index;
		double value = c + chi*b + chi*chi*a;
		return value;
	}

	private static String descriptor = "Quadratic interpolation.";
	@Override
	public String getDescriptor() {
		return descriptor;
	}

}
