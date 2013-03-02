package jaw.breaker.shenzerotemperature.drivers.interpolators;

public class LinearInterpolator extends Interpolator {

	@Override
	public double zeroPoint(double[] function) {
		return linearZeroPoint(function);
	}
	
	@Override
	public double interpolate(double i, double[] function) {
		int lowerIndex = floor(i);
		double chi = i - lowerIndex;
		double value = function[lowerIndex] + chi*(function[lowerIndex+1] - function[lowerIndex]);
		return value;
	}

	private final static String descriptor = "Linear interpolator.";
	public String getDescriptor() {
		return descriptor;
	}

}