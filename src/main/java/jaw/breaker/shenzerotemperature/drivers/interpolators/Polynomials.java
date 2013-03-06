package jaw.breaker.shenzerotemperature.drivers.interpolators;

import ta.tonita.math.linearalgebra.LinearAlgebra;

public class Polynomials {

    /**
     * Returns the coefficients of the interpolating polynomial of the data f.
     * Assumes that f[0] is at x=0.
     *
     * @param f The data to interpolate.
     * @return The coefficients of x^0, x^1, ..., x^n.
     */
    public static double[] interpolatingCoefficients(double[] f) {
        double[] coeffs = new double[f.length];
        double[][] matrix = new double[f.length][f.length];

        // Set up the data to interpolate.
        for (int iRow = 0; iRow < f.length; iRow++) {
            coeffs[iRow] = f[iRow];
            double value = 1;
            for (int iCol = 0; iCol < f.length; iCol++) {
                matrix[iRow][iCol] = value;
                value *= iRow;
            }
        }

        LinearAlgebra.solve(matrix, coeffs);

        return coeffs;
    }

    /**
     * Returns the coefficients of the interpolating polynomial of the data f.
     * Returns data so that f_interp(x[0]) = f[0].
     *
     * @param f The data to interpolate.
     * @param x The coordinates of the points.
     * @return The coefficients of x^0, x^1, ..., x^n.
     */
    public static double[] interpolatingCoefficients(double[] f, double[] x) {
        double[] coeffs = new double[f.length];
        double[][] matrix = new double[f.length][f.length];

        // Set up the data to interpolate.
        for (int iRow = 0; iRow < f.length; iRow++) {
            coeffs[iRow] = f[iRow];
            double value = 1;
            for (int iCol = 0; iCol < f.length; iCol++) {
                if (iCol > 0) {
                    value *= x[iRow];
                }
                matrix[iRow][iCol] = value;
            }
        }

        LinearAlgebra.solve(matrix, coeffs);

        return coeffs;
    }

    /**
     * Returns the coefficients of the interpolating polynomial of the data f.
     * Returns data so that f_interp(x[0]) = f[0].
     *
     * @param f The n data points to interpolate.
     * @param df The n derivatives of the data to interpolate.
     * @param x The coordinates of the points.
     * @return The coefficients of x^0, x^1, ..., x^2n.
     */
    public static double[] interpolatingCoefficients(double[] f, double[] df, double[] x) {
        double[] coeffs = new double[f.length * 2];
        double[][] matrix = new double[f.length * 2][f.length * 2];

        // Set up the function value rows..
        for (int iRow = 0; iRow < f.length; iRow++) {
            coeffs[iRow] = f[iRow];
            double value = 1;
            for (int iCol = 0; iCol < f.length * 2; iCol++) {
                if (iCol > 0) {
                    value *= x[iRow];
                }
                matrix[iRow][iCol] = value;
            }
        }
        // Set up the derivative value rows..
        for (int i = 0; i < f.length; i++) {
            int iRow = f.length + i;
            coeffs[iRow] = df[i];
            double coordValue = 0;
            for (int iCol = 0; iCol < f.length * 2; iCol++) {
                if (iCol > 1) {
                    coordValue *= x[i];
                } else if (iCol == 1) {
                    coordValue = iCol;
                }
                matrix[iRow][iCol] = iCol*coordValue;
            }
        }

        LinearAlgebra.solve(matrix, coeffs);

        return coeffs;
    }

    /**
     * Returns the interpolated values, given the coefficients and the position.
     *
     * @param coeffs The coefficients of the interpolant.
     * @param x The coordinate to evaluate the polynomial at.
     * @return The function evaluate at x.
     */
    public static double interpolate(double[] coeffs, double x) {
        double retval = 0;
        double powX = 1;
        for (int i = 0; i < coeffs.length; i++) {
            retval += coeffs[i] * powX;
            powX *= x;
        }
        return retval;
    }

    /**
     * Returns the differentiated value, given the coefficients and the
     * position.
     *
     * @param coeffs The coefficients of the interpolant.
     * @param x The coordinate to evaluate the polynomial at.
     * @return The function differentiated at x.
     */
    public static double differentiate(double[] coeffs, double x) {
        double retval = 0;
        double powX = 1;
        for (int i = 1; i < coeffs.length; i++) {
            retval += i * coeffs[i] * powX;
            powX *= x;
        }
        return retval;
    }

    /**
     * Returns the differentiated value, given the coefficients and the
     * position.
     *
     * @param coeffs The coefficients of the interpolant.
     * @param x The coordinate to evaluate the polynomial at.
     * @param dx The physical distance between points.
     * @return The function differentiated at x.
     */
    public static double differentiate(double[] coeffs, double x, double dx) {
        return differentiate(coeffs, x) / dx;
    }
}
