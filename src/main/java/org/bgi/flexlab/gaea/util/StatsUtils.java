package org.bgi.flexlab.gaea.util;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Locale;

public class StatsUtils {

    public static String perc(long numerator, long denominator) {
        DecimalFormat df = new DecimalFormat("0.00");
        df.setRoundingMode(RoundingMode.HALF_UP);
        final StringBuilder sb = new StringBuilder();
        if (denominator == 0) {
            sb.append("-");
        } else {
            sb.append(realFormat((double) numerator / denominator * 100, 1));
        }
        return sb.toString();
    }

    public static String percent(long numerator, long denominator) {
        DecimalFormat df = new DecimalFormat("0.00");
        df.setRoundingMode(RoundingMode.HALF_UP);
        final StringBuilder sb = new StringBuilder();
        if (denominator == 0) {
            sb.append("-");
        } else {
            sb.append(realFormat((double) numerator / denominator * 100, 1)).append("%");
        }
        sb.append(" (").append(numerator).append("/").append(denominator).append(")");
        return sb.toString();
    }

    public static String divide(long numerator, long denominator) {
        DecimalFormat df = new DecimalFormat("0.00");
        df.setRoundingMode(RoundingMode.HALF_UP);
        final StringBuilder sb = new StringBuilder();
        if (denominator == 0) {
            sb.append("-");
        } else {
            sb.append(realFormat((double) numerator / denominator, 2));
        }
        sb.append(" (").append(numerator).append("/").append(denominator).append(")");
        return sb.toString();
    }

    /**
     * Check if a string is of the form "-0.00".
     *
     * @param str the string being checked.
     * @param dp the number of zeros after the decimal point.
     * @return true iff of the the form being checked.
     */
    static boolean negZero(final String str, final int dp) {
        if (str.charAt(0) != '-' || str.charAt(1) != '0') {
            return false;
        } else if (dp == 0) {
            return true;
        } else if (str.charAt(2) != '.') {
            return false;
        }
        for (int i = 0; i < dp; ++i) {
            if (str.charAt(3 + i) != '0') {
                return false;
            }
        }
        return true;
    }

    /**
     * Format a real with specified number of digits after the decimal point.
     * @param x number to be formatted.
     * @param dp number of digits of precision.
     * @return the formatted string.
     */
    public static String realFormat(final double x, final int dp) {
        assert dp >= 0;
        if (Double.isNaN(x) || Double.isInfinite(x)) {
            return Double.toString(x);
        }
        final String fmt = "%1$01." + dp + "f";
        final String res = String.format(Locale.ROOT, fmt, x);
        if (x <= 0.0) {
            if (negZero(res, dp)) {
                return res.substring(1);
            } else {
                return res;
            }
        } else {
            return res;
        }
    }
}
