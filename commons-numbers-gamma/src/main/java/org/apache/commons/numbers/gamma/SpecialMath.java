/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.commons.numbers.gamma;

import org.apache.commons.numbers.core.DD;

/**
 * Special math functions.
 *
 * @since 1.1
 */
final class SpecialMath {
    /** Threshold for z in the Taylor series approximation of expm1(z * t). */
    static final double EXT_POW_LIMIT = 0x1p26;
    /** Minimum x for log1pmx(x). */
    private static final double X_MIN = -1;
    /** Low threshold to use log1p(x) - x. */
    private static final double X_LOW = -0.79149064;
    /** High threshold to use log1p(x) - x. */
    private static final double X_HIGH = 1;
    /** 2^-6. */
    private static final double TWO_POW_M6 = 0x1.0p-6;
    /** 2^-12. */
    private static final double TWO_POW_M12 = 0x1.0p-12;
    /** 2^-20. */
    private static final double TWO_POW_M20 = 0x1.0p-20;
    /** 2^-53. */
    private static final double TWO_POW_M53 = 0x1.0p-53;
    /** ln(0.5). */
    private static final double LN_0_5 = Math.log(0.5);
    /** ln(1.5). */
    private static final double LN_1_5 = Math.log(1.5);

    /** Private constructor. */
    private SpecialMath() {
        // intentionally empty.
    }

    /**
     * Returns {@code log(1 + x) - x}. This function is accurate when {@code x -> 0}.
     *
     * <p>This function uses a Taylor series expansion when x is small ({@code |x| < 0.01}):
     *
     * <pre>
     * ln(1 + x) - x = -x^2/2 + x^3/3 - x^4/4 + ...
     * </pre>
     *
     * <p>or around 0 ({@code -0.791 <= x <= 1}):
     *
     * <pre>
     * ln(x + a) = ln(a) + 2 [z + z^3/3 + z^5/5 + z^7/7 + ... ]
     *
     * z = x / (2a + x)
     * </pre>
     *
     * <p>For a = 1:
     *
     * <pre>
     * ln(x + 1) - x = -x + 2 [z + z^3/3 + z^5/5 + z^7/7 + ... ]
     *               = z * (-x + 2z^2 [ 1/3 + z^2/5 + z^4/7 + ... ])
     * </pre>
     *
     * <p>The code is based on the {@code log1pmx} documentation for the <a
     * href="https://rdrr.io/rforge/DPQ/man/log1pmx.html">R DPQ package</a> with addition of the
     * direct Taylor series for tiny x.
     *
     * <p>See Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical Functions. New York:
     * Dover. Formulas 4.1.24 and 4.2.29, p.68. <a
     * href="https://en.wikipedia.org/wiki/Abramowitz_and_Stegun">Wikipedia: Abramowitz_and_Stegun</a>
     * provides links to the full text which is in public domain.
     *
     * @param x Value x
     * @return {@code log(1 + x) - x}
     */
    static double log1pmx(double x) {
        // -1 is the minimum supported value
        if (x <= X_MIN) {
            return x == X_MIN ? Double.NEGATIVE_INFINITY : Double.NaN;
        }
        // Use the thresholds documented in the R implementation
        if (x < X_LOW || x > X_HIGH) {
            return Math.log1p(x) - x;
        }
        final double a = Math.abs(x);

        // Addition to the R version for small x.
        // Use a direct Taylor series:
        // ln(1 + x) = x - x^2/2 + x^3/3 - x^4/4 + ...
        if (a < TWO_POW_M6) {
            return log1pmxSmall(x, a);
        }

        // The use of the following series is fast converging:
        // ln(x + 1) - x = -x + 2 [z + z^3/3 + z^5/5 + z^7/7 + ... ]
        //               = z * (-x + 2z^2 [ 1/3 + z^2/5 + z^4/7 + ... ])
        // z = x / (2 + x)
        //
        // Tests show this is more accurate when |x| > 1e-4 than the direct Taylor series.
        // The direct series can be modified to sum multiple terms together for a small
        // increase in precision to a closer match to this variation but the direct series
        // takes approximately 3x longer to converge.

        final double z = x / (2 + x);
        final double zz = z * z;

        // Series sum
        // sum(k=0,...,Inf; zz^k/(3+k*2)) = 1/3 + zz/5 + zz^2/7 + zz^3/9 + ... )

        double sum = 1.0 / 3;
        double numerator = 1;
        int denominator = 3;
        for (;;) {
            numerator *= zz;
            denominator += 2;
            final double sum2 = sum + numerator / denominator;
            // Since |x| <= 1 the additional terms will reduce in magnitude.
            // Iterate until convergence. Expected iterations:
            // x      iterations
            // -0.79  38
            // -0.5   15
            // -0.1    5
            //  0.1    5
            //  0.5   10
            //  1.0   15
            if (sum2 == sum) {
                break;
            }
            sum = sum2;
        }
        return z * (2 * zz * sum - x);
    }

    /**
     * Returns {@code log(1 + x) - x}. This function is accurate when
     * {@code x -> 0}.
     *
     * <p>This function uses a Taylor series expansion when x is small
     * ({@code |x| < 0.01}):
     *
     * <pre>
     * ln(1 + x) - x = -x^2/2 + x^3/3 - x^4/4 + ...
     * </pre>
     *
     * <p>No loop iterations are used as the series is directly expanded
     * for a set number of terms based on the absolute value of x.
     *
     * @param x Value x (assumed to be small)
     * @param a Absolute value of x
     * @return {@code log(1 + x) - x}
     */
    private static double log1pmxSmall(double x, double a) {
        // Use a direct Taylor series:
        // ln(1 + x) = x - x^2/2 + x^3/3 - x^4/4 + ...
        // Reverse the summation (small to large) for a marginal increase in precision.
        // To stop the Taylor series the next term must be less than 1 ulp from the
        // answer.
        // x^n/n < |log(1+x)-x| * eps
        // eps = machine epsilon = 2^-53
        // x^n < |log(1+x)-x| * eps
        // n < (log(|log(1+x)-x|) + log(eps)) / log(x)
        // In practice this is a conservative limit.

        final double x2 = x * x;

        if (a < TWO_POW_M53) {
            // Below machine epsilon. Addition of x^3/3 is not possible.
            // Subtract from zero to prevent creating -0.0 for x=0.
            return 0 - x2 / 2;
        }

        final double x4 = x2 * x2;

        // +/-9.5367431640625e-07: log1pmx = -4.547470617660916e-13 :
        // -4.5474764000725028e-13
        // n = 4.69
        if (a < TWO_POW_M20) {
            // n=5
            return x * x4 / 5 -
                       x4 / 4 +
                   x * x2 / 3 -
                       x2 / 2;
        }

        // +/-2.44140625E-4: log1pmx = -2.9797472637290841e-08 : -2.9807173914456693e-08
        // n = 6.49
        if (a < TWO_POW_M12) {
            // n=7
            return x * x2 * x4 / 7 -
                       x2 * x4 / 6 +
                        x * x4 / 5 -
                            x4 / 4 +
                        x * x2 / 3 -
                            x2 / 2;
        }

        // Assume |x| < 2^-6
        // +/-0.015625: log1pmx = -0.00012081346403474586 : -0.00012335696813916864
        // n = 10.9974

        // n=11
        final double x8 = x4 * x4;
        return x * x2 * x8 / 11 -
                   x2 * x8 / 10 +
                    x * x8 /  9 -
                        x8 /  8 +
               x * x2 * x4 /  7 -
                   x2 * x4 /  6 +
                    x * x4 /  5 -
                        x4 /  4 +
                    x * x2 /  3 -
                        x2 /  2;
    }

    /**
     * Checks if the value is odd.
     *
     * <p>Helper function used in Boost gamma functions and the Hurwitz zeta function.
     *
     * @param v Value (assumed to be positive finite and an integer)
     * @return true if odd
     * @since 1.4
     */
    static boolean isOdd(double v) {
        // Note:
        // Any value larger than 2^53 should be even.
        // If the input is positive then truncation of extreme doubles (>2^63)
        // to the primitive long creates an odd value: 2^63-1.
        // This is corrected by inverting the sign of v and the extreme is even: -2^63.
        // The isOdd function is correct for all positive finite v.
        return (((long) -v) & 0x1) == 1;
    }

    /**
     * Extended precision power function {@code (x+n)^z}.
     *
     * <p>This method is used in the zeta function to compute the power function terms of
     * the infinite series {@code (a+k)^-s}. It is limited to the expected arguments
     * of the zeta function use case:
     *
     * <ul>
     * <li>{@code x > -1}
     * <li>{@code n >= 1}
     * <li>{@code x+n > 0} and finite
     * <li>{@code z < -1} and finite
     * </ul>
     *
     * <p>This method has no advantage over {@code Math.pow(x+n, z)} when there is no
     * round-off component to (x+n). In the event that {@code (x+n)} is exact this method
     * will return the same result as {@code Math.pow(x+n, z)}.
     *
     * <p>Calling with {@code n=0} will consume resources for no benefit. More generally
     * this method should not be called when {@code x + n} does not create a change of
     * exponent between {@code x} and {@code x + n}, and may offer no benefit even if the
     * exponent does change depending on the arguments. For example if {@code |x| < 2^52}
     * and an integer then addition of any {@code int} value will be exact.
     *
     * <p>Note: This method has been developed for use with integer {@code n} but can
     * be adapted to a {@code double} argument {@code n}; and for a larger domain of
     * input arguments.
     *
     * @param x Argument x ({@code 0 <= |x| < 2^53})
     * @param n Argument n ({@code > 0})
     * @param z Exponent z ({@code < -1})
     * @return the result
     * @since 1.4
     */
    static double powxpn(double x, int n, double z) {
        // (s+ss)^z = s^z * (1+ss/s)^z
        //          = s^z * exp(z*log1p(ss/s))
        //          = s^z * s^z * [ (1+ss/s)^z - 1 ]
        //          = s^z + s^z * expm1(z*log1p(ss/s))
        //
        // ss/s < 1 * machine epsilon : log1p(ss/s) = ss/s
        // ss/s < 2 * machine epsilon : log1p(ss/s) ~ ss/s
        //
        // Taylor series:
        // (1 + x)^n = 1 + nx + n(n-1)x^2/2! + n(n-1)(n-2)x^3/3! + ...
        // for |x| < 1 and real n

        // Compute the exact double-double sum.
        // Avoid generation of a DD intermediate by directly using twoSum.
        double s = x + n;
        double ss = twoSumLow(x, n, s);

        // Skip impossible rounding
        if (ss == 0) {
            // Handles (x+n) = 0 or 1 for the correct pow result
            return Math.pow(s, z);
        }

        double t;

        // Edge cases:
        // s < 1 can overflow to infinity
        // s = 1 creates 1.0
        // s > 1 can underflow to 0
        // Care must be taken to avoid 0 * infinity == NaN

        // Prevent overflow to infinity if s < 1 && ss > 0 by
        // increasing s by 1 ulp and changing ss to negative by subtracting the ulp.

        // Prevent underflow to zero if s > 1 && ss < 0 by
        // decreasing s by 1 ulp and changing ss to positive by adding the ulp.

        // Adjust s to be closer to 1.
        // This does increase accuracy in extreme cases but:
        // - the effect on zeta is not measurable
        // - the branch statements are unpredictable
//        if (s < 1 && ss > 0) {
//            t = Math.ulp(s);
//            s += t;
//            ss -= t;
//        } else if (s > 1 && ss < 0) {
//            t = Math.ulp(s);
//            s -= t;
//            ss += t;
//        }

        double r = Math.pow(s, z);
        // Do not adjust infinity.
        // NaNs or infinite are not expected from the calling zeta function.
        if (Double.isInfinite(r)) {
            return r;
        }

        // 0 < t < 2^-52
        t = ss / s;

        // Limit of Taylor series where it does not equal expm1(z * t)
        // is at approximately 2^-26.
        // z * t < 2^-26  =>  z < 2^-26 / 2^-52 = 2^26
        if (z < -EXT_POW_LIMIT) {
            // z is expected to be small (practical range for a zeta result is -(1, 1075))
            // and this path is unlikely. Use either:
            //   s^z * exp(z*t)
            //   s^z + s^z * expm1(z*t)
            // Pick expm1 when its result is in the range -0.5 to 0.5.
            // Otherwise use exp (when its result is in [0, 0.5] and [1.5, infinity].
            final double zt = z * t;
            if (zt < LN_0_5 || zt > LN_1_5) {
                // Result can be 0 or infinity.
                // Avoid 0 * infinity (infinite r has been handled).
                t = Math.exp(z * t);
                if (t == Double.POSITIVE_INFINITY && r == 0) {
                    // r == 0 when s < 1 and z is large negative
                    // t == infinity when ss < 0 and z is large negative
                    // assume (s,ss) < 1 and return zero
                    return 0;
                }
                return r * t;
            }
            // Fall-through with the (1 + x)^n - 1 result
            t = Math.expm1(z * t);
        } else {
            // Compute (1 + x)^n - 1 using a Taylor series.
            // This is the most common branch path when ss != 0.
            t = z * t * (1 + (z - 1) * t * 0.5);
        }

        // t in [-0.5, 0.5]
        return r + r * t;
    }

    /**
     * Compute the round-off of the sum of two numbers {@code a} and {@code b} using
     * Knuth two-sum algorithm. The values are not required to be ordered by magnitude,
     * i.e. the result is commutative {@code s = a + b == b + a}.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return the sum round-off
     */
    private static double twoSumLow(double a, double b, double x) {
        // (x, xx) = a + b
        // bVirtual = x - a
        // aVirtual = x - bVirtual
        // bRoundoff = b - bVirtual
        // aRoundoff = a - aVirtual
        // xx = aRoundoff + bRoundoff
        final double bVirtual = x - a;
        return (a - (x - bVirtual)) + (b - bVirtual);
    }
}
