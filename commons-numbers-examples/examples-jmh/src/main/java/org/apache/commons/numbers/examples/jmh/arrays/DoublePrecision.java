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
package org.apache.commons.numbers.examples.jmh.arrays;

/**
 * Computes double-length precision floating-point operations.
 *
 * <p>It is based on the 1971 paper
 * <a href="https://doi.org/10.1007/BF01397083">
 * Dekker (1971) A floating-point technique for extending the available precision</a>.
 */
final class DoublePrecision {
    /*
     * Caveat:
     *
     * The code below uses many additions/subtractions that may
     * appear redundant. However, they should NOT be simplified, as they
     * do use IEEE754 floating point arithmetic rounding properties.
     *
     * Algorithms are based on computing the product or sum of two values x and y in
     * extended precision. The standard result is stored using a double (high part z) and
     * the round-off error (or low part zz) is stored in a second double, e.g:
     * x * y = (z, zz); z + zz = x * y
     * x + y = (z, zz); z + zz = x + y
     *
     * To sum multiple (z, zz) results ideally the parts are sorted in order of
     * non-decreasing magnitude and summed. This is exact if each number's most significant
     * bit is below the least significant bit of the next (i.e. does not
     * overlap). Creating non-overlapping parts requires a rebalancing
     * of adjacent pairs using a summation z + zz = (z1, zz1) iteratively through the parts
     * (see Shewchuk (1997) Grow-Expansion and Expansion-Sum [1]).
     *
     * [1] Shewchuk (1997): Arbitrary Precision Floating-Point Arithmetic
     * http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
     */

    /**
     * The multiplier used to split the double value into high and low parts. From
     * Dekker (1971): "The constant should be chosen equal to 2^(p - p/2) + 1,
     * where p is the number of binary digits in the mantissa". Here p is 53
     * and the multiplier is {@code 2^27 + 1}.
     */
    private static final double MULTIPLIER = 1.34217729E8;

    /** The upper limit above which a number may overflow during the split into a high part.
     * Assuming the multiplier is above 2^27 and the maximum exponent is 1023 then a safe
     * limit is a value with an exponent of (1023 - 27) = 2^996. */
    private static final double SAFE_UPPER = 0x1.0p996;

    /** The scale to use when down-scaling during a split into a high part.
     * This must be smaller than the inverse of the multiplier and a power of 2 for exact scaling. */
    private static final double DOWN_SCALE = 0x1.0p-30;

    /** The scale to use when re-scaling during a split into a high part.
     * This is the inverse of {@link #DOWN_SCALE}. */
    private static final double UP_SCALE = 0x1.0p30;

    /** The mask to zero the lower 27-bits of a long . */
    private static final long ZERO_LOWER_27_BITS = 0xffff_ffff_f800_0000L;

    /**
     * Compute the low part of the double length number {@code (z,zz)} for the exact
     * product of {@code x} and {@code y}. The low part of the result {@code zz} is computed
     * by splitting the input numbers into parts to allow partial products to be computed.
     * Each instance uses a different split method.
     */
    enum Product {
        /**
         * Implement Dekker's method to split a value into two parts.
         * Uses scaling to avoid overflow in intermediate computations.
         */
        DEKKER {
            /**
             * {@inheritDoc}
             *
             * <p>Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) creates
             * a big value from which to derive the two split parts.
             * <pre>
             * c = (2^s + 1) * a
             * a_big = c - a
             * a_hi = c - a_big
             * a_lo = a - a_hi
             * a = a_hi + a_lo
             * </pre>
             *
             * <p>The multiplicand allows a p-bit value to be split into
             * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
             * Combined they have (p􏰔-1) bits of significand but the sign bit of {@code a_lo}
             * contains a bit of information. The constant is chosen so that s is ceil(p/2) where
             * the precision p for a double is 53-bits (1-bit of the mantissa is assumed to be
             * 1 for a non sub-normal number) and s is 27.
             *
             * <p>This conversion uses scaling to avoid overflow in intermediate computations.
             *
             * <p>Splitting a NaN or infinite value will return NaN. Any finite value will return
             * a finite value.
             */
            @Override
            double highPart(double a) {
                // Avoid overflow
                if (Math.abs(a) >= SAFE_UPPER) {
                    // Do scaling.
                    final double x = a * DOWN_SCALE;
                    final double c = MULTIPLIER * x;
                    final double hi = (c - (c - x)) * UP_SCALE;
                    if (Double.isInfinite(hi)) {
                        // Number is too large.
                        // This occurs if value is infinite or close to Double.MAX_VALUE.
                        // Note that Dekker's split creates an approximating 26-bit number which may
                        // have an exponent 1 greater than the input value. This will overflow if the
                        // exponent is already +1023. Revert to the raw upper 26 bits of the 53-bit
                        // mantissa (including the assumed leading 1 bit). This conversion will result in
                        // the low part being a 27-bit significand and the potential loss of bits during
                        // addition and multiplication. (Contrast to the Dekker split which creates two
                        // 26-bit numbers with a bit of information moved to the sign of low.)
                        // The conversion will maintain Infinite in the high part where the resulting
                        // low part a_lo = a - a_hi = inf - inf = NaN.
                        return Double.longBitsToDouble(Double.doubleToRawLongBits(a) & ZERO_LOWER_27_BITS);
                    }
                    return hi;
                }
                // normal conversion
                final double c = MULTIPLIER * a;
                return c - (c - a);
            }
        },
        /**
         * Implement Dekker's method to split a value into two parts.
         *
         * <p>Warning: This method does not perform scaling in Dekker's split and large
         * finite numbers can create NaN results.
         */
        FAST_DEKKER {
            /**
             * {@inheritDoc}
             *
             * <p>Implement Dekker's method to split a value into two parts (see {@link Product#DEKKER}).
             *
             * <p>Warning: This method does not perform scaling in Dekker's split and large
             * finite numbers can create NaN results. Overflow may occur when the exponent of
             * the input value is above 996.
             *
             * <p>Splitting a NaN or infinite value will return NaN. Any finite value will return
             * a finite value.
             */
            @Override
            double highPart(double a) {
                // normal conversion without overflow checks
                final double c = MULTIPLIER * a;
                return c - (c - a);
            }
        },
        /**
         * Implement a split using the upper and lower raw bits from the value.
         *
         * <p>Note: This method will not work for very small sub-normal numbers
         * ({@code <= 27} bits) as the high part will be zero and the low part will
         * have all the information. Methods that assume {@code hi > lo} will have
         * undefined behaviour.
         */
        SPLIT {
            /**
             * {@inheritDoc}
             *
             * <p>Split using the upper and lower raw bits from the double.
             *
             * <p>Splitting a NaN value will return NaN or infinite. Splitting an infinite
             * value will return infinite. Any finite value will return a finite value.
             */
            @Override
            double highPart(double a) {
                return Double.longBitsToDouble(Double.doubleToRawLongBits(a) & ZERO_LOWER_27_BITS);
            }
        };

        /**
         * Compute the low part of the double length number {@code (z,zz)} for the exact
         * product of {@code x} and {@code y}. The standard precision product {@code x*y}
         * must be provided. The numbers {@code x} and {@code y}
         * are split into high and low parts for an extended precision algorithm:
         * <pre>
         * lx * ly - (((xy - hx * hy) - lx * hy) - hx * ly)
         * </pre>
         *
         * <p>Note: This uses the high part of the result {@code (z,zz)} as {@code x * y} and not
         * {@code hx * hy + hx * ty + tx * hy} as specified in Dekker's original paper.
         * See Shewchuk (1997) for working examples.
         *
         * <p>Warning: Dekker's split can produce high parts that are larger in magnitude than
         * the input number as the high part is a 26-bit approximation of the number. Thus it is
         * possible that the standard product {@code x * y} does not overflow but the extended
         * precision sub-product {@code hx * hy} does overflow.
         *
         * @param x First factor.
         * @param y Second factor.
         * @param xy Product of the factors (x * y).
         * @return the low part of the product double length number
         * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
         * Shewchuk (1997) Theorum 18</a>
         */
        double low(double x, double y, double xy) {
            // Split the numbers using Dekker's algorithm
            final double hx = highPart(x);
            final double lx = x - hx;

            final double hy = highPart(y);
            final double ly = y - hy;

            // Compute the multiply low part:
            // err1 = xy - hx * hy
            // err2 = err1 - lx * hy
            // err3 = err2 - hx * ly
            // low = lx * ly - err3
            return lx * ly - (((xy - hx * hy) - lx * hy) - hx * ly);
        }

        /**
         * Split a value into the high-part {@code hi} of a two part number. Subtracting the
         * high part from the original number creates the low part.
         * <pre>
         * lo = a - hi
         * a = hi + lo
         * </pre>
         *
         * @param a Value.
         * @return the high part of the value.
         */
        abstract double highPart(double a);
    }

    /** Private constructor. */
    private DoublePrecision() {
        // intentionally empty.
    }

    /**
     * Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) creates
     * a big value from which to derive the two split parts.
     * <pre>
     * c = (2^s + 1) * a
     * a_big = c - a
     * a_hi = c - a_big
     * a_lo = a - a_hi
     * a = a_hi + a_lo
     * </pre>
     *
     * <p>The multiplicand allows a p-bit value to be split into
     * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
     * Combined they have (p􏰔-1) bits of significand but the sign bit of {@code a_lo}
     * contains a bit of information. The constant is chosen so that s is ceil(p/2) where
     * the precision p for a double is 53-bits (1-bit of the mantissa is assumed to be
     * 1 for a non sub-normal number) and s is 27.
     *
     * <p>This conversion uses scaling to avoid overflow in intermediate computations.
     *
     * <p>Splitting a NaN or infinite value will return NaN. Any finite value will return
     * a finite value.
     *
     * @param value Value.
     * @return the high part of the value.
     */
    static double highPart(double value) {
        // Avoid overflow
        if (Math.abs(value) >= SAFE_UPPER) {
            // Do scaling.
            final double hi = highPartUnscaled(value * DOWN_SCALE) * UP_SCALE;
            if (Double.isInfinite(hi)) {
                // Number is too large.
                // This occurs if value is infinite or close to Double.MAX_VALUE.
                // Note that Dekker's split creates an approximating 26-bit number which may
                // have an exponent 1 greater than the input value. This will overflow if the
                // exponent is already +1023. Revert to the raw upper 26 bits of the 53-bit
                // mantissa (including the assumed leading 1 bit). This conversion will result in
                // the low part being a 27-bit significand and the potential loss of bits during
                // addition and multiplication. (Contrast to the Dekker split which creates two
                // 26-bit numbers with a bit of information moved to the sign of low.)
                // The conversion will maintain Infinite in the high part where the resulting
                // low part a_lo = a - a_hi = inf - inf = NaN.
                return highPartSplit(value);
            }
            return hi;
        }
        // normal conversion
        return highPartUnscaled(value);
    }

    /**
     * Implement Dekker's method to split a value into two parts (see {@link #highPart(double)}).
     *
     * <p>This conversion does not use scaling and the result of overflow is NaN. Overflow
     * may occur when the exponent of the input value is above 996.
     *
     * <p>Splitting a NaN or infinite value will return NaN.
     *
     * @param value Value.
     * @return the high part of the value.
     * @see Math#getExponent(double)
     */
    static double highPartUnscaled(double value) {
        final double c = MULTIPLIER * value;
        return c - (c - value);
    }

    /**
     * Implement a split using the upper and lower raw bits from the value.
     *
     * <p>Note: This method will not work for very small sub-normal numbers
     * ({@code <= 27} bits) as the high part will be zero and the low part will
     * have all the information. Methods that assume {@code hi > lo} will have
     * undefined behaviour.
     *
     * <p>Splitting a NaN value will return NaN or infinite. Splitting an infinite
     * value will return infinite. Any finite value will return a finite value.
     *
     * @param value Value.
     * @return the high part of the value.
     * @see Math#getExponent(double)
     */
    static double highPartSplit(double value) {
        return Double.longBitsToDouble(Double.doubleToRawLongBits(value) & ZERO_LOWER_27_BITS);
    }

    /**
     * Compute the low part of the double length number {@code (z,zz)} for the exact
     * product of {@code x} and {@code y} using Dekker's mult12 algorithm. The standard
     * precision product {@code x*y} must be provided. The numbers {@code x} and {@code y}
     * are split into high and low parts using Dekker's algorithm.
     *
     * <p>Warning: Dekker's split can produce high parts that are larger in magnitude than
     * the input number as the high part is a 26-bit approximation of the number. Thus it is
     * possible that the standard product {@code x * y} does not overflow but the extended
     * precision sub-product {@code hx * hy} does overflow.
     *
     * @param x First factor.
     * @param y Second factor.
     * @param xy Product of the factors (x * y).
     * @return the low part of the product double length number
     * @see #highPart(double)
     * @see #productLow(double, double, double, double, double)
     */
    static double productLow(double x, double y, double xy) {
        // Split the numbers using Dekker's algorithm
        final double hx = highPart(x);
        final double lx = x - hx;

        final double hy = highPart(y);
        final double ly = y - hy;

        return productLow(hx, lx, hy, ly, xy);
    }

    /**
     * Compute the low part of the double length number {@code (z,zz)} for the exact
     * product of {@code x} and {@code y} using Dekker's mult12 algorithm. The standard
     * precision product {@code x*y} must be provided. The numbers {@code x} and {@code y}
     * are split into high and low parts using Dekker's algorithm.
     *
     * <p>Warning: This method does not perform scaling in Dekker's split and large
     * finite numbers can create NaN results.
     *
     * @param x First factor.
     * @param y Second factor.
     * @param xy Product of the factors (x * y).
     * @return the low part of the product double length number
     * @see #highPartUnscaled(double)
     * @see #productLow(double, double, double, double, double)
     */
    static double productLowUnscaled(double x, double y, double xy) {
        // Split the numbers using Dekker's algorithm without scaling
        final double hx = highPartUnscaled(x);
        final double lx = x - hx;

        final double hy = highPartUnscaled(y);
        final double ly = y - hy;

        return productLow(hx, lx, hy, ly, xy);
    }

    /**
     * Compute the low part of the double length number {@code (z,zz)} for the exact
     * product of {@code x} and {@code y} using Dekker's mult12 algorithm. The standard
     * precision product {@code x*y} must be provided. The numbers {@code x} and {@code y}
     * should already be split into low and high parts.
     *
     * <p>Note: This uses the high part of the result {@code (z,zz)} as {@code x * y} and not
     * {@code hx * hy + hx * ty + tx * hy} as specified in Dekker's original paper.
     * See Shewchuk (1997) for working examples.
     *
     * @param hx High part of first factor.
     * @param lx Low part of first factor.
     * @param hy High part of second factor.
     * @param ly Low part of second factor.
     * @param xy Product of the factors.
     * @return <code>lx * ly - (((xy - hx * hy) - lx * hy) - hx * ly)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    static double productLow(double hx, double lx, double hy, double ly, double xy) {
        // Compute the multiply low part:
        // err1 = xy - hx * hy
        // err2 = err1 - lx * hy
        // err3 = err2 - hx * ly
        // low = lx * ly - err3
        return lx * ly - (((xy - hx * hy) - lx * hy) - hx * ly);
    }

    /**
     * Compute the round-off {@code s} from the sum of two split numbers {@code (x, xx)}
     * and {@code (y, yy)} using Dekker's add2 algorithm. The values are not required to be
     * ordered by magnitude as an absolute comparison is made to determine the summation order.
     * The sum of the high parts {@code r} must be provided.
     *
     * <p>The result {@code (r, s)} must be re-balanced to create the split result {@code (z, zz)}:
     * <pre>
     * z = r + s
     * zz = r - z + s
     * </pre>
     *
     * @param x High part of first number.
     * @param xx Low part of first number.
     * @param y High part of second number.
     * @param yy Low part of second number.
     * @param r Sum of the parts (x + y) = r
     * @return The round-off from the sum (x + y) = s
     */
    static double sumLow(double x, double xx, double y, double yy, double r) {
        return Math.abs(x) > Math.abs(y) ?
                x - r + y + yy + xx :
                y - r + x + xx + yy;
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Dekker's two-sum algorithm. The values are required to be ordered by magnitude
     * {@code |a| >= |b|}. The standard precision sum must be provided.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param sum Sum of the parts (a + b).
     * @return <code>b - (sum - a)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 6</a>
     */
    static double fastTwoSumLow(double a, double b, double sum) {
        // bVitual = sum - a
        // b - bVirtual == b round-off
        return b - (sum - a);
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Knuth's two-sum algorithm. The values are not required to be ordered by magnitude.
     * The standard precision sum must be provided.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param sum Sum of the parts (a + b).
     * @return <code>(b - (sum - (sum - b))) + (a - (sum - b))</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 7</a>
     */
    static double twoSumLow(double a, double b, double sum) {
        final double bVirtual = sum - a;
        // sum - bVirtual == aVirtual.
        // a - aVirtual == a round-off
        // b - bVirtual == b round-off
        return (a - (sum - bVirtual)) + (b - bVirtual);
    }
}
