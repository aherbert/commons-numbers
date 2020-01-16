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
package org.apache.commons.numbers.arrays;

/**
 * Computes linear combinations accurately.
 * This class computes the sum of the products of two sequences of numbers
 * <code>a<sub>i</sub> b<sub>i</sub></code> to high accuracy.
 * It does so by using extended precision multiplication and addition algorithms to
 * preserve accuracy and reduce cancellation effects.
 *
 * <p>It is based on the 2005 paper
 * <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
 * Accurate Sum and Dot Product</a> by Takeshi Ogita, Siegfried M. Rump,
 * and Shin'ichi Oishi published in <em>SIAM J. Sci. Comput</em>.
 *
 * <p>Due to the methods used to increase precision it is possible that the computation
 * overflows when the standard dot product does not. This may occur when an individual
 * product exceeds the magnitude of {@link Double#MAX_VALUE}{@code  / 2}.
 * In this case the standard precision result will be returned.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Dot_product">Dot product</a>
 */
public final class LinearCombination {
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
     * This process is time consuming and the algorithms here by Ogita et al (2005) avoid
     * a full sort and perform a single rebalancing pass on the parts without sort;
     * however the parts are at least partially ordered by the algorithm processing which
     * maintains the order of the total sum round-off parts. The result is not expected to
     * be exact for all input but will out-perform a standard precision sum of products.
     * A single pass outputs the sum in 3-fold precision. Increasing the number of
     * rebalancing passes from 1 to (K-2) would increase to K-fold precision.
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
     * limit is a value with an exponent of (1023 - 28) = 2^995. */
    private static final double SAFE_UPPER = 0x1.0p995;

    /** The scale to use when down-scaling during a split into a high part.
     * This must be larger than the multiplier and a power of 2 for exact scaling. */
    private static final double DOWN_SCALE = 0x1.0p-30;

    /** The scale to use when re-scaling during a split into a high part.
     * This is the inverse of {@link #DOWN_SCALE}. */
    private static final double UP_SCALE = 0x1.0p30;

    /** Private constructor. */
    private LinearCombination() {
        // intentionally empty.
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a Factors.
     * @param b Factors.
     * @return \( \sum_i a_i b_i \).
     * @throws IllegalArgumentException if the sizes of the arrays are different.
     */
    public static double value(double[] a,
                               double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("Dimension mismatch: " + a.length + " != " + b.length);
        }

        final int len = a.length;

        if (len == 0) {
            return 0;
        }
        if (len == 1) {
            // Revert to scalar multiplication.
            return a[0] * b[0];
        }

        // Implement dotK (Algorithm 5.10) from Ogita et al (2005).
        // Store all round-off parts.
        // Round-off parts of each product are r[0 to (n-1)].
        // Round-off parts of each sum are r[n to (2n-2)].
        // The standard precision scalar product is term p which becomes r[2n-1].
        final double[] r = new double[len * 2];

        // p is the standard scalar product sum initialised with the first product
        double p = a[0] * b[0];
        r[0] = productLow(a[0], b[0], p);

        // Remaining split products added to the current sum and round-off stored
        for (int i = 1; i < len; i++) {
            final double h = a[i] * b[i];
            r[i] = productLow(a[i], b[i], h);

            final double x = p + h;
            r[i + len - 1] = sumLow(p, h, x);
            p = x;
        }

        // Sum the round-off with the standard sum as the final component.
        // Here the value passed to sumK is (K-1) for K-fold precision of the sum (dotK).
        // This method implements 3-fold precision (dot3). Increasing K may be of benefit
        // for longer arrays. The iteration is an error-free transform and higher K
        // never loses precision.
        r[r.length - 1] = p;
        return getSum(p, sumK(r, 2));
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @return \( a_1 b_1 + a_2 b_2 \)
     *
     * @see #value(double, double, double, double, double, double)
     * @see #value(double, double, double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2) {
        // Round-off parts of each product are r[0-1].
        // Round-off parts of each sum are r[2].
        // The standard precision scalar product is term p which becomes r3.
        // Working variables h (current product high part) and x (new sum).
        double h;
        double x;
        double p = a1 * b1;
        double r0 = productLow(a1, b1, p);
        h = a2 * b2;
        double r1 = productLow(a2, b2, h);
        x = p + h;
        double r2 = sumLow(p, h, x);
        p = x;

        double r3 = x;

        // In-line 1 round of vector sum for 3-fold precision (dot3).
        x = r0 + r1;
        r0 = sumLow(r1, r0, x);
        r1 = x;
        x = r1 + r2;
        r1 = sumLow(r2, r1, x);
        r2 = x;
        x = r2 + r3;
        r2 = sumLow(r3, r2, x);
        r3 = x;

        // Final summation
        return getSum(p, r0 + r1 + r2 + r3);
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param a3 First factor of the third term.
     * @param b3 Second factor of the third term.
     * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 \)
     *
     * @see #value(double, double, double, double)
     * @see #value(double, double, double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2,
                               double a3, double b3) {
        // Round-off parts of each product are r[0-2].
        // Round-off parts of each sum are r[3-4].
        // The standard precision scalar product is term p which becomes r5.
        // Working variables h (current product high part) and x (new sum).
        double h;
        double x;
        double p = a1 * b1;
        double r0 = productLow(a1, b1, p);
        h = a2 * b2;
        double r1 = productLow(a2, b2, h);
        x = p + h;
        double r3 = sumLow(p, h, x);
        p = x;
        h = a3 * b3;
        double r2 = productLow(a3, b3, h);
        x = p + h;
        double r4 = sumLow(p, h, x);
        p = x;

        double r5 = x;

        // In-line 1 round of vector sum for 3-fold precision (dot3).
        x = r0 + r1;
        r0 = sumLow(r1, r0, x);
        r1 = x;
        x = r1 + r2;
        r1 = sumLow(r2, r1, x);
        r2 = x;
        x = r2 + r3;
        r2 = sumLow(r3, r2, x);
        r3 = x;
        x = r3 + r4;
        r3 = sumLow(r4, r3, x);
        r4 = x;
        x = r4 + r5;
        r4 = sumLow(r5, r4, x);
        r5 = x;

        // Final summation
        return getSum(p, r0 + r1 + r2 + r3 + r4 + r5);
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param a3 First factor of the third term.
     * @param b3 Second factor of the third term.
     * @param a4 First factor of the fourth term.
     * @param b4 Second factor of the fourth term.
     * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 + a_4 b_4 \)
     *
     * @see #value(double, double, double, double)
     * @see #value(double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2,
                               double a3, double b3,
                               double a4, double b4) {
        // Round-off parts of each product are r[0-3].
        // Round-off parts of each sum are r[4-6].
        // The standard precision scalar product is term p which becomes r7.
        // Working variables h (current product high part) and x (new sum).
        double h;
        double x;
        double p = a1 * b1;
        double r0 = productLow(a1, b1, p);
        h = a2 * b2;
        double r1 = productLow(a2, b2, h);
        x = p + h;
        double r4 = sumLow(p, h, x);
        p = x;
        h = a3 * b3;
        double r2 = productLow(a3, b3, h);
        x = p + h;
        double r5 = sumLow(p, h, x);
        p = x;
        h = a4 * b4;
        double r3 = productLow(a4, b4, h);
        x = p + h;
        double r6 = sumLow(p, h, x);
        p = x;

        double r7 = x;

        // In-line 1 round of vector sum for 3-fold precision (dot3).
        x = r0 + r1;
        r0 = sumLow(r1, r0, x);
        r1 = x;
        x = r1 + r2;
        r1 = sumLow(r2, r1, x);
        r2 = x;
        x = r2 + r3;
        r2 = sumLow(r3, r2, x);
        r3 = x;
        x = r3 + r4;
        r3 = sumLow(r4, r3, x);
        r4 = x;
        x = r4 + r5;
        r4 = sumLow(r5, r4, x);
        r5 = x;
        x = r5 + r6;
        r5 = sumLow(r6, r5, x);
        r6 = x;
        x = r6 + r7;
        r6 = sumLow(r7, r6, x);
        r7 = x;

        // Final summation
        return getSum(p, r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7);
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
     * @param value Value.
     * @return the high part of the value.
     * @see <a href="https://doi.org/10.1007/BF01397083">
     * Dekker (1971) A floating-point technique for extending the available precision</a>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 17</a>
     */
    private static double highPart(double value) {
        // Avoid overflow
        if (value >= SAFE_UPPER || value <= -SAFE_UPPER) {
            // Do scaling.
            final double x = value * DOWN_SCALE;
            final double c = MULTIPLIER * x;
            final double hi = (c - (c - x)) * UP_SCALE;
            if (Double.isInfinite(hi)) {
                // Number is too large.
                // This occurs if value is infinite or close to Double.MAX_VALUE.
                // Note that multiplication by (2^s+1) can cause hi to have an exponent
                // 1 greater than input value which will overflow if the exponent is already +1023.
                // Revert to the raw upper 26 bits of the 53-bit mantissa (including the assumed
                // leading 1 bit). This conversion will result in the low part being a
                // 27-bit significand and the potential loss of bits during addition and
                // multiplication. (Contrast to the Dekker split which creates two 26-bit
                // numbers with a bit moved to the sign of low.)
                // The conversion will maintain Infinite in the high part where the resulting
                // low part (value - high) is NaN.
                return Double.longBitsToDouble(Double.doubleToRawLongBits(value) & ((-1L) << 27));
            }
            return hi;
        }
        // normal conversion
        final double c = MULTIPLIER * value;
        return c - (c - value);
    }

    /**
     * Compute the low part of the double length number {@code (z,zz)} for the exact
     * product of {@code x} and {@code y}. Each factor is split into high and low
     * parts to compute the lower part of the product result. The standard precision
     * product must be provided.
     *
     * <p>Note: This uses the high part as {@code x * y} and not
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
     * @see <a href="https://doi.org/10.1007/BF01397083">
     * Dekker (1971) A floating-point technique for extending the available precision</a>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    private static double productLow(double x, double y, double xy) {
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
    private static double sumLow(double a, double b, double sum) {
        final double bVirtual = sum - a;
        // sum - bVirtual == aVirtual.
        // a - aVirtual == a round-off
        // b - bVirtual == b round-off
        return (a - (sum - bVirtual)) + (b - bVirtual);
    }

    /**
     * Sum to K-fold precision.
     *
     * @param p Data to sum.
     * @param k The precision.
     * @return the sum
     */
    private static double sumK(double[] p, int k) {
        // k=1 will skip the vector transformation and sum in standard precision.
        for (int i = 1; i < k; i++) {
            vecSum(p);
        }
        double sum = 0;
        for (final double pi : p) {
            sum += pi;
        }
        return sum;
    }

    /**
     * Error free vector transformation for summation.
     *
     * @param p Data.
     */
    private static void vecSum(double[] p) {
        for (int i = 1; i < p.length; i++) {
            final double x = p[i] + p[i - 1];
            p[i - 1] = sumLow(p[i], p[i - 1], x);
            p[i] = x;
        }
    }

    /**
     * Gets the final sum. This checks the high precision sum is finite, otherwise
     * returns the standard precision sum for the IEEE754 result.
     *
     * <p>The standard of high precision sum may be non-finite due to input infinite
     * or NaN numbers or overflow in the summation. However the high precision sum
     * can also be non-finite when the standard dot sum is finite. This occurs when
     * the split product had a component high-part that overflowed during
     * computation of the hx * hy partial result. In all cases returning the
     * standard dot sum ensures the IEEE754 result.
     *
     * @param sum Standard dot sum.
     * @param dotKSum High precision dotK sum.
     * @return the sum
     */
    private static double getSum(double sum, double dotKSum) {
        if (!Double.isFinite(dotKSum)) {
            // Either we have split infinite numbers, some coefficients were NaNs,
            // or the split product overflowed.
            // Return the naive implementation for the IEEE754 result
            return sum;
        }
        return dotKSum;
    }
}
