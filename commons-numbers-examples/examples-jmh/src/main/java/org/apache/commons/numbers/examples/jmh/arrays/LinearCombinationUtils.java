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

import org.apache.commons.rng.UniformRandomProvider;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Utility class to create data for linear combinations.
 */
final class LinearCombinationUtils {
    /** ln(2). */
    private static final double LN_2 = Math.log(2);

    /** No construction. */
    private LinearCombinationUtils() {}
    /**
     * Generates ill conditioned dot products.
     * See {@link #genDot(double, UniformRandomProvider, double[], double[], double[], MathContext)}.
     *
     * <p>The exact dot product is computed using {@link MathContext#UNLIMITED}. This
     * will not scale to large lengths.
     *
     * @param c anticipated condition number of x’*y
     * @param rng source of randomness
     * @param x output array vector (length at least 6)
     * @param y output array vector (length at least 6)
     * @param computeC if not null compute the condition number and place at index 0
     * @return the exact dot product
     * @throws IllegalArgumentException If the vector length is below 6
     */
    static double genDot(double c, UniformRandomProvider rng,
            double[] x, double[] y, double[] computeC) {
        return genDot(c, rng, x, y, computeC, MathContext.UNLIMITED);
    }

    /**
     * Generates ill conditioned dot products. The length of the vectors should be
     * {@code n>=6}.
     *
     * <p>The condition number is defined as twice the sum of the absolute products
     * divided by the absolute of the dot product:
     *
     * <pre>
     * d = x'*y = DotExact(x, y)
     * C = 2*(abs(x’)*abs(y))/abs(d)
     * </pre>
     *
     * <p>A high condition number means that large values are added together to make
     * a small value, i.e. large cancellation occurs. Computation of the actual dot
     * product requires an exact method to avoid cancellation floating-point errors.
     * This method uses BigDecimal to compute DotExact(x, y). The dot product result
     * d is created in the range [-1, 1]. The actual condition number can be
     * obtained by passing an array {@code computeC}. This is computed in standard
     * precision thus if finite then the dot product of the factors does not
     * overflow. The routine has been tested with anticipated condition numbers up
     * to 1e300.
     *
     * <p>Uses the GenDot algorithm 6.1 from <a
     * href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.2.1547">
     * Accurate Sum and Dot Product</a> by Takeshi Ogita, Siegfried M. Rump, and
     * Shin'ichi Oishi published in <em>SIAM J. Sci. Comput</em>.
     *
     * <p>Note: Ogita et al state "in general the actual condition number is a
     * little larger than the anticipated one with quite some variation". This
     * implementation is computing a condition number approximately equal to the
     * requested for the parameters that have been tested.
     *
     * @param c anticipated condition number of x’*y
     * @param rng source of randomness
     * @param x output array vector (length at least 6)
     * @param y output array vector (length at least 6)
     * @param computeC if not null compute the condition number and place at index 0
     * @param context the context to use for rounding the exact sum
     * @return the exact dot product
     * @throws IllegalArgumentException If the vector length is below 6
     */
    static double genDot(double c, UniformRandomProvider rng,
            double[] x, double[] y, double[] computeC, MathContext context) {
        // Initialisation
        if (x.length < 6) {
            throw new IllegalArgumentException("Incorrect usage: n < 6: " + x.length);
        }
        final int n = x.length;
        final int n2 = n / 2;

        // log2(c)
        final double b = Math.log(c) / LN_2;
        final double b2 = b / 2;
        // e vector of exponents between 0 and b/2
        int[] e = new int[n2];
        // make sure exponents b/2 and 0 actually occur in e
        e[0] = (int) Math.round(b2) + 1;
        for (int i = 1; i < n2 - 1; i++) {
            e[i] = (int) Math.round(rng.nextDouble() * b2);
        }
        // e[end] = 0;

        // Generate first half vectors.
        // Maintain the exact dot product for use later
        BigDecimal exact = BigDecimal.ZERO;
        for (int i = 0; i < n2; i++) {
            x[i] = Math.scalb(m1p1(rng), e[i]);
            y[i] = Math.scalb(m1p1(rng), e[i]);
            exact = exact.add(new BigDecimal(x[i]).multiply(new BigDecimal(y[i])), context);
        }

        // for i=n2+1:n and v=1:i,
        // generate x(i), y(i) such that (*) x(v)’*y(v) ~ 2^e(i-n2)
        // i.e. the dot product up to position i is a value that will increasingly approach 0
        e = new int[n - n2];
        // exponents for second half as a linear vector of exponents from b/2 to 0
        // linspace(b/2, 0, n-n2)
        for (int i = 0; i < e.length - 1; i++) {
            e[i] = (int) Math.round(b2 * (e.length - i - 1) / (e.length - 1));
        }

        for (int i = n2; i < x.length; i++) {
            // x(i) random with generated exponent
            x[i] = Math.scalb(m1p1(rng), e[i - n2]);
            // y(i) according to (*).
            // sum(i) = xi * yi + sum(i-1)
            // yi = (sum(i) - sum(i-1)) / xi
            // Here the new sum(i) is a random number with the
            y[i] = (Math.scalb(m1p1(rng), e[i - n2]) - exact.doubleValue()) / x[i];
            // Maintain the exact dot product
            exact = exact.add(new BigDecimal(x[i]).multiply(new BigDecimal(y[i])), context);
        }

        // Shuffle x and y. Do a parallel Fisher-Yates shuffle.
        for (int i = n; i > 1; i--) {
            final int j = rng.nextInt(i);
            swap(x, i - 1, j);
            swap(y, i - 1, j);
        }

        // Compute condition number:
        // d = DotExact(x’,y);             % the true dot product rounded to nearest
        // C = 2*(abs(x’)*abs(y))/abs(d);  % the actual condition number
        final double d = exact.doubleValue();
        if (computeC != null) {
            double cc = 0;
            for (int i = 0; i < n; i++) {
                cc += Math.abs(x[i] * y[i]);
            }
            computeC[0] = 2 * cc / Math.abs(d);
        }
        return d;
    }

    /**
     * Create a double in the range [-1, 1).
     *
     * @param rng source of randomness
     * @return the double
     */
    private static double m1p1(UniformRandomProvider rng) {
        // Create in the range [0, 1) then randomly subtract 1.
        // This samples the 2^54 dyadic rationals in the range.
        return rng.nextDouble() - rng.nextInt(1);
    }

    /**
     * Swaps the two specified elements in the array.
     *
     * @param array Array.
     * @param i First index.
     * @param j Second index.
     */
    private static void swap(double[] array, int i, int j) {
        final double tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }
}
