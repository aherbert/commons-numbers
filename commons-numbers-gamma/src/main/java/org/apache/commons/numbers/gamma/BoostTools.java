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

//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

package org.apache.commons.numbers.gamma;

import java.util.function.DoubleSupplier;
import java.util.function.Supplier;
import org.apache.commons.numbers.core.Precision;

/**
 * Utility tools used by the Boost functions.
 *
 * <p>This code has been adapted from the <a href="https://www.boost.org/">Boost</a>
 * {@code c++} implementations in {@code <boost/math/tools/>}.
 * All work is copyright John Maddock 2006 and subject to the Boost Software License.
 */
public final class BoostTools {
    /**
     * The value for any number close to zero.
     *
     * <p>"The parameter small should be some non-zero number less than typical values of
     * eps * |b_n|, e.g., 1e-50". The Boost continued fraction evaluation uses
     * 16 * min normal = 3.5601181736115222e-307.
     */
    private static final double SMALL = 16 * Double.MIN_NORMAL;
    /**
     * The minimum epsilon value for relative error.
     * Equals to Math.ulp(1.0) or 2^-52.
     */
    private static final double EPSILON = 0x1.0p-52;
    /** Message for failure to converge. */
    private static final String MSG_FAILED_TO_CONVERGE = "Failed to converge within %d iterations";

    /** Private constructor. */
    private BoostTools() {
        // intentionally empty.
    }

    /**
     * Evaluate the continued fraction.
     * Evaluates:
     * <pre>
     *            a1
     *      ---------------
     *      b1 +     a2
     *           ----------
     *           b2 +   a3
     *                -----
     *                b3 + ...
     * </pre>
     * <p>Note that the first a1 and b1 returned by generator cf are both used.
     *
     * <p>Adapted from {@code boost/math/tools/fraction.hpp continued_fraction_a}.
     * This differs from {@code continued_fraction_b} which adds term b0 to the fraction
     * (and discards term a0). That implementation usage matches the usage of
     * {@link org.apache.commons.numbers.fraction.ContinuedFraction}.
     *
     * @param cf Continued fraction generator
     * @param eps Maximum error allowed
     * @param maxTerms Maximum number of terms
     * @return result
     */
    static double continuedFractionA(Supplier<double[]> cf, double eps, int maxTerms) {
        final double terminator = Math.abs(eps);

        double[] v = cf.get();
        double f = updateIfCloseToZero(v[1]);
        final double a0 = v[0];
        double c = f;
        double d = 0;

        int counter = maxTerms;

        double delta;
        do {
            v = cf.get();
            final double a = v[0];
            final double b = v[1];
            d = updateIfCloseToZero(b + a * d);
            c = updateIfCloseToZero(b + a / c);
            d = 1 / d;
            delta = c * d;
            f = f * delta;
        } while ((Math.abs(delta - 1) > terminator) && --counter > 0);

        if (counter <= 0) {
            throw new ArithmeticException(
               String.format(MSG_FAILED_TO_CONVERGE, maxTerms));
        }

        return a0 / f;
    }

    /**
     * Returns the value, or if close to zero returns a small epsilon.
     *
     * <p>This method is used in Thompson & Barnett to monitor both the numerator and denominator
     * ratios for approaches to zero.
     *
     * @param value the value
     * @return the value (or small epsilon)
     */
    private static double updateIfCloseToZero(double value) {
        return Precision.equals(value, 0.0, SMALL) ? SMALL : value;
    }

    /**
     * Sum the series.
     *
     * <p>Adapted from {@code boost/math/tools/series.hpp}.
     *
     * @param func Series generator
     * @param epsilon Maximum relative error allowed
     * @param maxTerms Maximum number of terms
     * @return result
     */
    static double sumSeries(DoubleSupplier func, double epsilon, int maxTerms) {
        return sumSeries(func, epsilon, maxTerms, 0);
    }

    /**
     * Sum the series.
     *
     * <p>Adapted from {@code boost/math/tools/series.hpp}.
     *
     * @param func Series generator
     * @param epsilon Maximum relative error allowed
     * @param maxTerms Maximum number of terms
     * @param initValue Initial value
     * @return result
     */
    static double sumSeries(DoubleSupplier func, double epsilon, int maxTerms, double initValue) {
        // Note:
        // The Boost code requires eps to be non-zero. It is created in the
        // <boost/math/policies/policy.hpp> as a non-zero relative error term.
        // An alternative termination condition with a divide is:
        // (eps < Math.abs(nextTerm / result))
        //
        // Here the argument is checked against the minimum epsilon for a double
        // to provide functional equivalence with the Boost policy.
        // In the min eps case the loop terminates if the most recently added term is
        // 0 or 1 ulp of the result. This condition is acceptable if the next
        // computed term will be at most half of the most recent term (thus
        // cannot be added to the current result).

        final double eps = getEpsilon(epsilon);

        int counter = maxTerms;

        double result = initValue;
        double nextTerm;
        do {
            nextTerm = func.getAsDouble();
            result += nextTerm;
        } while ((Math.abs(eps * result) < Math.abs(nextTerm)) && --counter > 0);

        if (counter <= 0) {
            throw new ArithmeticException(
               String.format(MSG_FAILED_TO_CONVERGE, maxTerms));
        }

        return result;
    }


    /**
     * Sum the series using Kahan summation.
     *
     * <p>Adapted from {@code boost/math/tools/series.hpp}.
     *
     * @param func Series generator
     * @param epsilon Maximum relative error allowed
     * @param maxTerms Maximum number of terms
     * @return result
     */
    static double kahanSumSeries(DoubleSupplier func, double epsilon, int maxTerms) {
        return kahanSumSeries(func, epsilon, maxTerms, 0);
    }

    /**
     * Sum the series using Kahan summation.
     *
     * <p>Adapted from {@code boost/math/tools/series.hpp}.
     *
     * @param func Series generator
     * @param epsilon Maximum relative error allowed
     * @param maxTerms Maximum number of terms
     * @param initValue Initial value
     * @return result
     */
    static double kahanSumSeries(DoubleSupplier func, double epsilon, int maxTerms, double initValue) {
        final double eps = getEpsilon(epsilon);

        int counter = maxTerms;

        // Kahan summation:
        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        // This summation is accurate if the term is smaller in magnitude
        // than the current sum. This is a condition required for the
        // series termination thus the extended precision sum need not
        // check magnitudes of terms to compute the carry.

        double result = initValue;
        double carry = 0;
        double nextTerm;
        do {
            nextTerm = func.getAsDouble();
            final double y = nextTerm - carry;
            final double t = result + y;
            carry = t - result;
            carry -= y;
            result = t;
        } while ((Math.abs(eps * result) < Math.abs(nextTerm)) && --counter > 0);

        if (counter <= 0) {
            throw new ArithmeticException(
               String.format(MSG_FAILED_TO_CONVERGE, maxTerms));
        }

        return result;
    }

    /**
     * Gets the epsilon ensuring it satisfies the minimum value for the relative error
     * of a {@code double} type.
     *
     * @param epsilon Configured epsilon
     * @return the epsilon
     */
    private static double getEpsilon(double epsilon) {
        return epsilon > EPSILON ? epsilon : EPSILON;
    }

    /**
     * Evaluate the polynomial using Horner's method.
     *
     * @param c Polynomial coefficients (must have length > 0)
     * @param x Argument x
     * @return polynomial value
     */
    static double evaluatePolynomial(double[] c, double x) {
        final int count = c.length;
        double sum = c[count - 1];
        for (int i = count - 2; i >= 0; --i) {
            sum *= x;
            sum += c[i];
        }
        return sum;
    }
}
