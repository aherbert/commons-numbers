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

//  Copyright John Maddock 2006.
//  Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

package org.apache.commons.numbers.gamma;

/**
 * Implementation of the
 * <a href="https://en.wikipedia.org/wiki/Digamma_function">Digamma function</a>.
 *
 * <p>This code has been adapted from the <a href="https://www.boost.org/">Boost</a>
 * {@code c++} implementation {@code <boost/math/special_functions/digamma.hpp>}.
 * All work is copyright to the original authors and subject to the Boost Software License.
 *
 * @see
 * <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/sf_gamma/digamma.html">
 * Boost Digamma Function</a>
 * @since 1.4
 */
final class BoostDigamma {
    /** 1/2. */
    private static final double HALF = 0.5;
    /** Limit for asymptotic expansion. */
    private static final double DIGAMMA_LARGE_LIMIT = 10;
    /** Constant for the digamma rational function. */
    private static final double Y = 0.99558162689208984;
    /** Constant for the digamma rational function. Value = 1.4616321446374059.
     * This is close to the root of the digamma function. 1.461632144968362341262660 (25 digits).
     * Note: ROOT1 + ROOT2 + ROOT3 ~ root of digamma in extended precision:
     * 1.46163214496836234126265954232572132846820517356040000000000 (added as BigDecimal). */
    private static final double ROOT1 = 1569415565.0 / 1073741824;
    /** Constant for the digamma rational function. Value = 3.309564688275257e-10. */
    private static final double ROOT2 = (381566830.0 / 1073741824) / 1073741824;
    /** Constant for the digamma rational function. */
    private static final double ROOT3 = 0.9016312093258695918615325266959189453125e-19;

    /** No instances. */
    private BoostDigamma() {}

    /**
     * Returns the digamma function.
     *
     * @param x Argument.
     * @return digamma(x)
     */
    static double digamma(double x) {
        if (!Double.isFinite(x)) {
            // Function oscillates between infinite poles at integer values for
            // negative x; the -infinity limit is undefined
            return x == Double.POSITIVE_INFINITY ? x : Double.NaN;
        }
        return digammaImp(x);
    }

    /**
     * Implementation for the digamma function.
     * Adapted from {@code boost::math::detail::digamma_imp}.
     *
     * @param x Argument.
     * @return digamma(x)
     */
    private static double digammaImp(double x) {
        double result = 0;
        //
        // Check for negative arguments and use reflection:
        //
        if (x <= -1) {
            // Reflect:
            x = 1 - x;
            // Argument reduction for tan:
            double remainder = x - Math.floor(x);
            // Shift to negative if > 0.5:
            if (remainder > HALF) {
                remainder -= 1;
            } else if (remainder == 0) {
                // Digamma is undefined at negative integers.
                // As x -> -n: digamma(x) -> +infty from one side and -infty from the other.
                return Double.NaN;
            }
            result = Math.PI / Math.tan(Math.PI * remainder);
        } else if (x == 0) {
            // Undefined at the pole:
            // As -x -> 0: digamma(x) -> +infty
            // As x -> 0: digamma(x) -> -infty
            return Double.NaN;
        }
        //
        // If we're above the lower-limit for the
        // asymptotic expansion then use it:
        //
        if (x >= DIGAMMA_LARGE_LIMIT) {
            result += digammaImpLarge(x);
        } else {
            //
            // If x > 2 reduce to the interval [1,2]:
            //
            while (x > 2) {
                x -= 1;
                result += 1 / x;
            }
            //
            // If x < 1 use recurrence to shift to > 1:
            //
            while (x < 1) {
                result -= 1 / x;
                x += 1;
            }
            result += digammaImp1to2(x);
        }
        return result;
    }

    /**
     * Implementation for the digamma function using the asymptotic expansion.
     * Adapted from {@code boost::math::detail::digamma::digamma_imp_large} for
     * 53-bit precision.
     *
     * @param x Argument {@code >= 10}
     * @return digamma(x)
     */
    private static double digammaImpLarge(double x) {
        x -= 1;
        double result = Math.log(x);
        // Change from boost is to subtract the smaller terms first:
        // original: result += 1 / (2 * x)
        // boost::tools::evaluate_polynomial
        final double z = 1 / (x * x);
        double p;
        p =   -0.4432598039215686;
        p =   0.08333333333333333 + p * z;
        p = -0.021092796092796094 + p * z;
        p =  0.007575757575757576 + p * z;
        p = -0.004166666666666667 + p * z;
        p =  0.003968253968253968 + p * z;
        p = -0.008333333333333333 + p * z;
        p =   0.08333333333333333 + p * z;
        // Combine smaller terms first, then add log(x - 1)
        // original: result -= z * p
        result += 1 / (2 * x) - z * p;
        return result;
    }

    /**
     * Implementation for the digamma function with {@code 1 <= x <= 2}.
     * Adapted from {@code boost::math::detail::digamma_imp_1_2} for
     * 53-bit precision.
     *
     * @param x Argument in [1, 2]
     * @return digamma(x)
     */
    private static double digammaImp1to2(double x) {
        // Root 1 - 3 are successively smaller. Subtract sequentially.
        double g = x - ROOT1;
        g -= ROOT2;
        g -= ROOT3;
        final double xm1 = x - 1;
        // boost::tools::evaluate_polynomial P and Q
        double p;
        p = -0.0020713321167745952;
        p =  -0.045251321448739056 + p * xm1;
        p =   -0.28919126444774784 + p * xm1;
        p =   -0.65031853770896507 + p * xm1;
        p =   -0.32555031186804491 + p * xm1;
        p =    0.25479851061131551 + p * xm1;
        double q;
        q = -0.55789841321675513e-6;
        q =   0.0021284987017821144 + q * xm1;
        q =    0.054151797245674225 + q * xm1;
        q =     0.43593529692665969 + q * xm1;
        q =      1.4606242909763515 + q * xm1;
        q =      2.0767117023730469 + q * xm1;
        q =                     1.0 + q * xm1;
        final double r = p / q;
        return g * Y + g * r;
    }
}
