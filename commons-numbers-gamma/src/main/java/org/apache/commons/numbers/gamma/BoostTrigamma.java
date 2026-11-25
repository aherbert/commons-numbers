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
 * <a href="https://en.wikipedia.org/wiki/Trigamma_function">Trigamma function</a>.
 *
 * <p>This code has been adapted from the <a href="https://www.boost.org/">Boost</a>
 * {@code c++} implementation {@code <boost/math/special_functions/trigamma.hpp>}.
 * All work is copyright to the original authors and subject to the Boost Software License.
 *
 * @see
 * <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/sf_gamma/trigamma.html">
 * Boost trigamma Function</a>
 * @since 1.4
 */
final class BoostTrigamma {
    /** 1/2. */
    private static final double HALF = 0.5;
    /** Threshold for small |x| where sin(pi * x) = pi * x. */
    private static final double SMALL_PI_X = 1e-9;
    /** pi^2 is correctly rounded (from pi with 50 digits) using Math.PI.
     * It has a relative error of 0.0357 ulp. */
    private static final double PI_SQUARED = Math.PI * Math.PI;

    /** No instances. */
    private BoostTrigamma() {}

    /**
     * Returns the trigamma function.
     *
     * @param x Argument.
     * @return trigamma(x)
     */
    static double trigamma(double x) {
        if (!Double.isFinite(x)) {
            // Function oscillates between infinite poles at integer values for
            // negative x; the -infinity limit is undefined
            return x == Double.POSITIVE_INFINITY ? 0 : Double.NaN;
        }
        return trigammaImp(x);
    }

    /**
     * Implementation for the trigamma function.
     * Adapted from {@code boost::math::detail::trigamma_imp}.
     *
     * @param x Argument.
     * @return trigamma(x)
     */
    private static double trigammaImp(double x) {
        //
        // Check for negative arguments and use reflection:
        //
        if (x <= 0) {
            final double cx = Math.ceil(x);
            if (cx == x) {
                // policies::raise_pole_error
                // Note that approach to the pole from above or below integer x
                // tend to infinity.
                return Double.POSITIVE_INFINITY;
            }
            // Reflect:
            // pi^2 / sin(pi x)^2  - trigamma(1 - x)
            return pi2sin2pi(x, cx) - trigammaPrec53(1 - x);
        }
        double result = 0;
        if (x < 1) {
            result = 1 / (x * x);
            x += 1;
        }
        return result + trigammaPrec53(x);
    }

    /**
     * Implementation for the trigamma function with {@code 1 <= x}.
     * Adapted from {@code boost::math::detail::trigamma_prec} for
     * 53-bit precision.
     *
     * @param x Argument in [1, 2]
     * @return trigamma(x)
     */
    private static double trigammaPrec53(double x) {
        if (x <= 2) {
            // Max error in interpolated form: 3.736e-017
            double p;
            p = 0.64468386819102836;
            p =  1.6638069578676164 + p * x;
            p = 0.28080574467981213 + p * x;
            p = -3.3703848401898283 + p * x;
            p = -3.8310674472619321 + p * x;
            p = -1.1093280605946045 + p * x;
            double q;
            q = -0.20314516859987728e-6;
            q =     0.64468798399785611 + q * x;
            q =      2.7012734178351534 + q * x;
            q =      4.5208926987851437 + q * x;
            q =      3.4535389668541151 + q * x;
            q =                     1.0 + q * x;
            return (2.1093254089355469 + p / q) / (x * x);
        }
        final double y = 1 / x;
        if (x <= 4) {
            // Max error in interpolated form: 1.159e-017
            double p;
            p =     0.74566981111565923;
            p =      2.0534873203680393 + p * y;
            p =      2.5645435828098254 + p * y;
            p =      1.6077979838469348 + p * y;
            p =     0.50000049158540261 + p * y;
            p = -0.13803835004508849e-7 + p * y;
            double q;
            q = -0.00057069112416246805;
            q =     0.74967671848044792 + q * y;
            q =      2.7853527819234466 + q * y;
            q =      4.1681660554090917 + q * y;
            q =      2.8822787662376169 + q * y;
            q =                     1.0 + q * y;
            return (1 + p / q) / x;
        }
        // Maximum Deviation Found:                     6.896e-018
        // Expected Error Term :                       -6.895e-018
        // Maximum Relative Change in Control Points :  8.497e-004
        double p;
        p =     0.40154388356961734;
        p =      1.5897035272532764 + p * y;
        p =      2.1921221359427595 + p * y;
        p =       2.498208511343429 + p * y;
        p =      1.0177274392923795 + p * y;
        p =     0.49999999999998975 + p * y;
        p = 0.68947581948701249e-17 + p * y;
        double q;
        q = 0.022892987908906897;
        q =  0.28360399799075752 + q * y;
        q =   2.3013614809773616 + q * y;
        q =   2.9745631894384922 + q * y;
        q =   4.4290431747556469 + q * y;
        q =   1.7021215452463932 + q * y;
        q =                  1.0 + q * y;
        return (1 + p / q) / x;
    }

    /**
     * Ad hoc function calculates pi^2 / sin(pi * x)^2, taking extra care near when x is
     * near a whole number.
     *
     * @param x Value (assumed to be negative)
     * @param cx Ceiling x.
     * @return pi^2 / sin(pi * x)^2
     */
    private static double pi2sin2pi(double x, double cx) {
        // This function exploits the periodic sine function.
        // Sine is periodic in [0, 2pi). However the magnitude of x
        // which sets sign of sin(pi x) if ceil(x) is odd or even is not
        // required as the result is squared. So this function only
        // requires the half-period result in [0, pi).

        // x < 0 : z in [0, 1)
        final double z = cx - x;
        // Note: sin(pi z) ~ pi z when z is small
        // pi^2 / (pi z)^2 = 1 / z^2
        if (z < SMALL_PI_X) {
            // This result will reduce the rounding error as x -> integer
            // and the trigamma function is a pole: it avoids squaring
            // the rounding error of pi * z, and using the inexact pi^2 constant.
            return 1 / (z * z);
        }
        // sin is symmetric in (0, pi/2) and (pi/2, pi).
        final double s = Math.sin(Math.PI * (z < HALF ? z : 1 - z));
        return PI_SQUARED / (s * s);
    }
}
