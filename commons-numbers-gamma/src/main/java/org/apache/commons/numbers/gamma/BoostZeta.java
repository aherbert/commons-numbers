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

//  Copyright John Maddock 2007, 2014.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

package org.apache.commons.numbers.gamma;

/**
 * Implementation of the
 * <a href="https://en.wikipedia.org/wiki/Riemann_zeta_function">
 * Riemann zeta</a> function..
 *
 * <p>This code has been adapted from the <a href="https://www.boost.org/">Boost</a>
 * {@code c++} implementation {@code <boost/math/special_functions/gamma.hpp>}.
 * All work is copyright to the original authors and subject to the Boost Software License.
 *
 * @see
 * <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/zetas/zeta.html">
 * Boost C++ Riemann Zeta Function</a>
 */
public final class BoostZeta {
    /** zeta(s) values for odd integer s excusing s=1. Indexed using (s - 3) / 2.
     * Length 26; max s = 105. */
    private static final double[] ZETA_ODD_INTEGER = {
        1.2020569031595942853997381615114500,
        1.0369277551433699263313654864570342,
        1.0083492773819228268397975498497968,
        1.0020083928260822144178527692324121,
        1.0004941886041194645587022825264699,
        1.0001227133475784891467518365263574,
        1.0000305882363070204935517285106451,
        1.0000076371976378997622736002935630,
        1.0000019082127165539389256569577951,
        1.0000004769329867878064631167196044,
        1.0000001192199259653110730677887189,
        1.0000000298035035146522801860637051,
        1.0000000074507117898354294919810042,
        1.0000000018626597235130490064039099,
        1.0000000004656629065033784072989233,
        1.0000000001164155017270051977592974,
        1.0000000000291038504449709968692943,
        1.0000000000072759598350574810145209,
        1.0000000000018189896503070659475848,
        1.0000000000004547473783042154026799,
        1.0000000000001136868407680227849349,
        1.0000000000000284217097688930185546,
        1.0000000000000071054273952108527129,
        1.0000000000000017763568435791203275,
        1.0000000000000004440892103143813364,
        1.0000000000000001110223025141066134,
    };

    /** No instances. */
    private BoostZeta() {}

    /**
     * Returns the Riemann zeta function.
     *
     * @param s Argument.
     * @return zeta(s)
     */
    static double zeta(double s) {
        if (Double.isNaN(s)) {
            return Double.NaN;
        }
        return zetaImp(s, 1 - s);
    }

    /**
     * Implementation for the zeta function.
     * Adapted from {@code boost::math::special_functions::zeta::zeta_imp}.
     *
     * @param s Argument.
     * @param sc Complement argument (1 - s).
     * @return zeta(s)
     */
    private static double zetaImp(double s, double sc) {
        if (sc == 0) {
            return Double.POSITIVE_INFINITY;
        }
        //
        // Trivial case:
        //
        if (s > 53) {
            return 1;
        }
        double result;
        //
        // Start by seeing if we have a simple closed form:
        //
        if (Math.floor(s) == s) {
            // Special handling for small integer s
            int v = (int) s;
            if (v == s) {
                if (v < 0) {
                    if ((v & 1) == 0) {
                        // Negative even integer
                        return 0;
                    }
                    // Negative odd integer
                    // TODO: Get a cache of the bernoulli numbers for 2n
                    // Q. What is the max value?
//                   int n = (-v + 1) / 2;
//                   if(n <= (int)boost::math::max_bernoulli_b2n<double>::value) {
//                      return double((-v & 1) ? -1 : 1) * boost::math::unchecked_bernoulli_b2n<double>(n) / (1 - v);
//                   }
//                } else if ((v & 1) == 0) {
                    // Positive even integer
                    // TODO: Get a cache of the bernoulli numbers for 2n
//                   if(((v / 2) <= (int)boost::math::max_bernoulli_b2n<double>::value) &&
//                    (v <= (int)boost::math::max_factorial<double>::value))
//                      return double(((v / 2 - 1) & 1) ? -1 : 1) * ldexp(double 1, v - 1) *
//                    static_cast<double>(pow(constants::pi<double, Policy>(), double(v))) *
//                         boost::math::unchecked_bernoulli_b2n<double>(v / 2) /
//                    boost::math::unchecked_factorial<double>(v);
                    // This computes the bernoulli number.
                    // Given s < 53 we can use a cache
//                    return double(((v / 2 - 1) & 1) ? -1 : 1) * ldexp(double 1, v - 1) *
//                    static_cast<double>(pow(constants::pi<double, Policy>(), double(v))) *
//                      boost::math::bernoulli_b2n<double>(v / 2) / boost::math::factorial<double>(v, pol);
                } else {
                    // Positive odd integer with s != 1
                    final int i = (v - 3) >>> 1;
                    return i < ZETA_ODD_INTEGER.length ? ZETA_ODD_INTEGER[i] : 1;
                }
            }
        }

        if (Math.abs(s) < BoostGamma.ROOT_EPSILON) {
            result = -0.5f - BoostGamma.LOG_ROOT_TWO_PI * s;
        } else if (s < 0) {
            // Swap: s is now positive; sc = 1 - s
            double tmp = s;
            s = sc;
            sc = tmp;
            if (Math.floor(sc * 0.5) == sc * 0.5) {
                // Negative even integer
                result = 0;
            } else {
                result = 0;
                if (s > BoostGamma.MAX_FACTORIAL) {
                    // This has been simplified from the Boost implementation which will
                    // catch overflow conditions when compiled with an appropriate evaluation
                    // policy and return signed infinity. Java floating-point arithmetic does
                    // not create overflow exceptions and will return infinity.
                    double mult = BoostGamma.sinp(0.5 * sc) * 2 * zetaImp53(s, sc);
                    result = LogGamma.value(s);
                    result -= s * Math.log(2 * Math.PI);
                    // Possible overflow if result > 709
                    result = Math.exp(result);
                    // Possible overflow.
                    // Needs result to be just on the verge of overflow when /s/ is
                    // very close to a half integer.
                    result *= mult;
                } else {
                    result = BoostGamma.sinp(0.5 * sc) *
                        2 * Math.pow(2 * Math.PI, -s) *
                        Gamma.value(s) *
                        zetaImp53(s, sc);
                }
            }
        } else {
            result = zetaImp53(s, sc);
        }
        return result;
    }

    /**
     * Implementation for the zeta function.
     * Adapted from {@code boost::math::special_functions::zeta::zeta_imp_prec} for
     * 53-bit precision.
     *
     * <p>This is package private to allow usage from the {@link HurwitzZeta} function.
     *
     * @param s Argument (assumed to be positive).
     * @param sc Complement argument (1 - s).
     * @return zeta(s)
     */
    static double zetaImp53(double s, double sc) {
        double result;
        if (s < 1) {
            // Rational Approximation
            // Maximum Deviation Found:                     2.020e-18
            // Expected Error Term:                        -2.020e-18
            // Max error found at double precision:         3.994987e-17
            // LCOV_EXCL_START
            double[] P = {
                0.24339294433593750202,
                -0.49092470516353571651,
                0.0557616214776046784287,
                -0.00320912498879085894856,
                0.000451534528645796438704,
                -0.933241270357061460782e-5,
            };
            double[] Q = {
                1,
                -0.279960334310344432495,
                0.0419676223309986037706,
                -0.00413421406552171059003,
                0.00024978985622317935355,
                -0.101855788418564031874e-4,
            };
            // LCOV_EXCL_STOP
            result = evaluatePolynomial(P, sc) / evaluatePolynomial(Q, sc);
            result -= 1.2433929443359375F;
            result += sc;
            result /= sc;
        } else if (s <= 2) {
            // Maximum Deviation Found:                     9.007e-20
            // Expected Error Term:                         9.007e-20
            // LCOV_EXCL_START
            double[] P = {
                0.577215664901532860516,
                0.243210646940107164097,
                0.0417364673988216497593,
                0.00390252087072843288378,
                0.000249606367151877175456,
                0.110108440976732897969e-4,
            };
            double[] Q = {
                1.0,
                0.295201277126631761737,
                0.043460910607305495864,
                0.00434930582085826330659,
                0.000255784226140488490982,
                0.10991819782396112081e-4,
            };
            // LCOV_EXCL_STOP
            result = evaluatePolynomial(P, -sc) / evaluatePolynomial(Q, -sc);
            result += 1 / -sc;
        } else if (s <= 4) {
            // Maximum Deviation Found:                     5.946e-22
            // Expected Error Term:                        -5.946e-22
            // LCOV_EXCL_START
            double Y = 0.6986598968505859375;
            double[] P = {
                -0.0537258300023595030676,
                0.0445163473292365591906,
                0.0128677673534519952905,
                0.00097541770457391752726,
                0.769875101573654070925e-4,
                0.328032510000383084155e-5,
            };
            double[] Q = {
                1.0f,
                0.33383194553034051422,
                0.0487798431291407621462,
                0.00479039708573558490716,
                0.000270776703956336357707,
                0.106951867532057341359e-4,
                0.236276623974978646399e-7,
            };
            // LCOV_EXCL_STOP
            result = evaluatePolynomial(P, s - 2) / evaluatePolynomial(Q, s - 2);
            result += Y + 1 / -sc;
        } else if (s <= 7) {
            // Maximum Deviation Found:                     2.955e-17
            // Expected Error Term:                         2.955e-17
            // Max error found at double precision:         2.009135e-16
            // LCOV_EXCL_START
            double[] P = {
                -2.49710190602259410021,
                -2.60013301809475665334,
                -0.939260435377109939261,
                -0.138448617995741530935,
                -0.00701721240549802377623,
                -0.229257310594893932383e-4,
            };
            double[] Q = {
                1.0f,
                0.706039025937745133628,
                0.15739599649558626358,
                0.0106117950976845084417,
                -0.36910273311764618902e-4,
                0.493409563927590008943e-5,
                -0.234055487025287216506e-6,
                0.718833729365459760664e-8,
                -0.1129200113474947419e-9,
            };
            // LCOV_EXCL_STOP
            result = evaluatePolynomial(P, s - 4) / evaluatePolynomial(Q, s - 4);
            result = 1 + Math.exp(result);
        } else if (s < 15) {
            // Maximum Deviation Found:                     7.117e-16
            // Expected Error Term:                         7.117e-16
            // Max error found at double precision:         9.387771e-16
            // LCOV_EXCL_START
            double[] P = {
                -4.78558028495135619286,
                -1.89197364881972536382,
                -0.211407134874412820099,
                -0.000189204758260076688518,
                0.00115140923889178742086,
                0.639949204213164496988e-4,
                0.139348932445324888343e-5,
            };
            double[] Q = {
                1.0f,
                0.244345337378188557777,
                0.00873370754492288653669,
                -0.00117592765334434471562,
                -0.743743682899933180415e-4,
                -0.21750464515767984778e-5,
                0.471001264003076486547e-8,
                -0.833378440625385520576e-10,
                0.699841545204845636531e-12,
            };
            // LCOV_EXCL_STOP
            result = evaluatePolynomial(P, s - 7) / evaluatePolynomial(Q, s - 7);
            result = 1 + Math.exp(result);
        } else if (s < 36) {
            // Max error in interpolated form:              1.668e-17
            // Max error found at long double precision:    1.669714e-17
            // LCOV_EXCL_START
            double[] P = {
                -10.3948950573308896825,
                -2.85827219671106697179,
                -0.347728266539245787271,
                -0.0251156064655346341766,
                -0.00119459173416968685689,
                -0.382529323507967522614e-4,
                -0.785523633796723466968e-6,
                -0.821465709095465524192e-8,
            };
            double[] Q = {
                1.0f,
                0.208196333572671890965,
                0.0195687657317205033485,
                0.00111079638102485921877,
                0.408507746266039256231e-4,
                0.955561123065693483991e-6,
                0.118507153474022900583e-7,
                0.222609483627352615142e-14,
            };
            // LCOV_EXCL_STOP
            result = evaluatePolynomial(P, s - 15) / evaluatePolynomial(Q, s - 15);
            result = 1 + Math.exp(result);
        } else {
            result = 1 + Math.pow(2, -s);
        }
        return result;
    }

    /**
     * Evaluate the polynomial using Horner's method.
     * The coefficients are used in descending order, for example a polynomial of order
     * 3 requires 4 coefficients:
     * <pre>
     * f(x) = c[3] * x^3 + c[2] * x^2 + c[1] * x + c[0]
     * </pre>
     *
     * @param c Polynomial coefficients (must have {@code length > 0})
     * @param x Argument x
     * @return polynomial value
     */
    private static double evaluatePolynomial(double[] c, double x) {
        final int count = c.length;
        double sum = c[count - 1];
        for (int i = count - 2; i >= 0; --i) {
            sum *= x;
            sum += c[i];
        }
        return sum;
    }
}
