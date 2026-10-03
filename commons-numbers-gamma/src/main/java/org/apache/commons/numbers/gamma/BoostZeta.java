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
    /** ln(2pi). Computed to 30-digits precision. */
    private static final double LOG_2PI = 1.83787706640934548356065947281;
    /** Value for {@code s < 0} where zeta(s) is +/- infinity for all
     * {@code s != -2n}. Note the zeta function for negative arguments oscillates
     * above and below zero with roots at even {@code s}. The value is increasingly
     * large with larger negative {@code s} until all non-zero values cannot be represented
     * as a double. */
    private static final double LARGE_NEGATIVE_S = -267;
    /** Exact result for zeta(-2n + 1). The Boost implementation uses precomputed
     * Bernoulli numbers for B_2n divided by -2n. Here we tabulate the zeta
     * result as the table size is the same. All other odd negative s are infinite.
     * Computed using mpmath (1.14.1). Max n = 129. */
    private static final double[] ZN = {
        -0.0833333333333333333333333333333,
        0.00833333333333333333333333333333,
        -0.00396825396825396825396825396825,
        0.00416666666666666666666666666667,
        -0.00757575757575757575757575757576,
        0.0210927960927960927960927960928,
        -0.0833333333333333333333333333333,
        0.443259803921568627450980392157,
        -3.05395433027011974380395433027,
        26.4562121212121212121212121212,
        -281.460144927536231884057971014,
        3607.5105463980463980463980464,
        -54827.5833333333333333333333333,
        974936.82385057471264367816092,
        -20052695.7966880789461434622725,
        472384867.721629901960784313725,
        -12635724795.9166666666666666667,
        380879311252.453688115530220793,
        -12850850499305.0833333333333333,
        482414483548501.703715816703622,
        -20040310656516252.7381084216632,
        916774360319533077.569927536232,
        -45979888343656503490.4379432624,
        2518047192145109569708.90233202,
        -150017334921539287337114.401515,
        9689957887463594065649794.28946,
        -676458823792928209909452423.018,
        50890659468662289689766332915.9,
        -4.11472887925579786976654860676e+30,
        3.56665820953755561096845746087e+32,
        -3.30660898765775767256802146704e+34,
        3.27156342364787162642112270157e+36,
        -3.447378255827805387825645508e+38,
        3.86142798327052588930927202002e+40,
        -4.58929744324543321688639890061e+42,
        5.77753863427704318248848256879e+44,
        -7.69198587595071351674100759718e+46,
        1.08136354499716546963540333511e+49,
        -1.60293645220089654060671023458e+51,
        2.50194790415604628436566614985e+53,
        -4.10670523358102124797520450041e+55,
        7.07987744084945806174529724334e+57,
        -1.28045468879395087901908497563e+60,
        2.42673403923335240780208920671e+62,
        -4.8143218874045769355129570066e+64,
        9.98755741757275306806527774082e+66,
        -2.16456348684351856313351361598e+69,
        4.89623270396205532068492245156e+71,
        -1.15490239239635196639542716916e+74,
        2.83822495706937069592641563365e+76,
        -7.26120088036067163036772815107e+78,
        1.93235142334198120033323266084e+81,
        -5.34501604252886240053956094628e+83,
        1.53560288464224230702071420133e+86,
        -4.57898726822657976538994744683e+88,
        1.41620252121948092583601799759e+91,
        -4.54006522960926552491870532303e+93,
        1.50766567588078597755948498945e+96,
        -5.18309491482645637761224790375e+98,
        1.84356474272565291185736028806e+101,
        -6.78055547530909588969025102131e+103,
        2.577332670275460450289647933e+106,
        -1.0119112875704597605007955796e+109,
        4.10163461615422921089084567379e+111,
        -1.71552445340320193922071524606e+114,
        7.40034257052690942716920556053e+116,
        -3.29092253570544434867706140857e+119,
        1.50798315341647712056833365644e+122,
        -7.11698791882545486286760649291e+124,
        3.45804291415777717919922833643e+127,
        -1.72909076066767483167489207071e+130,
        8.89369916950329690887674533236e+132,
        -4.7038470619636014515138279862e+135,
        2.55719382310602058749858077138e+138,
        -1.42840675004435277005808820901e+141,
        8.1952152218313782940918703695e+143,
        -4.82764854227273717816101742819e+146,
        2.91896123747703236500405982201e+149,
        -1.81089321625689040160530678804e+152,
        1.15235772200211685798051266585e+155,
        -7.51923119519817697500081265839e+157,
        5.02940165764110497246840522742e+160,
        -3.4473420444477676704609427599e+163,
        2.42074586458685147183142674899e+166,
        -1.74094659203776765075736879892e+169,
        1.28194898634822427378088228066e+172,
        -9.66241211085609184243169684778e+174,
        7.45269103043008957309390945203e+177,
        -5.8808393311674371248220704455e+180,
        4.74627186549076153992212525722e+183,
        -3.91691325947728254682903333391e+186,
        3.30450714432260322283069086243e+189,
        -2.84928905509945827581152102174e+192,
        2.51033293450775865129599057986e+195,
        -2.25939019954752532049562261337e+198,
        2.07691380042876080434623778929e+201,
        -1.94947321749272591308734114001e+204,
        1.86807314712659138998396980689e+207,
        -1.8270752662814576943866394409e+210,
        1.82353863225956771810691544328e+213,
        -1.85686908101259450981907133715e+216,
        1.92871898511956020928868278693e+219,
        -2.04311704602864475750762704423e+222,
        2.20684116445278455076828336814e+225,
        -2.4300821796490274251390389767e+228,
        2.72748878790834695290272297756e+231,
        -3.11974215737550845945157847856e+234,
        3.63589387242826001493479749833e+237,
        -4.31683000307608832681396003788e+240,
        5.22042448793871999720448178213e+243,
        -6.42926069497693048519893333726e+246,
        8.06230338701308438135204143382e+249,
        -1.02927147379030111575795568666e+253,
        1.33753296997805240211378429501e+256,
        -1.76894809027973797575657445272e+259,
        2.38064790180923972522551743144e+262,
        -3.2597127947194184823502123513e+265,
        4.54049623716012131918616627123e+268,
        -6.43285751931478506106894605686e+271,
        9.26870486757493111152509786938e+274,
        -1.35796195002851814738921378693e+278,
        2.02278397360493216811767187783e+281,
        -3.06299069922083360655392703162e+284,
        4.71430853007426521282616631982e+287,
        -7.37410458713557576506584806391e+290,
        1.17209627670508265765085284663e+294,
        -1.89288666446856573885385316552e+297,
        3.10555175960489268960251418622e+300,
        -5.1754977470366797965163888379e+303,
        8.76015634462292151490407301349e+306,
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
        if (Double.isNaN(s) || s == Double.NEGATIVE_INFINITY) {
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
        // Changed threshold from 53 to a higher value where the extended precision result is 1.0.
        //
        if (s > 53.0000000001) {
            return 1;
        }
        //
        // Start by seeing if we have a simple closed form for integer s:
        //
        if (Math.floor(s) == s && s < 0) {
            // Change from Boost implementation to only handle negative s.
            // Handle any negative even integer. Finite values below -(2^63) are clipped
            // to -(2^63) and will be even.
            if (((long) s & 1) == 0) {
                // Negative even integer
                return 0;
            }

            // Special handling for small integer s using Bernoulli number B_{2n} (B2N).
            // Note: Change from the Boost implementation:
            //
            // Remove positive even integer case:
            //   2^(v-1) * pow(pi, v) * abs(B2N[v/2]) / v!
            // Code to handle positive even integers is less accurate than the rational
            // function approximation. zetaImp53 is exact except for 1 ULP at s=2.
            //
            // Remove positive odd integer case:
            // For odd integers zetaImp53 is exact except for 1 ULP at s=53. This
            // is exact if the asymptote also includes pow(3, -s).
            final int v = (int) s;
            if (v == s) {
                // Negative odd small integer.
                // Result:
                // n = (1 - v) / 2
                // -B2N[n] / (1 - v)
                final int n = -v / 2;
                if (n < ZN.length) {
                    // Tabulated result to avoid B2N division
                    return ZN[n];
                }
            }
        }

        double result;
        if (Math.abs(s) < BoostGamma.ROOT_EPSILON) {
            result = -0.5 - BoostGamma.LOG_ROOT_TWO_PI * s;
        } else if (s < 0) {
            // Negative; |s| > small; and not an even integer (all even integers handled above).
            // This ensures s/2 is not odd and avoids sin(pi * s/2) = 0.
            // This would generate 0 * infinity = NaN for large |s|.

            // Change from Boost implementation.
            // zeta(s) where the value either side of even s overflows
            if (s <= LARGE_NEGATIVE_S) {
                // Set the sign without using sinp.
                // Bypasses infinity * sinp(0.5 * s) when the reflection formula will overflow.
                // Assumes 0.5 * s is non-integer since zeta(-2n) = 0.0.
                return SpecialMath.isOdd(Math.floor(0.5 * s)) ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            }

            // Swap: s is now positive; sc = 1 - s
            final double tmp = s;
            s = sc;
            sc = tmp;

            if (s > BoostGamma.MAX_FACTORIAL) {
                // This has been simplified from the Boost implementation which will
                // catch overflow conditions when compiled with an appropriate evaluation
                // policy and return signed infinity, or raise an error. Java floating-point
                // arithmetic does not create overflow exceptions and will return infinity.
                // Note that zeta(s > 170) == 1.0
                final double mult = BoostGamma.sinp(0.5 * sc) * 2; // * zeta(s)
                result = LogGamma.value(s);
                result -= s * LOG_2PI;
                // Possible overflow if result > 709
                // Use exp(2x) = exp(x) * exp(x) knowing |mult| in (0, 2]
                if (result > BoostGamma.LOG_MAX_VALUE) {
                    result = Math.exp(result * 0.5);
                    result = result * mult * result;
                } else {
                    result = Math.exp(result) * mult;
                }
            } else {
                result = BoostGamma.sinp(0.5 * sc) *
                    2 * Math.pow(2 * Math.PI, -s) *
                    Gamma.value(s) *
                    zetaImp(s, sc);
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
     * @param s Argument (assumed to be positive finite).
     * @param sc Complement argument (1 - s).
     * @return zeta(s)
     */
    static double zetaImp53(double s, double sc) {
        double result;
        if (s < 1) {
            // [-0.5, -9007199254740991.42278433509847]
            // Rational Approximation
            // Maximum Deviation Found:                     2.020e-18
            // Expected Error Term:                        -2.020e-18
            // Max error found at double precision:         3.994987e-17
            double P;
            P = -0.933241270357061460782e-5;
            P =  0.000451534528645796438704 + P * sc;
            P =  -0.00320912498879085894856 + P * sc;
            P =    0.0557616214776046784287 + P * sc;
            P =     -0.49092470516353571651 + P * sc;
            P =      0.24339294433593750202 + P * sc;
            double Q;
            Q = -0.101855788418564031874e-4;
            Q =   0.00024978985622317935355 + Q * sc;
            Q =  -0.00413421406552171059003 + Q * sc;
            Q =    0.0419676223309986037706 + Q * sc;
            Q =    -0.279960334310344432495 + Q * sc;
            Q =                           1 + Q * sc;
            result = P / Q;
            result -= 1.2433929443359375F;
            result += sc;
            result /= sc;
        } else if (s <= 2) {
            // [4503599627370496.57721566490153, 1.64493406684822643647241516665]
            // Maximum Deviation Found:                     9.007e-20
            // Expected Error Term:                         9.007e-20
            double P;
            P = 0.110108440976732897969e-4;
            P = 0.000249606367151877175456 + P * -sc;
            P =  0.00390252087072843288378 + P * -sc;
            P =   0.0417364673988216497593 + P * -sc;
            P =    0.243210646940107164097 + P * -sc;
            P =    0.577215664901532860516 + P * -sc;
            double Q;
            Q =  0.10991819782396112081e-4;
            Q = 0.000255784226140488490982 + Q * -sc;
            Q =  0.00434930582085826330659 + Q * -sc;
            Q =    0.043460910607305495864 + Q * -sc;
            Q =    0.295201277126631761737 + Q * -sc;
            Q =                        1.0 + Q * -sc;
            result = P / Q;
            result += 1 / -sc;
        } else if (s <= 4) {
            // [1.64493406684822602011735171122, 1.08232323371113819151600369654]
            // Maximum Deviation Found:                     5.946e-22
            // Expected Error Term:                        -5.946e-22
            final double Y = 0.6986598968505859375;
            final double x = s - 2;
            double P;
            P = 0.328032510000383084155e-5;
            P = 0.769875101573654070925e-4 + P * x;
            P =  0.00097541770457391752726 + P * x;
            P =   0.0128677673534519952905 + P * x;
            P =   0.0445163473292365591906 + P * x;
            P =  -0.0537258300023595030676 + P * x;
            double Q;
            Q = 0.236276623974978646399e-7;
            Q = 0.106951867532057341359e-4 + Q * x;
            Q = 0.000270776703956336357707 + Q * x;
            Q =  0.00479039708573558490716 + Q * x;
            Q =   0.0487798431291407621462 + Q * x;
            Q =     0.33383194553034051422 + Q * x;
            Q =                        1.0 + Q * x;
            result = P / Q;
            result += Y + 1 / -sc;
        } else if (s <= 7) {
            // [1.08232323371113813031050445339, 1.00834927738192282683979754985]
            // Maximum Deviation Found:                     2.955e-17
            // Expected Error Term:                         2.955e-17
            // Max error found at double precision:         2.009135e-16
            final double x = s - 4;
            double P;
            P = -0.229257310594893932383e-4;
            P =  -0.00701721240549802377623 + P * x;
            P =    -0.138448617995741530935 + P * x;
            P =    -0.939260435377109939261 + P * x;
            P =     -2.60013301809475665334 + P * x;
            P =     -2.49710190602259410021 + P * x;
            double Q;
            Q =   -0.1129200113474947419e-9;
            Q =  0.718833729365459760664e-8 + Q * x;
            Q = -0.234055487025287216506e-6 + Q * x;
            Q =  0.493409563927590008943e-5 + Q * x;
            Q =  -0.36910273311764618902e-4 + Q * x;
            Q =    0.0106117950976845084417 + Q * x;
            Q =      0.15739599649558626358 + Q * x;
            Q =     0.706039025937745133628 + Q * x;
            Q =                         1.0 + Q * x;
            result = P / Q;
            result = 1 + Math.exp(result);
        } else if (s < 15) {
            // [1.00834927738192282148095799031, 1.00003058823630702053126571251]
            // Maximum Deviation Found:                     7.117e-16
            // Expected Error Term:                         7.117e-16
            // Max error found at double precision:         9.387771e-16
            final double x = s - 7;
            double P;
            P =  0.139348932445324888343e-5;
            P =  0.639949204213164496988e-4 + P * x;
            P =   0.00115140923889178742086 + P * x;
            P = -0.000189204758260076688518 + P * x;
            P =    -0.211407134874412820099 + P * x;
            P =     -1.89197364881972536382 + P * x;
            P =     -4.78558028495135619286 + P * x;
            double Q;
            Q =  0.699841545204845636531e-12;
            Q = -0.833378440625385520576e-10 + Q * x;
            Q =   0.471001264003076486547e-8 + Q * x;
            Q =   -0.21750464515767984778e-5 + Q * x;
            Q =  -0.743743682899933180415e-4 + Q * x;
            Q =   -0.00117592765334434471562 + Q * x;
            Q =    0.00873370754492288653669 + Q * x;
            Q =      0.244345337378188557777 + Q * x;
            Q =                          1.0 + Q * x;
            result = P / Q;
            result = 1 + Math.exp(result);
        } else if (s < 36) {
            // [1.00003058823630702049355172851, 1.00000000001455192189104205591]
            // Max error in interpolated form:              1.668e-17
            // Max error found at long double precision:    1.669714e-17
            final double x = s - 15;
            double P;
            P = -0.821465709095465524192e-8;
            P = -0.785523633796723466968e-6 + P * x;
            P = -0.382529323507967522614e-4 + P * x;
            P =  -0.00119459173416968685689 + P * x;
            P =   -0.0251156064655346341766 + P * x;
            P =    -0.347728266539245787271 + P * x;
            P =     -2.85827219671106697179 + P * x;
            P =     -10.3948950573308896825 + P * x;
            double Q;
            Q = 0.222609483627352615142e-14;
            Q =  0.118507153474022900583e-7 + Q * x;
            Q =  0.955561123065693483991e-6 + Q * x;
            Q =  0.408507746266039256231e-4 + Q * x;
            Q =   0.00111079638102485921877 + Q * x;
            Q =    0.0195687657317205033485 + Q * x;
            Q =     0.208196333572671890965 + Q * x;
            Q =                         1.0 + Q * x;
            result = P / Q;
            result = 1 + Math.exp(result);
        } else {
            // [1.00000000001455192189104198424, 1.0]
            // Change from: 1 + Math.pow(2, -s);
            // Adding 3^-s increases ULP accuracy as the result approaches 1.0
            result = Math.pow(3, -s) + Math.pow(2, -s) + 1;
        }
        return result;
    }
}
