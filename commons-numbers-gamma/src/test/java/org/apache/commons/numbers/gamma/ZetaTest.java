/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.apache.commons.numbers.gamma;

import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.SplittableRandom;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.apache.commons.numbers.core.DD;
import org.apache.commons.numbers.fraction.BigFraction;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.MethodOrderer;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestMethodOrder;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * Test the {@link HurwitzZeta} function when {@code a=1}, i.e. the Riemann zeta function.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class ZetaTest {
    /** Flag set when the JVM version is printed. Used for testing. */
    private static boolean jvm = false;

    /**
     * Pre-computed table of Borwein coefficients (-1)^k * (d_n - d_k) / d_n for n = 24.
     * Borwein's algorithm requires approximately n = 1.3d for d digits of precision.
     *
     * <p>Borwein coefficient:
     * d_k = n * sum_{i=0}^k ( (n + i - 1)! * 4^i / ( (n - i)! * (2i)! ) )
     */
    private static final double[] DN_DK = {
        1.0,
        -0.999999999999999,
        0.999999999999812,
        -0.9999999999855517,
        0.9999999994080065,
        -0.9999999850335489,
        0.9999997450236661,
        -0.9999968965547265,
        0.999971877502541,
        -0.9998044297284367,
        0.998931938694946,
        -0.995336218072075,
        0.9834807623952175,
        -0.9519634894573569,
        0.8840929599033391,
        -0.7655145634411474,
        0.5976878813515132,
        -0.4062278518731426,
        0.23178649168173823,
        -0.1067246914876162,
        0.03778036573957456,
        -0.00959406764049133,
        0.0015493525382159912,
        -1.1918096447815317E-4,
    };

    /** ln(2). */
    private static final double LN2 = Math.log(2);


    /** Define the expected error for a test. */
    private interface TestError {
        /**
         * @return maximum allowed error
         */
        double getTolerance();

        /**
         * @return maximum allowed RMS error
         */
        double getRmsTolerance();
    }

    /**
     * Define the test cases for each resource file.
     * This encapsulates the function to test, the expected maximum and RMS error, and
     * the resource file containing the data.
     *
     * <p>The Boost functions use the default policy of internal promotion
     * of double to long double if it offers more precision. Code comments
     * in the implementations for the maximum error are using the defaults with
     * promotion enabled where the error is 'effectively zero'. Java does not
     * support long double computation. Tolerances have been set to allow tests to
     * pass. Spot checks on larger errors have been verified against the reference
     * implementation compiled with promotion of double <strong>disabled</strong>.
     *
     * @see <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/relative_error.html">Relative error</a>
     * @see <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/pol_tutorial/policy_tut_defaults.html">Policy defaults</a>
     */
    private enum TestCase implements TestError {
        // s in (1, 32); a in [1, 2^31)
        BORWEIN_ZETA_INT(ZetaTest::borweinZeta, "zeta.csv", 2.8, 0.66),
        HZETA_INT(s -> HurwitzZeta.value(s, 1), "zeta.csv", 1.68, 0.5),
        ZETA_INT(ZetaTest::zeta_imp_prec, "zeta.csv", 1.5, 0.34)
        ;

        /** The function. */
        private final DoubleUnaryOperator fun;

        /** The filename containing the test data. */
        private final String filename;

        /** The field containing the expected value. */
        private final int expected;

        /** The maximum allowed ulp. */
        private final double maxUlp;

        /** The maximum allowed RMS ulp. */
        private final double rmsUlp;

        /**
         * Instantiates a new test case.
         *
         * @param fun function to test
         * @param filename Filename of the test data
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         */
        TestCase(DoubleUnaryOperator fun, String filename, double maxUlp, double rmsUlp) {
            this(fun, filename, 1, maxUlp, rmsUlp);
        }

        /**
         * Instantiates a new test case.
         *
         * @param fun function to test
         * @param filename Filename of the test data
         * @param expected Expected result field index
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         */
        TestCase(DoubleUnaryOperator fun, String filename, int expected, double maxUlp, double rmsUlp) {
            this.fun = fun;
            this.filename = filename;
            this.expected = expected;
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
        }

        /**
         * @return function to test
         */
        DoubleUnaryOperator getFunction() {
            return fun;
        }

        /**
         * @return Filename of the test data
         */
        String getFilename() {
            return filename;
        }

        /**
         * @return Expected result field index
         */
        int getExpectedField() {
            return expected;
        }

        @Override
        public double getTolerance() {
            return maxUlp;
        }

        @Override
        public double getRmsTolerance() {
            return rmsUlp;
        }
    }

    // Test zeta implementation

    /**
     * Compute the value of the Riemann zeta function {@code zeta(s)}. Uses the Borwein summation.
     * Equations provided in Belovas, eq 1.1 - 1.2.
     *
     * <ol>
     * <li>Borwein, P (1995)
     * An efficient algorithm for the Riemann zeta function.
     * Constructive experimental and nonlinear analysis, CMS Conference Proceedings 27, pp. 29–34.
     * <li>Belovas, I (2019)
     * A central limit theorem for coefficients of the modified Borwein method for the
     * calculation of the Riemann zeta-function.
     * Lith Math J, pp. 17–23.
     * <a href="https://doi.org/10.1007/s10986-019-09421-4">10.1007/s10986-019-09421-4</a>
     * </ol>
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     *
     * @param s Argument {@code s > 1}
     * @return zeta(s)
     */
    static double borweinZeta(double s) {
        // Asymptote of Riemann zeta function: 1 + 2^-s [+ 3^-s + ...]
        // When 3^-s and later terms cannot be added.
        if (s > 36) {
            return 1 + Math.pow(2, -s);
        }

        // Evaluate the accelerated Dirichlet eta series alternation
        // sum_{k=0}^{n-1} (-1)^k * (d_n - d_k) / (k + 1)^s

        DD sum = DD.ZERO;
        // sum terms by magnitude
        for (int k = DN_DK.length - 1; k > 0; k--) {
            sum = sum.add(
                DD.ofProduct(DN_DK[k], Math.pow(k + 1, -s)));
        }
        // Final term
        sum = sum.add(DD.ONE);

        // Normalise by 1 / (d_n * (1 - 2^(1-s)))
        // The d_n factor has been incorporated into the coefficients
        return sum.divide(-Math.expm1(LN2 * (1 - s))).hi();
    }

    static double zeta_polynomial_series(double s) // , double sc)
    {
        //
        // This is algorithm 3 from:
        //
        // "An Efficient Algorithm for the Riemann Zeta Function", P. Borwein,
        // Canadian Mathematical Society, Conference Proceedings.
        // See: http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
        //
        int n = 18; // itrunc(T(log(boost::math::tools::epsilon<T>()) / -2));
        double sum = 0;
        double two_n = Math.scalb(1, n); // ldexp(T(1), n);
        int ej_sign = 1;
        for (int j = 0; j < n; ++j) {
            sum += ej_sign * -two_n / Math.pow(j + 1, s);
            ej_sign = -ej_sign;
        }
        double ej_sum = 1;
        double ej_term = 1;
        for (int j = n; j <= 2 * n - 1; ++j) {
            // XXX: These could be precomputed.
            // Do not know what they are.
            sum += ej_sign * (ej_sum - two_n) / Math.pow(j + 1, s);
            ej_sign = -ej_sign;
            ej_term *= 2 * n - j;
            ej_term /= j - n + 1;
            ej_sum += ej_term;
        }

//        // No difference
//        BigFraction ej_sum = BigFraction.ONE;
//        BigFraction ej_term = BigFraction.ONE;
//        BigFraction two_nf = BigFraction.from(two_n);
//        for (int j = n; j <= 2 * n - 1; ++j) {
//            sum += ej_sign * ej_sum.subtract(two_nf).doubleValue() / Math.pow(j + 1, s);
//            ej_sign = -ej_sign;
//            ej_term = ej_term.multiply(BigFraction.of(2 * n - j, j - n + 1));
//            ej_sum = ej_sum.add(ej_term);
//        }

        // return -sum / (two_n * (-powm1(T(2), sc)));
        return -sum / (two_n * (-Math.expm1(LN2 * (1 - s))));
    }

    static double zeta_imp_prec(double s) {
        return zeta_imp_prec(s, 1-s);
    }

    static double zeta_imp_prec(double s, double sc)
    {
       double result;
       if(s < 1)
       {
          // Rational Approximation
          // Maximum Deviation Found:                     2.020e-18
          // Expected Error Term:                         -2.020e-18
          // Max error found at double precision:         3.994987e-17
          // LCOV_EXCL_START
          double[] P = {
             (0.24339294433593750202),
             (-0.49092470516353571651),
             (0.0557616214776046784287),
             (-0.00320912498879085894856),
             (0.000451534528645796438704),
             (-0.933241270357061460782e-5),
            };
          double[] Q = {
             (1),
             (-0.279960334310344432495),
             (0.0419676223309986037706),
             (-0.00413421406552171059003),
             (0.00024978985622317935355),
             (-0.101855788418564031874e-4),
          };
          // LCOV_EXCL_STOP
          result = evaluatePolynomial(P, sc) / evaluatePolynomial(Q, sc);
          result -= 1.2433929443359375F;
          result += (sc);
          result /= (sc);
       }
       else if(s <= 2)
       {
          // Maximum Deviation Found:        9.007e-20
          // Expected Error Term:            9.007e-20
          // LCOV_EXCL_START
          double[] P = {
             (0.577215664901532860516),
             (0.243210646940107164097),
             (0.0417364673988216497593),
             (0.00390252087072843288378),
             (0.000249606367151877175456),
             (0.110108440976732897969e-4),
          };
          double[] Q = {
             (1.0),
             (0.295201277126631761737),
             (0.043460910607305495864),
             (0.00434930582085826330659),
             (0.000255784226140488490982),
             (0.10991819782396112081e-4),
          };
          // LCOV_EXCL_STOP
          result = evaluatePolynomial(P, (-sc)) / evaluatePolynomial(Q, (-sc));
          result += 1 / (-sc);
       }
       else if(s <= 4)
       {
          // Maximum Deviation Found:          5.946e-22
          // Expected Error Term:              -5.946e-22
          // LCOV_EXCL_START
          double Y = 0.6986598968505859375;
          double[] P = {
             (-0.0537258300023595030676),
             (0.0445163473292365591906),
             (0.0128677673534519952905),
             (0.00097541770457391752726),
             (0.769875101573654070925e-4),
             (0.328032510000383084155e-5),
          };
          double[] Q = {
             1.0f,
             (0.33383194553034051422),
             (0.0487798431291407621462),
             (0.00479039708573558490716),
             (0.000270776703956336357707),
             (0.106951867532057341359e-4),
             (0.236276623974978646399e-7),
          };
          // LCOV_EXCL_STOP
          result = evaluatePolynomial(P, (s - 2)) / evaluatePolynomial(Q, (s - 2));
          result += Y + 1 / (-sc);
       }
       else if(s <= 7)
       {
          // Maximum Deviation Found:                     2.955e-17
          // Expected Error Term:                         2.955e-17
          // Max error found at double precision:         2.009135e-16
          // LCOV_EXCL_START
          double[] P = {
             (-2.49710190602259410021),
             (-2.60013301809475665334),
             (-0.939260435377109939261),
             (-0.138448617995741530935),
             (-0.00701721240549802377623),
             (-0.229257310594893932383e-4),
          };
          double[] Q = {
             1.0f,
             (0.706039025937745133628),
             (0.15739599649558626358),
             (0.0106117950976845084417),
             (-0.36910273311764618902e-4),
             (0.493409563927590008943e-5),
             (-0.234055487025287216506e-6),
             (0.718833729365459760664e-8),
             (-0.1129200113474947419e-9),
          };
          // LCOV_EXCL_STOP
          result = evaluatePolynomial(P, (s - 4)) / evaluatePolynomial(Q, (s - 4));
          result = 1 + Math.exp(result);
       }
       else if(s < 15)
       {
          // Maximum Deviation Found:                     7.117e-16
          // Expected Error Term:                         7.117e-16
          // Max error found at double precision:         9.387771e-16
          // LCOV_EXCL_START
          double[] P = {
             (-4.78558028495135619286),
             (-1.89197364881972536382),
             (-0.211407134874412820099),
             (-0.000189204758260076688518),
             (0.00115140923889178742086),
             (0.639949204213164496988e-4),
             (0.139348932445324888343e-5),
            };
          double[] Q = {
             1.0f,
             (0.244345337378188557777),
             (0.00873370754492288653669),
             (-0.00117592765334434471562),
             (-0.743743682899933180415e-4),
             (-0.21750464515767984778e-5),
             (0.471001264003076486547e-8),
             (-0.833378440625385520576e-10),
             (0.699841545204845636531e-12),
            };
          // LCOV_EXCL_STOP
          result = evaluatePolynomial(P, (s - 7)) / evaluatePolynomial(Q, (s - 7));
          result = 1 + Math.exp(result);
       }
       else if(s < 36)
       {
          // Max error in interpolated form:             1.668e-17
          // Max error found at long double precision:   1.669714e-17
          // LCOV_EXCL_START
          double[] P = {
             (-10.3948950573308896825),
             (-2.85827219671106697179),
             (-0.347728266539245787271),
             (-0.0251156064655346341766),
             (-0.00119459173416968685689),
             (-0.382529323507967522614e-4),
             (-0.785523633796723466968e-6),
             (-0.821465709095465524192e-8),
          };
          double[] Q = {
             1.0f,
             (0.208196333572671890965),
             (0.0195687657317205033485),
             (0.00111079638102485921877),
             (0.408507746266039256231e-4),
             (0.955561123065693483991e-6),
             (0.118507153474022900583e-7),
             (0.222609483627352615142e-14),
          };
          // LCOV_EXCL_STOP
          result = evaluatePolynomial(P, (s - 15)) / evaluatePolynomial(Q, (s - 15));
          result = 1 + Math.exp(result);
       }
       else
       {
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
    static double evaluatePolynomial(double[] c, double x) {
        final int count = c.length;
        double sum = c[count - 1];
        for (int i = count - 2; i >= 0; --i) {
            sum *= x;
            sum += c[i];
        }
        return sum;
    }

    @Test
    void testK() {
        // Require n ~ 1.3d terms for d digits of precision
        final int n = 24;
        // Factorial up to 2n
        final BigInteger[] factorial = new BigInteger[2 * n + 1];
        BigInteger f = BigInteger.ONE;
        factorial[0] = f;
        factorial[1] = f;
        for (int i = 2; i < factorial.length; i++) {
            f = f.multiply(BigInteger.valueOf(i));
            factorial[i] = f;
        }
        // Compute d_k
        final BigFraction[] dk = new BigFraction[n + 1];
        final BigInteger four = BigInteger.valueOf(4);
        // d_k = n * sum_{i=0}^k ( (n + i - 1)! * 4^i / ( (n - i)! * (2i)! ) )
        for (int k = 0; k <= n; k++) {
            BigFraction sum = BigFraction.ZERO;
            for (int i = 0; i <= k; i++) {
                final BigInteger num = factorial[n + i - 1].multiply(four.pow(i));
                final BigInteger den = factorial[n - i].multiply(factorial[2 * i]);
                final BigFraction term = BigFraction.of(num, den);
                sum = sum.add(term);
            }
            dk[k] = sum.multiply(n);
            // The n factor cancels
            //dk[k] = sum;
        }
        // Create the table
        final BigFraction dn = dk[n];
        for (int k = 0; k < n; k++) {
            final double sign = Math.pow(-1, k);
            final double dn_dk = sign * dn.subtract(dk[k]).divide(dn).doubleValue();
            //System.out.printf("%s,%n", dn_dk);
            Assertions.assertEquals(DN_DK[k], dn_dk);
        }
        //System.out.printf("DN %s%n", dn.doubleValue());
    }

    /**
     * Spot tests for the zeta function to check various points in the domain and extreme values.
     */
    @ParameterizedTest
    @MethodSource
    void testZetaSpot(double s, double z, int ulp) {
        assertClose(ZetaTest::zeta_imp_prec, s, z, ulp);
    }

    static Stream<Arguments> testZetaSpot() {
        return Stream.of(
            // Reference values from mpmath version 1.4.1.
            // from mpmath import mp, zeta
            // mp.dps = 30; mp.pretty = True
            // def f(s):
            //   print(f'Arguments.of({s}, {zeta(s)}, 0),')
            // f(1.001) etc.

            // Borwein is suitable with Re(s) >= 0.5
            Arguments.of(0.75, -3.44128538694522289439513996071, 1),

            Arguments.of(1.0000000000000002, 4503599627370496.57721566490153, 0),
            Arguments.of(1.00000000001, 99999991726.5408003527075723256, 1),
            Arguments.of(1.000000001, 999999917.836851511852401825175, 1),
            Arguments.of(1.0000001, 10000000.5713770004182384783263, 1),
            Arguments.of(1.001, 1000.57728847601162684806668989, 0),
            Arguments.of(1.1678, 6.54877176355186372373480299218, 0),
            Arguments.of(1.3567, 3.40603490577277492861753138946, 0),
            Arguments.of(1.123456, 8.68618258920102287221702184772, 0),
            Arguments.of(1.268894, 4.31537649881889766400411555278, 0),
            Arguments.of(1.347949, 3.475936638461954139465311868, 1),
            Arguments.of(1.56235747, 2.39480872059997655621537802997, 0),
            Arguments.of(1.623897243, 2.22351825201319065543403824013, 1),
            Arguments.of(1.7237423, 2.00898039161124865158508485803, 0),
            Arguments.of(1.83254674, 1.83545995242142981354351854655, 0),
            Arguments.of(1.9236825, 1.72275978825438842136896455342, 0),
            Arguments.of(2.123794, 1.54242584690869940905284285109, 0),
            Arguments.of(2.32579, 1.41897842729817658361567621351, 1),
            Arguments.of(2.467548, 1.35437062072512909687849604628, 0),
            Arguments.of(2.66723684, 1.28401689949720753221761208015, 1),
            Arguments.of(2.728939234, 1.26600562506285409433958148773, 1),
            Arguments.of(2.926374752, 1.2173196265702590082342636002, 0),
            Arguments.of(3.0253674, 1.19710708981398055292662080543, 0),
            Arguments.of(3.1263478146, 1.17881947588107530387481077309, 0),
            Arguments.of(3.2347624, 1.16142938699197343431801274069, 1),
            Arguments.of(3.78979232, 1.09836739792891183083027283164, 0),
            Arguments.of(4.1231432789, 1.0743088573933574620870357475, 0),
            Arguments.of(4.26783234, 1.06598689761243309972313126874, 0),
            Arguments.of(4.4524634634, 1.05683164467801864447672646289, 0),
            Arguments.of(4.78787234, 1.04355994239549749069888177376, 0),
            Arguments.of(5.056958305089068, 1.03533807109327536672677745261, 0),
            Arguments.of(5.408082506656524, 1.02701946655787519577208051479, 0),
            Arguments.of(5.804759627466787, 1.0200524353843260063985878335, 0),
            Arguments.of(7.659570635240448, 1.00519791192315097478738463442, 0),
            Arguments.of(7.69218964845904, 1.00507813829841149853117735458, 0),
            Arguments.of(7.709385977606417, 1.00501613136233265276328382046, 0),
            Arguments.of(8.260524059512388, 1.0033882099196978922536586679, 0),
            Arguments.of(8.53828107111228, 1.00278281265556351136140118904, 0),
            Arguments.of(14.63270073489883, 1.00003947166440522694300950522, 0),
            Arguments.of(14.679887595148118, 1.00003819956745528232273286406, 0),
            Arguments.of(15.141806830586567, 1.00002772105376629598554243201, 0),
            Arguments.of(18.228630385704083, 1.00000325765188598303313821034, 0),
            Arguments.of(19.500685366621884, 1.00000134855680639452004367183, 0),
            Arguments.of(20.47205107115528, 1.00000068771214939325961623929, 0),
            Arguments.of(21.79946922232697, 1.00000027401160294459031905157, 0),
            Arguments.of(21.97996251516234, 1.00000024178569421473939141652, 0),
            Arguments.of(22.582265507388463, 1.00000015925996627732636438428, 0),
            Arguments.of(23.42093612093932, 1.00000008904886011938510467284, 0),
            Arguments.of(24.50652339609889, 1.00000004195873582436774985101, 0),
            Arguments.of(24.955890659420607, 1.00000003072881872070260829246, 0),

            // Asymptote of Riemann zeta function: 1 + 2^-s [ + 3^-s + ... ]
            // Threshold 3^-s < 2^-52 : s = 52 * log(2) / log(3) = 33.43
            Arguments.of(25.67, 1.00000001873152459253171691257, 0),
            Arguments.of(29.67, 1.00000000117069191296622655835, 0),
            Arguments.of(34.67, 1.0000000000365839328567015569213203977130966754644, 0),
            Arguments.of(35.67, 1.0000000000182919616411105681150690755793311651307, 0),
            Arguments.of(36.67, 1.0000000000091459792248364536695784035262447502029, 0),
            Arguments.of(39.67, 1.00000000000114324712238181431, 0),
            Arguments.of(47.2342, 1.00000000000000604072382727681, 0),
            Arguments.of(50.2342, 1.00000000000000075509047585184, 0),
            Arguments.of(51.2342, 1.00000000000000037754523774643, 0),
            Arguments.of(52.2342, 1.00000000000000018877261881338, 0),
            Arguments.of(53.2342, 1.00000000000000009438630938675, 0),
            Arguments.of(59.67, 1.00000000000000000109028530523, 0),
            Arguments.of(69.67, 1.00000000000000000000106473174, 0),
            Arguments.of(89.67, 1.00000000000000000000000000102, 0),

            // s -> large, a == 1
            Arguments.of(1001, 1.0, 0),
            Arguments.of(1e+18, 1.0, 0),
            Arguments.of(1e+19, 1.0, 0)
        );
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class)
    void testZeta(TestCase tc) {
        assertFunction(tc);
    }

    /**
     * Assert the function is close to the expected value.
     *
     * @param fun Function
     * @param x Input value
     * @param expected Expected value
     * @param tolerance the tolerance
     */
    private static void assertClose(DoubleUnaryOperator fun, double x, double expected, int tolerance) {
        final double actual = fun.applyAsDouble(x);
        TestUtils.assertEquals(expected, actual, tolerance, null, () -> Double.toString(x));
    }

    /**
     * Assert the function using extended precision.
     *
     * @param tc Test case
     */
    private static void assertFunction(TestCase tc) {
        final TestUtils.ErrorStatistics stats = new TestUtils.ErrorStatistics();
        try (DataReader in = new DataReader(tc.getFilename())) {
            while (in.next()) {
                try {
                    final double x = in.getDouble(0);
                    final BigDecimal expected = in.getBigDecimal(tc.getExpectedField());
                    final double actual = tc.getFunction().applyAsDouble(x);
                    TestUtils.assertEquals(expected, actual, tc.getTolerance(), stats::add,
                        () -> tc + " x=" + x);
                } catch (final NumberFormatException ex) {
                    Assertions.fail("Failed to load data: " + Arrays.toString(in.getFields()), ex);
                }
            }
        } catch (final IOException ex) {
            Assertions.fail("Failed to load data: " + tc.getFilename(), ex);
        }

        assertRms(tc, stats);
    }

    /**
     * Assert the Root Mean Square (RMS) error of the function is below the allowed
     * maximum for the specified TestError.
     *
     * @param te Test error
     * @param stats Error statistics
     */
    private static void assertRms(TestError te, TestUtils.ErrorStatistics stats) {
        final double rms = stats.getRMS();
        debugRms(te.toString(), stats.getMaxAbs(), rms, stats.getMean(), stats.size());
        Assertions.assertTrue(rms <= te.getRmsTolerance(),
            () -> String.format("%s RMS %s < %s", te, rms, te.getRmsTolerance()));
    }

    /**
     * Output the maximum and RMS ulp for the named test. Used for reporting the
     * errors and setting appropriate test tolerances. This is relevant across
     * different JDK implementations where the java.util.Math functions used in
     * BoostGamma may compute to different accuracy.
     *
     * @param name Test name
     * @param maxAbsUlp Maximum |ulp|
     * @param rmsUlp RMS ulp
     * @param meanUlp Mean ulp
     * @param size Number of measurements
     */
    private static void debugRms(String name, double maxAbsUlp, double rmsUlp, double meanUlp, int size) {
        // CHECKSTYLE: stop regexp
        if (!jvm) {
            jvm = true;
                System.out.printf("// %s %s%n",
                System.getProperty("java.vm.vendor"),
                System.getProperty("java.vm.version")
            );
        }
        System.out.printf("%-35s   max %10.6g   RMS %10.6g   mean %14.6g  n %4d%n",
            name, maxAbsUlp, rmsUlp, meanUlp, size);
    }

//    @Test
    void testSample() throws IOException {
        SplittableRandom rng = new SplittableRandom();
        // Create node points using 1 + 2^-b
        // Final node point is large s.
        double[] nodes = new double[58];
        double b = 0x1p-52;
        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = 1 + b;
            b *= 2;
        }

        int[] size = new int[nodes.length - 1];
        Arrays.fill(size, 10000 / size.length);

        try (PrintStream out = new PrintStream(Files.newOutputStream(Paths.get("zeta.txt")))) {
            for (int j = 0; j < nodes.length - 1; j++) {
                sample(rng, size[j], nodes[j], nodes[j + 1], out);
            }
        }
    }

    /**
     * Sample uniformly from the double values within the given range. This has a a log-uniform
     * distribution as the limiting distribution. If the range is smaller than the maximum
     * sample size then the sample is all double values in the enumerated range.
     *
     * @param rng the source of randomness
     * @param n the maximum number of samples
     * @param lo the low bound (inclusive)
     * @param hi the hi bound (exclusive)
     */
    private static void sample(SplittableRandom rng, int n, double lo, double hi,
            PrintStream out) {
        assert lo >= 0;
        assert hi > lo;
        long lb = Double.doubleToRawLongBits(lo);
        long hb = Double.doubleToRawLongBits(hi);
        long range = hb - lb;
        if (range < n) {
            out.printf("# [%s, %s) %d / %d%n", lo, hi, range, range);
            for (int i = 0; i < range; i++) {
                out.println(Double.longBitsToDouble(lb + i));
            }
            return;
        }
        long[] s = new long[n];
        boolean check = true;
        if (range < Integer.MAX_VALUE) {
            int m = (int) range;
            if (range < 20L * n) {
                // Avoid duplicates
                check = false;
                int[] natural = IntStream.range(0, m).toArray();
                // Partial Fisher-Yates shuffle
                for (int i = 0; i < n; i++) {
                    int k = m - i - 1;
                    int j = rng.nextInt(m - k);
                    int t = natural[j];
                    natural[j] = natural[k];
                    // No required: natural[k] = t;
                    s[i] = lb + t;
                }
            } else {
                for (int i = 0; i < n; i++) {
                    s[i] = lb + rng.nextInt(m);
                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                s[i] = lb + rng.nextLong(range);
            }
        }
        Arrays.sort(s);
        if (check) {
            // Eliminate duplicates
            n = 0;
            long last = -1;
            for (int i = 0; i < s.length; i++) {
                if (last == s[i]) {
                    continue;
                }
                last = s[i];
                s[n++] = last;
            }
        }
        out.printf("# [%s, %s) %d / %d%n", lo, hi, n, range);
        for (int i = 0; i < n; i++) {
            out.println(Double.longBitsToDouble(s[i]));
        }
    }
}
