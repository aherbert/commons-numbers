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
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.MethodOrderer;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestMethodOrder;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Test the {@link RiemannZeta} function.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class RiemannZetaTest {
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
        // s in (1, 32)
        BORWEIN_ZETA_ABOVE1(RiemannZetaTest::borweinZeta, "zeta_above1.csv", 2.8, 0.66),
        HURWITZ_ZETA_ABOVE1(s -> HurwitzZeta.value(s, 1), "zeta_above1.csv", 1.68, 0.5),
        ZETA_ABOVE1(RiemannZeta::value, "zeta_above1.csv", 1.46, 0.34),
        BORWEIN_ZETA_0_1(RiemannZetaTest::borweinZeta, "zeta_0_1.csv", 4.05, 1.28),
        HURWITZ_ZETA_0_1(s -> HurwitzZeta.value(s, 1), "zeta_0_1.csv", 19.5, 4.10),
        ZETA_0_1(RiemannZeta::value, "zeta_0_1.csv", 1.79, 0.69),
        ;

     // Temurin 25.492-b09
//        BORWEIN_ZETA_ABOVE1                   max    2.76154   RMS   0.641891   mean       0.311163  n 8827
//        HURWITZ_ZETA_ABOVE1                   max    1.66474   RMS   0.486392   mean    -0.00900755  n 8827
//        ZETA_ABOVE1                           max    1.44866   RMS   0.329523   mean    -0.00187909  n 8827
//        BORWEIN_ZETA_0_1                      max    3.99744   RMS    1.26149   mean      -0.463086  n  500
//        HURWITZ_ZETA_0_1                      max    19.0470   RMS    4.07605   mean      0.0769531  n  500
//        ZETA_0_1                              max    1.77152   RMS   0.679988   mean      0.0238564  n  500

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
        }
        // Create the table
        final BigFraction dn = dk[n];
        for (int k = 0; k < n; k++) {
            final double sign = Math.pow(-1, k);
            final double dndk = sign * dn.subtract(dk[k]).divide(dn).doubleValue();
            Assertions.assertEquals(DN_DK[k], dndk);
        }
    }


    @ParameterizedTest
    @CsvSource({
        "1.0, Infinity",
        "-262384234, 0.0",
        "Infinity, 1.0",
        "-Infinity, NaN",
        "NaN, NaN",
    })
    void testZetaSpecial(double s, double z) {
        Assertions.assertEquals(z, RiemannZeta.value(s));
    }

    @ParameterizedTest
    @CsvSource({
        // Created using mpmath (1.14.1) zeta
        // from mpmath import zeta, mp
        // mp.pretty = True
        // for i in range(0, 28):
        //   print(f'"{2*i}, {zeta(2*i)}",')
        "0, -0.5",
        "2, 1.6449340668482264364724151666460251892189499012068",
        "4, 1.0823232337111381915160036965411679027747509519187",
        "6, 1.0173430619844491397145179297909205279018174900329",
        "8, 1.0040773561979443393786852385086524652589607906499",
        "10, 1.0009945751278180853371459589003190170060195315645",
        "12, 1.000246086553308048298637998047739670960416088458",
        "14, 1.0000612481350587048292585451051353337474816961692",
        "16, 1.0000152822594086518717325714876367220232373889905",
        "18, 1.000003817293264999839856461644621939730454697219",
        "20, 1.0000009539620338727961131520386834493459437941874",
        "22, 1.0000002384505027277329900036481867529949350418218",
        "24, 1.0000000596081890512594796124402079358012275039188",
        "26, 1.0000000149015548283650412346585066306986288647882",
        "28, 1.0000000037253340247884570548192040184024232328931",
        "30, 1.000000000931327432419668182871764735021219813568",
        "32, 1.000000000232831183367650549200145597594049502483",
        "34, 1.0000000000582077208790270088924368598910630541731",
        "36, 1.0000000000145519218910419842359296322453184209838",
        "38, 1.0000000000036379795473786511902372363558732735126",
        "40, 1.0000000000009094947840263889282533118386949087539",
        "42, 1.0000000000002273736845824652515226821577978691214",
        "44, 1.0000000000000568434198762758560927718296752406855",
        "46, 1.0000000000000142108548280316067698343071417395377",
        "48, 1.000000000000003552713691337113673298469534059343",
        "50, 1.0000000000000008881784210930815903096091386391386",
        "52, 1.0000000000000002220446050798041983999320094204654",
        // Effectively 1.0
        "54, 1.0000000000000000555111512484548124372373659050943",
    })
    void testZetaEvenInteger(int s, double z) {
        // Exact except at s=2
        assertClose(RiemannZeta::value, s, z, s == 2 ? 1 : 0);
    }

    @ParameterizedTest
    @CsvSource({
        // Created using mpmath (1.14.1) zeta
        // from mpmath import zeta, mp
        // mp.pretty = True
        // for i in range(1, 28):
        //   print(f'"{2*i+1}, {zeta(2*i+1)}",')
        "3, 1.2020569031595942853997381615114499907649862923405",
        "5, 1.0369277551433699263313654864570341680570809195019",
        "7, 1.0083492773819228268397975498497967595998635605652",
        "9, 1.0020083928260822144178527692324120604856058513949",
        "11, 1.0004941886041194645587022825264699364686064357582",
        "13, 1.0001227133475784891467518365263573957142751058955",
        "15, 1.0000305882363070204935517285106450625876279487069",
        "17, 1.0000076371976378997622736002935630292130882490903",
        "19, 1.0000019082127165539389256569577951013532585711448",
        "21, 1.0000004769329867878064631167196043730459664466948",
        "23, 1.0000001192199259653110730677887188823263872549978",
        "25, 1.0000000298035035146522801860637050693660118447309",
        "27, 1.000000007450711789835429491981004170604119454719",
        "29, 1.0000000018626597235130490064039099454169480616653",
        "31, 1.0000000004656629065033784072989233251220071062692",
        "33, 1.0000000001164155017270051977592973835456309516522",
        "35, 1.000000000029103850444970996869294252278840464107",
        "37, 1.0000000000072759598350574810145208690123380592649",
        "39, 1.0000000000018189896503070659475848321007300850306",
        "41, 1.0000000000004547473783042154026799112029488570339",
        "43, 1.0000000000001136868407680227849349104838025906437",
        "45, 1.0000000000000284217097688930185545507370494266207",
        "47, 1.0000000000000071054273952108527128773544799568",
        "49, 1.0000000000000017763568435791203274733490144002796",
        "51, 1.0000000000000004440892103143813364197770940268121",
        "53, 1.0000000000000001110223025141066133720544569921383",
        // Effectively 1.0
        "55, 1.0000000000000000277555756213612417258163245385407",
    })
    void testZetaOddInteger(int s, double z) {
        // Boost uses pre-computed values
        // "as these are of great benefit to some infinite series calculations".
        // Check this is exact.
        assertClose(RiemannZeta::value, s, z, 0);
        // Check the function is monotonic in s (relevant when precomputed values are used).
        Assertions.assertTrue(z >= RiemannZeta.value(Math.nextUp((double) s)));
        Assertions.assertTrue(z <= RiemannZeta.value(Math.nextDown((double) s)));
    }

    @ParameterizedTest
    @ValueSource(doubles = {-2, -4, -6, -235648, Integer.MIN_VALUE, Integer.MIN_VALUE - 2L,
        // Anything below -(2^53) is even except -infinity
        -0x1p53, -0x1p53 - 1, -0x1p53 - 2,
        Long.MIN_VALUE, -Double.MAX_VALUE})
    void testZetaNegativeEvenInteger(double s) {
        Assertions.assertEquals(0.0, RiemannZeta.value(s));
    }

    /**
     * Spot tests for the zeta function to check various points in the domain and extreme values.
     */
    @ParameterizedTest
    @MethodSource
    void testZetaSpot(double s, double z, int ulp) {
        assertClose(RiemannZeta::value, s, z, ulp);
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
            Arguments.of(53.0000000001, 1.00000000000000011102230250641, 0),
            // Effectively 1.0
            Arguments.of(53.000000001, 1.00000000000000011102230243715, 0),
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
        // CHECKSTYLE: resume regexp
    }

    /**
     * Create test data for {@code s in [0, 1)}. Note that the zeta function is easily
     * computed for {@code s -> 0} so using dyadic rationals in [0, 1) is sufficient.
     *
     * @throws IOException Signals that an I/O exception has occurred.
     */
    @Test
    @Disabled("Used to generate test data")
    void testSample01() throws IOException {
        final double[] sample = new SplittableRandom().doubles(500).toArray();
        Arrays.sort(sample);
        try (PrintStream out = new PrintStream(Files.newOutputStream(Paths.get("target", "zeta_0_1.txt")))) {
            out.printf("# Dyadic doubles in [0, 1) : [%s, %s] n=%d%n",
                sample[0], sample[sample.length - 1], sample.length);
            for (double s : sample) {
                out.println(s);
            }
        }
    }

    /**
     * Create test data for {@code s in (1, 32)}. Note that the zeta function is easily
     * computed for {@code s > 32} using 3^-s + 2^-s + 1.
     *
     * @throws IOException Signals that an I/O exception has occurred.
     */
    @Test
    @Disabled("Used to generate test data")
    void testSampleAbove1() throws IOException {
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

        try (PrintStream out = new PrintStream(Files.newOutputStream(Paths.get("target", "zeta_above1.txt")))) {
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
