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
    /** Seed for data generation. */
    private static final long SEED = 6516587940839803692L;
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
        // Hurwitz zeta is accurate for all s > 1 including s -> 1.
        // Require by passing public API which can call the zeta function when a==1.
        HURWITZ_ZETA_1_32(s -> HurwitzZeta.zetaImp(s, 1), "zeta_1_32.csv", 1.62, 0.43),
        HURWITZ_ZETA_ABOVE_1(s -> HurwitzZeta.zetaImp(s, 1), "zeta_above1.csv", 1.7, 0.52),
//        // Require by-passing s <= 1
//        HURWITZ_ZETA_BELOW_1(s -> HurwitzZeta.zetaImp(s, 1), "zeta_below1.csv", 2.13, 0.63),
//        HURWITZ_ZETA_0_1(s -> HurwitzZeta.zetaImp(s, 1), "zeta_0_1.csv", 27, 4.3),
        // Borwein zeta has no support for negative s using reflection.
        // It is worse than BoostZeta
        BORWEIN_ZETA_1_32(RiemannZetaTest::borweinZeta, "zeta_1_32.csv", 2.7, 0.4),
        BORWEIN_ZETA_ABOVE_1(RiemannZetaTest::borweinZeta, "zeta_above1.csv", 2.7, 0.71),
        BORWEIN_ZETA_BELOW_1(RiemannZetaTest::borweinZeta, "zeta_below1.csv", 3.35, 0.72),
        BORWEIN_ZETA_0_1(RiemannZetaTest::borweinZeta, "zeta_0_1.csv", 6, 1.3),
        // BoostZeta is better than Hurwitz zeta for the domain s in (1, 32).
        // The method has increasing error with larger negative s.
        ZETA_1_32(RiemannZeta::value, "zeta_1_32.csv", 1.5, 0.27),
        ZETA_ABOVE_1(RiemannZeta::value, "zeta_above1.csv", 1.0, 0.34),
        ZETA_BELOW_1(RiemannZeta::value, "zeta_below1.csv", 1.48, 0.46),
        ZETA_0_1(RiemannZeta::value, "zeta_0_1.csv", 2.32, 0.69),
        ZETA_N_0_1(RiemannZeta::value, "zeta_N_0_1.csv", 6, 1.5),
        ZETA_N_1_4(RiemannZeta::value, "zeta_N_1_4.csv", 7, 1.6),
        ZETA_N_4_16(RiemannZeta::value, "zeta_N_4_16.csv", 23, 3.6),
        ZETA_N_16_64(RiemannZeta::value, "zeta_N_16_64.csv", 170, 15.5);

//        JDK Temurin 25.492-b09
//        HURWITZ_ZETA_1_32                     max    1.60678   RMS   0.414903   mean     0.00784791  n 3000
//        HURWITZ_ZETA_ABOVE_1                  max    1.65995   RMS   0.491651   mean     -0.0154402  n 2606
//        HURWITZ_ZETA_BELOW_1                  max    2.10653   RMS   0.613495   mean     -0.0393119  n 2535
//        HURWITZ_ZETA_0_1                      max    26.6667   RMS    4.18436   mean      0.0446286  n 2000
//        BORWEIN_ZETA_1_32                     max    2.56933   RMS   0.363840   mean      0.0850378  n 3000
//        BORWEIN_ZETA_ABOVE_1                  max    2.62981   RMS   0.686734   mean       0.354976  n 2606
//        BORWEIN_ZETA_BELOW_1                  max    3.31413   RMS   0.698234   mean      -0.334494  n 2535
//        BORWEIN_ZETA_0_1                      max    5.18607   RMS    1.16183   mean      -0.429437  n 2000
//        ZETA_1_32                             max    1.44796   RMS   0.258206   mean     0.00387724  n 3000
//        ZETA_ABOVE_1                          max   0.978747   RMS   0.329409   mean     0.00618317  n 2606
//        ZETA_BELOW_1                          max    1.46877   RMS   0.449973   mean    -0.00828687  n 2535
//        ZETA_0_1                              max    2.24354   RMS   0.676258   mean      0.0551843  n 2000
//        ZETA_N_0_1                            max    5.66801   RMS    1.43564   mean      -0.249907  n 2000
//        ZETA_N_1_4                            max    6.75694   RMS    1.52330   mean       0.154205  n 3000
//        ZETA_N_4_16                           max    22.1416   RMS    3.47762   mean       0.301526  n 3000
//        ZETA_N_16_64                          max    165.252   RMS    14.9244   mean       0.301025  n 3000

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

    @ParameterizedTest
    @CsvSource({
        // Created using mpmath (1.14.1) zeta
        // from mpmath import zeta, mp
        // mp.pretty = True
        // for i in range(1, 131): # double limit of Bernoulli B_2n
        //   print(f'"{-2*i+1}, {zeta(-2*i+1)}",')
        "-1, -0.0833333333333333333333333333333",
        "-3, 0.00833333333333333333333333333333",
        "-5, -0.00396825396825396825396825396825",
        "-7, 0.00416666666666666666666666666667",
        "-9, -0.00757575757575757575757575757576",
        "-11, 0.0210927960927960927960927960928",
        "-13, -0.0833333333333333333333333333333",
        "-15, 0.443259803921568627450980392157",
        "-17, -3.05395433027011974380395433027",
        "-19, 26.4562121212121212121212121212",
        "-21, -281.460144927536231884057971014",
        "-23, 3607.5105463980463980463980464",
        "-25, -54827.5833333333333333333333333",
        "-27, 974936.82385057471264367816092",
        "-29, -20052695.7966880789461434622725",
        "-31, 472384867.721629901960784313725",
        "-33, -12635724795.9166666666666666667",
        "-35, 380879311252.453688115530220793",
        "-37, -12850850499305.0833333333333333",
        "-39, 482414483548501.703715816703622",
        "-41, -20040310656516252.7381084216632",
        "-43, 916774360319533077.569927536232",
        "-45, -45979888343656503490.4379432624",
        "-47, 2518047192145109569708.90233202",
        "-49, -150017334921539287337114.401515",
        "-51, 9689957887463594065649794.28946",
        "-53, -676458823792928209909452423.018",
        "-55, 50890659468662289689766332915.9",
        "-57, -4.11472887925579786976654860676e+30",
        "-59, 3.56665820953755561096845746087e+32",
        "-61, -3.30660898765775767256802146704e+34",
        "-63, 3.27156342364787162642112270157e+36",
        "-65, -3.447378255827805387825645508e+38",
        "-67, 3.86142798327052588930927202002e+40",
        "-69, -4.58929744324543321688639890061e+42",
        "-71, 5.77753863427704318248848256879e+44",
        "-73, -7.69198587595071351674100759718e+46",
        "-75, 1.08136354499716546963540333511e+49",
        "-77, -1.60293645220089654060671023458e+51",
        "-79, 2.50194790415604628436566614985e+53",
        "-81, -4.10670523358102124797520450041e+55",
        "-83, 7.07987744084945806174529724334e+57",
        "-85, -1.28045468879395087901908497563e+60",
        "-87, 2.42673403923335240780208920671e+62",
        "-89, -4.8143218874045769355129570066e+64",
        "-91, 9.98755741757275306806527774082e+66",
        "-93, -2.16456348684351856313351361598e+69",
        "-95, 4.89623270396205532068492245156e+71",
        "-97, -1.15490239239635196639542716916e+74",
        "-99, 2.83822495706937069592641563365e+76",
        "-101, -7.26120088036067163036772815107e+78",
        "-103, 1.93235142334198120033323266084e+81",
        "-105, -5.34501604252886240053956094628e+83",
        "-107, 1.53560288464224230702071420133e+86",
        "-109, -4.57898726822657976538994744683e+88",
        "-111, 1.41620252121948092583601799759e+91",
        "-113, -4.54006522960926552491870532303e+93",
        "-115, 1.50766567588078597755948498945e+96",
        "-117, -5.18309491482645637761224790375e+98",
        "-119, 1.84356474272565291185736028806e+101",
        "-121, -6.78055547530909588969025102131e+103",
        "-123, 2.577332670275460450289647933e+106",
        "-125, -1.0119112875704597605007955796e+109",
        "-127, 4.10163461615422921089084567379e+111",
        "-129, -1.71552445340320193922071524606e+114",
        "-131, 7.40034257052690942716920556053e+116",
        "-133, -3.29092253570544434867706140857e+119",
        "-135, 1.50798315341647712056833365644e+122",
        "-137, -7.11698791882545486286760649291e+124",
        "-139, 3.45804291415777717919922833643e+127",
        "-141, -1.72909076066767483167489207071e+130",
        "-143, 8.89369916950329690887674533236e+132",
        "-145, -4.7038470619636014515138279862e+135",
        "-147, 2.55719382310602058749858077138e+138",
        "-149, -1.42840675004435277005808820901e+141",
        "-151, 8.1952152218313782940918703695e+143",
        "-153, -4.82764854227273717816101742819e+146",
        "-155, 2.91896123747703236500405982201e+149",
        "-157, -1.81089321625689040160530678804e+152",
        "-159, 1.15235772200211685798051266585e+155",
        "-161, -7.51923119519817697500081265839e+157",
        "-163, 5.02940165764110497246840522742e+160",
        "-165, -3.4473420444477676704609427599e+163",
        "-167, 2.42074586458685147183142674899e+166",
        "-169, -1.74094659203776765075736879892e+169",
        "-171, 1.28194898634822427378088228066e+172",
        "-173, -9.66241211085609184243169684778e+174",
        "-175, 7.45269103043008957309390945203e+177",
        "-177, -5.8808393311674371248220704455e+180",
        "-179, 4.74627186549076153992212525722e+183",
        "-181, -3.91691325947728254682903333391e+186",
        "-183, 3.30450714432260322283069086243e+189",
        "-185, -2.84928905509945827581152102174e+192",
        "-187, 2.51033293450775865129599057986e+195",
        "-189, -2.25939019954752532049562261337e+198",
        "-191, 2.07691380042876080434623778929e+201",
        "-193, -1.94947321749272591308734114001e+204",
        "-195, 1.86807314712659138998396980689e+207",
        "-197, -1.8270752662814576943866394409e+210",
        "-199, 1.82353863225956771810691544328e+213",
        "-201, -1.85686908101259450981907133715e+216",
        "-203, 1.92871898511956020928868278693e+219",
        "-205, -2.04311704602864475750762704423e+222",
        "-207, 2.20684116445278455076828336814e+225",
        "-209, -2.4300821796490274251390389767e+228",
        "-211, 2.72748878790834695290272297756e+231",
        "-213, -3.11974215737550845945157847856e+234",
        "-215, 3.63589387242826001493479749833e+237",
        "-217, -4.31683000307608832681396003788e+240",
        "-219, 5.22042448793871999720448178213e+243",
        "-221, -6.42926069497693048519893333726e+246",
        "-223, 8.06230338701308438135204143382e+249",
        "-225, -1.02927147379030111575795568666e+253",
        "-227, 1.33753296997805240211378429501e+256",
        "-229, -1.76894809027973797575657445272e+259",
        "-231, 2.38064790180923972522551743144e+262",
        "-233, -3.2597127947194184823502123513e+265",
        "-235, 4.54049623716012131918616627123e+268",
        "-237, -6.43285751931478506106894605686e+271",
        "-239, 9.26870486757493111152509786938e+274",
        "-241, -1.35796195002851814738921378693e+278",
        "-243, 2.02278397360493216811767187783e+281",
        "-245, -3.06299069922083360655392703162e+284",
        "-247, 4.71430853007426521282616631982e+287",
        "-249, -7.37410458713557576506584806391e+290",
        "-251, 1.17209627670508265765085284663e+294",
        "-253, -1.89288666446856573885385316552e+297",
        "-255, 3.10555175960489268960251418622e+300",
        "-257, -5.1754977470366797965163888379e+303",
        // B_2n for n=130 is not a double
        // "-259, 8.76015634462292151490407301349e+306",
    })
    void testZetaNegativeOddInteger(int s, double z) {
        // Uses the Bernoulli numbers divided by an integer
        assertClose(RiemannZeta::value, s, z, 1);
        // Check the function is monotonic in s (relevant when precomputed values are used).
        // This works for close s but not next s
        double zl;
        double zu;
        // ULP  Fails
        // 3    -7
        // 2    -5, -7, -225, -227, -251
        // 3    16 cases
        zl = RiemannZeta.value(Double.longBitsToDouble(Double.doubleToRawLongBits(s) + 4));
        zu = RiemannZeta.value(Double.longBitsToDouble(Double.doubleToRawLongBits(s) - 4));
        // Values alternate sign so test [zl, zu] contains z
        if (z <= zu) {
            Assertions.assertTrue(z >= zl, () -> String.format("%s, %s, %s", zl, z, zu));
        } else {
            Assertions.assertTrue(z <= zl, () -> String.format("%s, %s, %s", zl, z, zu));
        }
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
            Arguments.of(1e+19, 1.0, 0),

            // First negative odd integer where the Bernoulli number B_2n is not a double
            Arguments.of(-259, 8.76015634462292151490407301349e+306, 230),
            // Large integer that is not an int:
            //Arguments.of(-2147483649.0, -9.10272023221439406602025581451e+17393448289, 0)
            Arguments.of(-2147483649.0, Double.NEGATIVE_INFINITY, 0)
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
            System.out.printf("JDK %s %s%n",
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
     */
    @Test
    @Disabled("Used to generate test data")
    void testDataSample0to1() throws IOException {
        final double[] sample = new SplittableRandom(SEED).doubles(2000).toArray();
        Arrays.sort(sample);
        try (PrintStream out = getPrintStream("zeta_0_1.txt")) {
            out.printf("# Dyadic doubles in [0, 1) : [%s, %s] n=%d%n",
                sample[0], sample[sample.length - 1], sample.length);
            for (final double s : sample) {
                out.println(s);
            }
        }
    }

    /**
     * Create test data for close to and above 1.
     */
    @Test
    @Disabled("Used to generate test data")
    void testDataSampleCloseAbove1() throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Create node points using 1 + 2^-b
        final double[] nodes = new double[42];
        double b = 0x1p-52;
        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = 1 + b;
            b *= 2;
        }

        final int size = 3000 / (nodes.length - 1);
        try (PrintStream out = getPrintStream("zeta_above1.txt")) {
            for (int j = 0; j < nodes.length - 1; j++) {
                sample(rng, size, nodes[j], nodes[j + 1], out);
            }
        }
    }

    /**
     * Create test data for close to and below 1.
     */
    @Test
    @Disabled("Used to generate test data")
    void testDataSampleCloseBelow1() throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Create node points using 1 - 2^-b
        final double[] nodes = new double[42];
        double b = 0x1p-53;
        // Nodes must be ascending magnitude so fill from the end.
        // This is an exclusive upper bound.
        nodes[nodes.length - 1] = 1;
        for (int i = 1; i < nodes.length; i++) {
            nodes[nodes.length - i - 1] = 1 - b;
            b *= 2;
        }

        final int size = 3000 / (nodes.length - 1);
        try (PrintStream out = getPrintStream("zeta_below1.txt")) {
            for (int j = 0; j < nodes.length - 1; j++) {
                sample(rng, size, nodes[j], nodes[j + 1], out);
            }
        }
    }

    /**
     * Create test data for {@code s in [2^lb, 2^ub)}.
     */
    @ParameterizedTest
    @CsvSource({
        // For negative s in [1, 4); [4, 16); [16; 64)
        "0, 2",
        "2, 4",
        "4, 6",
        // For positive s in [1, 32)
        // Note that the zeta function is easily
        // computed for {@code s > 32} using 3^-s + 2^-s + 1.
        "0, 5",
    })
    @Disabled("Used to generate test data")
    void testDataSamplePowerOf2(int lb, int ub) throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Create node points using 2^b
        final double[] nodes = new double[ub - lb + 1];
        double b = Math.scalb(1, lb);
        final int lower = (int) b;
        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = b;
            b *= 2;
        }
        final int upper = (int) (b / 2);

        final int size = 3000 / (nodes.length - 1);
        try (PrintStream out = getPrintStream(String.format("zeta_%d_%d.txt", lower, upper))) {
            for (int j = 0; j < nodes.length - 1; j++) {
                sample(rng, size, nodes[j], nodes[j + 1], out);
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
     * @param hi the high bound (exclusive)
     */
    private static void sample(SplittableRandom rng, int n, double lo, double hi,
            PrintStream out) {
        assert lo >= 0;
        assert hi > lo;
        final long lb = Double.doubleToRawLongBits(lo);
        final long hb = Double.doubleToRawLongBits(hi);
        final long range = hb - lb;
        if (range < n) {
            out.printf("# [%s, %s) %d / %d%n", lo, hi, range, range);
            for (int i = 0; i < range; i++) {
                out.println(Double.longBitsToDouble(lb + i));
            }
            return;
        }
        final long[] s = new long[n];
        boolean check = true;
        if (range < Integer.MAX_VALUE) {
            final int m = (int) range;
            if (range < 20L * n) {
                // Avoid duplicates when sampling a small range
                check = false;
                final int[] natural = IntStream.range(0, m).toArray();
                // Partial Fisher-Yates shuffle
                for (int i = 0; i < n; i++) {
                    // Index into the remaining array length
                    final int k = m - i - 1;
                    // Swap index k with any position down to 0 (including itself)
                    final int j = rng.nextInt(k + 1);
                    final int t = natural[j];
                    natural[j] = natural[k];
                    // No required: natural[k] = t
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

    /**
     * Gets the prints the stream.
     * Adds a header line to the output indicating how the data was created.
     *
     * @param filename the filename
     * @return the stream
     */
    private PrintStream getPrintStream(String filename) throws IOException {
        final PrintStream out = new PrintStream(Files.newOutputStream(Paths.get("target", filename)));
        out.printf("# Generated by %s%n", getClass().getSimpleName());
        return out;
    }
}
