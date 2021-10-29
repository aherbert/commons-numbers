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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;
import java.util.regex.Pattern;
import org.apache.commons.numbers.gamma.BoostGamma.Lanczos;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.EnumSource.Mode;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Tests for {@link BoostGamma}. Special functions from {@link BoostMath} and {@link SpecialMath}
 * are also tested as these are used within the {@link BoostGamma} class.
 *
 * <p>Note: Some resource data files used in these tests have been extracted
 * from the Boost test files for the gamma functions.
 */
class BoostGammaTest {
    /** All representable factorials. */
    private static final double[] FACTORIAL = BoostGamma.getFactorials();
    /** Value for the sqrt of the epsilon for relative error.
     * This is equal to the Boost constant {@code boost::math::tools::root_epsilon<double>()}. */
    private static final double ROOT_EPSILON = 1.4901161193847656E-8;
    /** Approximate value for ln(Double.MAX_VALUE).
     * This is equal to the Boost constant {@code boost::math::tools::log_max_value<double>()}. */
    private static final int LOG_MAX_VALUE = 709;
    /** The largest factorial that can be represented as a double.
     * This is equal to the Boost constant {@code boost::math::max_factorial<double>::value}. */
    private static final int MAX_FACTORIAL = 170;
    /** Euler's constant. */
    private static final double EULER = 0.5772156649015328606065120900824024310;

    /**
     * Define the expected error for a test.
     */
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
     * <h2>Note on accuracy</h2>
     *
     * <p>The Boost functions use the default policy of internal promotion
     * of double to long double if it offers more precision. Code comments
     * in the implementations for the maximum error are using the defaults with
     * promotion enabled where the error is 'effectively zero'. Java does not
     * support long double computation. Tolerances have been set to allow tests to
     * pass. Spot checks on larger errors have been verified against the reference
     * implementation compiled with promotion of double <strong>disabled</strong>.
     *
     * @see <a href="https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/relative_error.html">Relative error</a>
     * @see <a href="https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/policy_tut_defaults.html">Policy defaults</a>
     */
    private enum TestCase implements TestError {
        // Note:
        // The original Boost tgamma function is not as accurate as the
        // NSWC Library of Mathematics Subroutines in the range [-20, 20].
        // The default implementation uses long double to compute and
        // performs a narrowing cast to the double result.
        // The code is here for testing.

        /** gamma Boost near 0 data. */
        TGAMMAO_NEAR_0(BoostGammaTest::tgammaOriginal, "gamma_near_0_data.csv", 3.3, 1.25),
        /** gamma Boost near 1 data. */
        TGAMMAO_NEAR_1(BoostGammaTest::tgammaOriginal, "gamma_near_1_data.csv", 3.3, 1.25),
        /** gamma Boost near 2 data. */
        TGAMMAO_NEAR_2(BoostGammaTest::tgammaOriginal, "gamma_near_2_data.csv", 3, 1.25),
        /** gamma Boost near -10 data. */
        TGAMMAO_NEAR_M10(BoostGammaTest::tgammaOriginal, "gamma_near_m10_data.csv", 2.5, 1.25),
        /** gamma -20 to 0 data. */
        TGAMMAO_M20_0(BoostGammaTest::tgammaOriginal, "gamma_m20_0_data.csv", 4.5, 1.5),
        /** gamma 0 to 20 data. */
        TGAMMAO_0_20(BoostGammaTest::tgammaOriginal, "gamma_0_20_data.csv", 3.25, 1.25),
        /** gamma very near 0 data. */
        TGAMMAO_VERY_NEAR_0(BoostGamma::tgamma, "gamma_very_near_0_data.csv", 4, 0.75),

        /** gamma Boost factorial data. */
        TGAMMA_FACTORIALS(BoostGamma::tgamma, "gamma_factorials_data.csv", 2.5, 0.99),
        /** gamma Boost near 0 data. */
        TGAMMA_NEAR_0(BoostGamma::tgamma, "gamma_near_0_data.csv", 1.6, 0.75),
        /** gamma Boost near 1 data. */
        TGAMMA_NEAR_1(BoostGamma::tgamma, "gamma_near_1_data.csv", 1.5, 0.75),
        /** gamma Boost near 2 data. */
        TGAMMA_NEAR_2(BoostGamma::tgamma, "gamma_near_2_data.csv", 1.1, 0.75),
        /** gamma Boost near -10 data. */
        TGAMMA_NEAR_M10(BoostGamma::tgamma, "gamma_near_m10_data.csv", 1.8, 0.75),
        /** gamma Boost near -55 data. */
        TGAMMA_NEAR_M55(BoostGamma::tgamma, "gamma_near_m55_data.csv", 2.5, 1.25),
        /** gamma -20 to 0 data. */
        TGAMMA_M20_0(BoostGamma::tgamma, "gamma_m20_0_data.csv", 3, 0.75),
        /** gamma 0 to 20 data. */
        TGAMMA_0_20(BoostGamma::tgamma, "gamma_0_20_data.csv", 3, 0.75),
        /** gamma 20 to 150 data. */
        TGAMMA_20_150(BoostGamma::tgamma, "gamma_20_150_data.csv", 4, 1.25),
        /** gamma very near 0 data. */
        TGAMMA_VERY_NEAR_0(BoostGamma::tgamma, "gamma_very_near_0_data.csv", 4, 0.75),

        /** gamma Boost factorial data. */
        LGAMMA_FACTORIALS(BoostGamma::lgamma, "gamma_factorials_data.csv", 2, 1.5, 0.3),
        /** gamma Boost near 0 data. */
        LGAMMA_NEAR_0(BoostGamma::lgamma, "gamma_near_0_data.csv", 2, 1.25, 0.5),
        /** gamma Boost near 1 data. */
        LGAMMA_NEAR_1(BoostGamma::lgamma, "gamma_near_1_data.csv", 2, 1.5, 0.75),
        /** gamma Boost near 2 data. */
        LGAMMA_NEAR_2(BoostGamma::lgamma, "gamma_near_2_data.csv", 2, 0.9, 0.25),
        /** gamma Boost near -10 data. */
        // Negative z is lower precision as logs are used with values approaching 1
        LGAMMA_NEAR_M10(BoostGamma::lgamma, "gamma_near_m10_data.csv", 2, 8, 2.5),
        /** gamma Boost near -55 data. */
        LGAMMA_NEAR_M55(BoostGamma::lgamma, "gamma_near_m55_data.csv", 2, 1.5, 0.75),
        /** gamma -20 to 0 data. */
        // The value -2.75 is low precision
        LGAMMA_M20_0(BoostGamma::lgamma, "gamma_m20_0_data.csv", 2, 100, 10),
        /** gamma 0 to 20 data. */
        LGAMMA_0_20(BoostGamma::lgamma, "gamma_0_20_data.csv", 2, 1.8, 0.75),
        /** gamma 20 to 150 data. */
        LGAMMA_20_150(BoostGamma::lgamma, "gamma_20_150_data.csv", 2, 1.5, 0.5),
        /** gamma very near 0 data. */
        LGAMMA_VERY_NEAR_0(BoostGamma::lgamma, "gamma_very_near_0_data.csv", 2, 1.8, 0.5),

        /** gamma(1+x) - 1 Boost data. */
        TGAMMAP1M1(BoostGamma::tgamma1pm1, "gamma1pm1_data.csv", 1.8, 0.75),

        /** log(1+x) - 1  data. */
        LOG1PMX(SpecialMath::log1pmx, "log1pmx_data.csv", 0.9, 0.15);

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

    /**
     * Define the test cases for each resource file for two argument functions.
     * This encapsulates the function to test, the expected maximum and RMS error, and
     * the resource file containing the data.
     */
    private enum BiTestCase implements TestError {
        /** pow(x, y) - 1 Boost data. */
        POWM1(BoostMath::powm1, "powm1_data.csv", 2.5, 0.4),
        /** igamma Boost int data. */
        IGAMMA_UPPER_INT(BoostGamma::tgamma, "igamma_int_data.csv", 10, 2),
        /** igamma Boost small data. */
        IGAMMA_UPPER_SMALL(BoostGamma::tgamma, "igamma_small_data.csv", 4.5, 0.9),
        /** igamma Boost med data. */
        IGAMMA_UPPER_MED(BoostGamma::tgamma, "igamma_med_data.csv", 16, 2.8),
        /** igamma Boost big data. */
        IGAMMA_UPPER_BIG(BoostGamma::tgamma, "igamma_big_data.csv", 8, 1.5),
        /** igamma extra data containing edge cases. */
        IGAMMA_UPPER_EXTRA(BoostGamma::tgamma, "igamma_extra_data.csv", 15, 4),
        /** igamma Boost int data. */
        IGAMMA_Q_INT(BoostGamma::gammaQ, "igamma_int_data.csv", 3, 14, 2.4),
        /** igamma Boost small data. */
        IGAMMA_Q_SMALL(BoostGamma::gammaQ, "igamma_small_data.csv", 3, 4, 1.1),
        /** igamma Boost med data. */
        IGAMMA_Q_MED(BoostGamma::gammaQ, "igamma_med_data.csv", 3, 43, 6),
        /** igamma Boost big data. */
        IGAMMA_Q_BIG(BoostGamma::gammaQ, "igamma_big_data.csv", 3, 550, 62),
        /** igamma extra data containing edge cases. */
        IGAMMA_Q_EXTRA(BoostGamma::gammaQ, "igamma_extra_data.csv", 3, 60, 20),
        /** igamma Boost int data. */
        IGAMMA_LOWER_INT(BoostGamma::tgammaLower, "igamma_int_data.csv", 4, 6, 1.5),
        /** igamma Boost small data. */
        IGAMMA_LOWER_SMALL(BoostGamma::tgammaLower, "igamma_small_data.csv", 4, 3, 0.75),
        /** igamma Boost med data. */
        IGAMMA_LOWER_MED(BoostGamma::tgammaLower, "igamma_med_data.csv", 4, 12, 2.4),
        /** igamma Boost big data. */
        IGAMMA_LOWER_BIG(BoostGamma::tgammaLower, "igamma_big_data.csv", 4, 10, 1.5),
        /** igamma extra data containing edge cases. */
        IGAMMA_LOWER_EXTRA(BoostGamma::tgammaLower, "igamma_extra_data.csv", 4, 5, 1.5),
        /** igamma Boost int data. */
        IGAMMA_P_INT(BoostGamma::gammaP, "igamma_int_data.csv", 5, 21, 4.6),
        /** igamma Boost small data. */
        IGAMMA_P_SMALL(BoostGamma::gammaP, "igamma_small_data.csv", 5, 3.3, 0.9),
        /** igamma Boost med data. */
        IGAMMA_P_MED(BoostGamma::gammaP, "igamma_med_data.csv", 5, 60, 11),
        /** igamma Boost big data. */
        IGAMMA_P_BIG(BoostGamma::gammaP, "igamma_big_data.csv", 5, 430, 55),
        /** igamma extra data containing edge cases. */
        IGAMMA_P_EXTRA(BoostGamma::gammaP, "igamma_extra_data.csv", 5, 4, 1.5);

        /** The function. */
        private final DoubleBinaryOperator fun;

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
        BiTestCase(DoubleBinaryOperator fun, String filename, double maxUlp, double rmsUlp) {
            this(fun, filename, 2, maxUlp, rmsUlp);
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
        BiTestCase(DoubleBinaryOperator fun, String filename, int expected, double maxUlp, double rmsUlp) {
            this.fun = fun;
            this.filename = filename;
            this.expected = expected;
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
        }

        /**
         * @return function to test
         */
        DoubleBinaryOperator getFunction() {
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

    @ParameterizedTest
    @CsvSource({
        // Pole errors
        "0, NaN",
        "-1, NaN",
        "-2, NaN",
        // Factorials: gamma(n+1) = n!
        "1, 1",
        "2, 1",
        "3, 2",
        "4, 6",
        "5, 24",
        "171, 0.7257415615307998967396728211129263114717e307",
        "172, Infinity",
    })
    void testTGammaEdgeCases(double z, double p) {
        Assertions.assertEquals(p, BoostGamma.tgamma(z));
    }

    /**
     * tgamma spot tests extracted from
     * {@code boost/libs/math/test/test_gamma.hpp}.
     */
    @Test
    void testTGammaSpotTests() {
        final int tolerance = 50;
        assertClose(BoostGamma::tgamma, 3.5, 3.3233509704478425511840640312646472177454052302295, tolerance);
        assertClose(BoostGamma::tgamma, 0.125, 7.5339415987976119046992298412151336246104195881491, tolerance);
        assertClose(BoostGamma::tgamma, -0.125, -8.7172188593831756100190140408231437691829605421405, tolerance);
        assertClose(BoostGamma::tgamma, -3.125, 1.1668538708507675587790157356605097019141636072094, tolerance);
        // Lower tolerance on this one, is only really needed on Linux x86 systems, result is mostly down to std lib accuracy:
        assertClose(BoostGamma::tgamma, -53249.0 / 1024, -1.2646559519067605488251406578743995122462767733517e-65, tolerance * 3);

        // Very small values, from a bug report by Rocco Romeo:
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -12), 4095.42302574977164107280305038926932586783813167844235368772, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -14), 16383.4228446989052821887834066513143241996925504706815681204, tolerance * 2);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -25), 3.35544314227843645746319656372890833248893111091576093784981e7, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -27), 1.34217727422784342467508497080056807355928046680073490038257e8, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -29), 5.36870911422784336940727488260481582524683632281496706906706e8, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -35), 3.43597383674227843351272524573929605605651956475300480712955e10, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -54), 1.80143985094819834227843350984671942971248427509141008005685e16, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -64), 1.84467440737095516154227843350984671394471047428598176073616e19, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -66), 7.37869762948382064634227843350984671394068921181531525785592922800e19, tolerance);
        assertClose(BoostGamma::tgamma, Math.scalb(1.0, -33), 8.58993459142278433521360841138215453639282914047157884932317481977e9, tolerance);
        assertClose(BoostGamma::tgamma, 4 / Double.MAX_VALUE, Double.MAX_VALUE / 4, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -12), -4096.57745718775464971331294488248972086965434176847741450728, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -14), -16384.5772760354695939336148831283410381037202353359487504624, tolerance * 2);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -25), -3.35544325772156943776992988569766723938420508937071533029983e7, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -27), -1.34217728577215672270574319043497450577151370942651414968627e8, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -29), -5.36870912577215666743793215770406791630514293641886249382012e8, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -34), -1.71798691845772156649591034966100693794360502123447124928244e10, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -54), -1.80143985094819845772156649015329155101490229157245556564920e16, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -64), -1.84467440737095516165772156649015328606601289230246224694513e19, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -66), -7.37869762948382064645772156649015328606199162983179574406439e19, tolerance);
        assertClose(BoostGamma::tgamma, -Math.scalb(1.0, -33), -8.58993459257721566501667413261977598620193488449233402857632e9, tolerance);
        assertClose(BoostGamma::tgamma, -4 / Double.MAX_VALUE, -Double.MAX_VALUE / 4, tolerance);
        assertClose(BoostGamma::tgamma, -1 + Math.scalb(1.0, -22), -4.19430442278467170746130758391572421252211886167956799318843e6, tolerance);
        assertClose(BoostGamma::tgamma, -1 - Math.scalb(1.0, -22), 4.19430357721600151046968956086404748206205391186399889108944e6, tolerance);
        assertClose(BoostGamma::tgamma, -4 + Math.scalb(1.0, -20), 43690.7294216755534842491085530510391932288379640970386378756, tolerance);
        assertClose(BoostGamma::tgamma, -4 - Math.scalb(1.0, -20), -43690.6039118698506165317137699180871126338425941292693705533, tolerance);
        assertClose(BoostGamma::tgamma, -1 + Math.scalb(1.0, -44), -1.75921860444164227843350985473932247549232492467032584051825e13, tolerance);
        assertClose(BoostGamma::tgamma, -1 - Math.scalb(1.0, -44), 1.75921860444155772156649016131144377791001546933519242218430e13, tolerance);
        assertClose(BoostGamma::tgamma, -4 + Math.scalb(1.0, -44), 7.33007751850729421569517998006564998020333048893618664936994e11, tolerance);
        assertClose(BoostGamma::tgamma, -4 - Math.scalb(1.0, -44), -7.33007751850603911763815347967171096249288790373790093559568e11, tolerance);
        // Test bug fixes in tgamma:
        assertClose(BoostGamma::tgamma, 142.75, 7.8029496083318133344429227511387928576820621466e244, tolerance * 4);
        assertClose(BoostGamma::tgamma, -Double.MIN_VALUE, Double.NEGATIVE_INFINITY, 0);
        assertClose(BoostGamma::tgamma, Double.MIN_VALUE, Double.POSITIVE_INFINITY, 0);
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {"TGAMMAO_.*"})
    void testTGammaOriginal(TestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {"TGAMMA_.*"})
    void testTGamma(TestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @CsvSource({
        // Pole errors
        "0, NaN",
        "-1, NaN",
        "-2, NaN",
        // Factorials: gamma(n+1) = n!
        "1, 0",
        "2, 0",
    })
    void testLGammaEdgeCases(double z, double p) {
        Assertions.assertEquals(p, BoostGamma.lgamma(z));
    }

    /**
     * lgamma spot tests extracted from
     * {@code boost/libs/math/test/test_gamma.hpp}.
     */
    @Test
    void testLGammaSpotTests() {
        final int tolerance = 1;
        final int[] sign = {0};
        final DoubleUnaryOperator fun = z -> BoostGamma.lgamma(z, sign);

        assertClose(fun, 3.5, 1.2009736023470742248160218814507129957702389154682, tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, 0.125, 2.0194183575537963453202905211670995899482809521344, tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -0.125, 2.1653002489051702517540619481440174064962195287626, tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -3.125, 0.1543111276840418242676072830970532952413339012367, tolerance * 2);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -53249.0 / 1024, -149.43323093420259741100038126078721302600128285894, tolerance);
        Assertions.assertEquals(-1, sign[0]);
        // Very small values, from a bug report by Rocco Romeo:
        assertClose(fun, Math.scalb(1.0, -12), Math.log(4095.42302574977164107280305038926932586783813167844235368772), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -14), Math.log(16383.4228446989052821887834066513143241996925504706815681204), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -25), Math.log(3.35544314227843645746319656372890833248893111091576093784981e7), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -27), Math.log(1.34217727422784342467508497080056807355928046680073490038257e8), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -29), Math.log(5.36870911422784336940727488260481582524683632281496706906706e8), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -35), Math.log(3.43597383674227843351272524573929605605651956475300480712955e10), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -54), Math.log(1.80143985094819834227843350984671942971248427509141008005685e16), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -64), Math.log(1.84467440737095516154227843350984671394471047428598176073616e19), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -66), Math.log(7.37869762948382064634227843350984671394068921181531525785592922800e19), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, Math.scalb(1.0, -33), Math.log(8.58993459142278433521360841138215453639282914047157884932317481977e9), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, 4 / Double.MAX_VALUE, Math.log(Double.MAX_VALUE / 4), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -12), Math.log(4096.57745718775464971331294488248972086965434176847741450728), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -14), Math.log(16384.5772760354695939336148831283410381037202353359487504624), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -25), Math.log(3.35544325772156943776992988569766723938420508937071533029983e7), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -27), Math.log(1.34217728577215672270574319043497450577151370942651414968627e8), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -29), Math.log(5.36870912577215666743793215770406791630514293641886249382012e8), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -34), Math.log(1.71798691845772156649591034966100693794360502123447124928244e10), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -54), Math.log(1.80143985094819845772156649015329155101490229157245556564920e16), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -64), Math.log(1.84467440737095516165772156649015328606601289230246224694513e19), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -66), Math.log(7.37869762948382064645772156649015328606199162983179574406439e19), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -Math.scalb(1.0, -33), Math.log(8.58993459257721566501667413261977598620193488449233402857632e9), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -4 / Double.MAX_VALUE, Math.log(Double.MAX_VALUE / 4), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -1 + Math.scalb(1.0, -22), Math.log(4.19430442278467170746130758391572421252211886167956799318843e6), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -1 - Math.scalb(1.0, -22), Math.log(4.19430357721600151046968956086404748206205391186399889108944e6), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -4 + Math.scalb(1.0, -20), Math.log(43690.7294216755534842491085530510391932288379640970386378756), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -4 - Math.scalb(1.0, -20), Math.log(43690.6039118698506165317137699180871126338425941292693705533), tolerance);
        Assertions.assertEquals(-1, sign[0]);

        assertClose(fun, -1 + Math.scalb(1.0, -44), Math.log(1.75921860444164227843350985473932247549232492467032584051825e13), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        assertClose(fun, -1 - Math.scalb(1.0, -44), Math.log(1.75921860444155772156649016131144377791001546933519242218430e13), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -4 + Math.scalb(1.0, -44), Math.log(7.33007751850729421569517998006564998020333048893618664936994e11), tolerance);
        Assertions.assertEquals(1, sign[0]);
        assertClose(fun, -4 - Math.scalb(1.0, -44), Math.log(7.33007751850603911763815347967171096249288790373790093559568e11), tolerance);
        Assertions.assertEquals(-1, sign[0]);
        //
        // Extra large values for lgamma, see https://github.com/boostorg/math/issues/242
        //
        assertClose(fun, Math.scalb(11103367432951928.0, 32), 2.7719825960021351251696385101478518546793793286704974382373670822285114741208958e27, tolerance);
        assertClose(fun, Math.scalb(11103367432951928.0, 62), 4.0411767712186990905102512019058204792570873633363159e36, tolerance);
        assertClose(fun, Math.scalb(11103367432951928.0, 326), 3.9754720509185529233002820161357111676582583112671658e116, tolerance);
        //
        // Super small values may cause spurious overflow:
        //
        double value = Double.MIN_NORMAL;
        while (value != 0) {
            Assertions.assertTrue(Double.isFinite(fun.applyAsDouble(value)));
            value /= 2;
        }

        // Simple check to ensure a zero length array is ignored
        final int[] signEmpty = {};
        for (final double z : new double[] {3.5, 6.76, 8.12}) {
            Assertions.assertEquals(BoostGamma.lgamma(z), BoostGamma.lgamma(z, signEmpty));
        }
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {"LGAMMA_.*"})
    void testLGamma(TestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @CsvSource({
        // Pole errors
        "-1, NaN",
        "-2, NaN",
        // Factorials: gamma(n+1)-1 = n! - 1
        "0, 0",
        "1, 0",
        "2, 1",
        "3, 5",
        "4, 23",
        "5, 119",
    })
    void testTGammap1m1EdgeCases(double z, double p) {
        Assertions.assertEquals(p, BoostGamma.tgamma1pm1(z));
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {"TGAMMAP1M1.*"})
    void testTGammap1m1(TestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @CsvSource({
        "0, 1, -1",
        "0, 0, 0",
        "0, -1, Infinity",
        "2, -2, -0.75",
        "2, 1024, Infinity",
        "2, -1075, -1",
        "NaN, 1, NaN",
        "1, NaN, NaN",
        // Negative x, even integer y
        "-2, 2, 3",
        // Negative x, non (even integer) y
        "-2, 2.1, NaN",
    })
    void testPowm1EdgeCases(double x, double y, double expected) {
        Assertions.assertEquals(expected, BoostMath.powm1(x, y));
    }

    @ParameterizedTest
    @EnumSource(value = BiTestCase.class, mode = Mode.MATCH_ANY, names = {"POWM1.*"})
    void testPowm1(BiTestCase tc) {
        assertFunction(tc);
    }

    /**
     * Test the log1pmx function with values that do not require high precision.
     *
     * @param x Argument x
     */
    @ParameterizedTest
    @ValueSource(doubles = {-1.1, -1, 0, 1, 1.5, 2, 3})
    void testLog1pmxStandard(double x) {
        Assertions.assertEquals(Math.log1p(x) - x, SpecialMath.log1pmx(x));
    }

    /**
     * Test the log1pmx function. The function is not a direct port of the Boost log1pmx
     * function so resides in the {@link SpecialMath} class. It is only used in {@link BoostGamma}
     * so tested here using the same test framework.
     */
    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {"LOG1PMX.*"})
    void testLog1pmx(TestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @CsvSource({
        // Argument a > 0
        "NaN, 1, NaN",
        "0, 1, NaN",
        "-1, 1, NaN",
        // Argument z >= 0
        "1, NaN, NaN",
        "1, -1, NaN",
    })
    void testIGammaEdgeCases(double a, double z, double expected) {
        // All functions have the same invalid domain for a and z
        Assertions.assertEquals(expected, BoostGamma.tgamma(a, z), "tgamma");
        Assertions.assertEquals(expected, BoostGamma.tgammaLower(a, z), "tgammaLower");
        Assertions.assertEquals(expected, BoostGamma.gammaP(a, z), "gammaP");
        Assertions.assertEquals(expected, BoostGamma.gammaQ(a, z), "gammaQ");
        Assertions.assertEquals(expected, BoostGamma.gammaPDerivative(a, z), "gammaPDerivative");
    }

    @ParameterizedTest
    @CsvSource({
        // z==0
        "2, 0, 0",
        "1, 0, 1",
        "0.5, 0, Infinity",
    })
    void testGammaPDerivativeEdgeCases(double a, double z, double expected) {
        Assertions.assertEquals(expected, BoostGamma.gammaPDerivative(a, z));
    }

    /**
     * tgamma spot tests extracted from
     * {@code boost/libs/math/test/test_igamma.hpp}.
     */
    @Test
    void testIGammaSpotTests() {
        int tolerance = 10;
        assertClose(BoostGamma::tgamma, 5, 1, 23.912163676143750903709045060494956383977723517065, tolerance);
        assertClose(BoostGamma::tgamma, 5, 5, 10.571838841565097874621959975919877646444998907920, tolerance);
        assertClose(BoostGamma::tgamma, 5, 10, 0.70206451384706574414638719662835463671916532623256, tolerance);
        assertClose(BoostGamma::tgamma, 5, 100, 3.8734332808745531496973774140085644548465762343719e-36, tolerance);
        assertClose(BoostGamma::tgamma, 0.5, 0.5, 0.56241823159440712427949495730204306902676756479651, tolerance);
        assertClose(BoostGamma::tgamma, 0.5, 9.0 / 10, 0.31853210360412109873859360390443790076576777747449, tolerance * 10);
        assertClose(BoostGamma::tgamma, 0.5, 5, 0.0027746032604128093194908357272603294120210079791437, tolerance);
        assertClose(BoostGamma::tgamma, 0.5, 100, 3.7017478604082789202535664481339075721362102520338e-45, tolerance);

        assertClose(BoostGamma::tgammaLower, 5, 1, 0.087836323856249096290954939505043616022276482935091, tolerance);
        assertClose(BoostGamma::tgammaLower, 5, 5, 13.428161158434902125378040024080122353555001092080, tolerance);
        assertClose(BoostGamma::tgammaLower, 5, 10, 23.297935486152934255853612803371645363280834673767, tolerance);
        assertClose(BoostGamma::tgammaLower, 5, 100, 23.999999999999999999999999999999999996126566719125, tolerance);

        assertClose(BoostGamma::gammaQ, 5, 1, 0.99634015317265628765454354418728984933240514654437, tolerance);
        assertClose(BoostGamma::gammaQ, 5, 5, 0.44049328506521241144258166566332823526854162116334, tolerance);
        assertClose(BoostGamma::gammaQ, 5, 10, 0.029252688076961072672766133192848109863298555259690, tolerance);
        assertClose(BoostGamma::gammaQ, 5, 100, 1.6139305336977304790405739225035685228527400976549e-37, tolerance);
        assertClose(BoostGamma::gammaQ, 1.5, 2, 0.26146412994911062220282207597592120190281060919079, tolerance);
        assertClose(BoostGamma::gammaQ, 20.5, 22, 0.34575332043467326814971590879658406632570278929072, tolerance);

        assertClose(BoostGamma::gammaP, 5, 1, 0.0036598468273437123454564558127101506675948534556288, tolerance);
        assertClose(BoostGamma::gammaP, 5, 5, 0.55950671493478758855741833433667176473145837883666, tolerance);
        assertClose(BoostGamma::gammaP, 5, 10, 0.97074731192303892732723386680715189013670144474031, tolerance);
        assertClose(BoostGamma::gammaP, 5, 100, 0.9999999999999999999999999999999999998386069466302, tolerance);
        assertClose(BoostGamma::gammaP, 1.5, 2, 0.73853587005088937779717792402407879809718939080921, tolerance);
        assertClose(BoostGamma::gammaP, 20.5, 22, 0.65424667956532673185028409120341593367429721070928, tolerance);

        // naive check on derivative function:
        tolerance = 50;
        assertClose(BoostGamma::gammaPDerivative, 20.5, 22,
            Math.exp(-22) * Math.pow(22, 19.5) / BoostGamma.tgamma(20.5), tolerance);

        // Bug reports from Rocco Romeo:
        assertClose(BoostGamma::tgamma, 20, Math.scalb(1.0, -40), 1.21645100408832000000e17, tolerance);
        assertClose(BoostGamma::tgammaLower, 20, Math.scalb(1.0, -40), 7.498484069471659696438206828760307317022658816757448882e-243, tolerance);
        assertClose(BoostGamma::gammaP, 20, Math.scalb(1.0, -40), 6.164230243774976473534975936127139110276824507876192062e-260, tolerance);

        assertClose(BoostGamma::tgamma, 30, Math.scalb(1.0, -30), 8.841761993739701954543616000000e30, tolerance);
        assertClose(BoostGamma::tgammaLower, 30, Math.scalb(1.0, -30), 3.943507283668378474979245322638092813837393749566146974e-273, tolerance);
        assertClose(BoostGamma::gammaP, 30, Math.scalb(1.0, -30), 4.460092102072560946444018923090222645613009128135650652e-304, tolerance);
        assertClose(BoostGamma::gammaPDerivative, 2, Math.scalb(1.0, -575), 8.08634922390438981326119906687585206568664784377654648227177e-174, tolerance);

        Assertions.assertEquals(BoostGamma.tgamma(176, 100), Double.POSITIVE_INFINITY);
        Assertions.assertEquals(BoostGamma.tgamma(530, 2000), Double.POSITIVE_INFINITY);
        Assertions.assertEquals(BoostGamma.tgamma(740, 2500), Double.POSITIVE_INFINITY);
        Assertions.assertEquals(BoostGamma.tgamma(530.5, 2000), Double.POSITIVE_INFINITY);
        Assertions.assertEquals(BoostGamma.tgamma(740.5, 2500), Double.POSITIVE_INFINITY);
        Assertions.assertEquals(BoostGamma.tgammaLower(10000.0f, 10000.0f / 4), Double.POSITIVE_INFINITY);
        assertClose(BoostGamma::tgamma, 170, 165, 2.737338337642022829223832094019477918166996032112404370e304, tolerance);
        assertClose(BoostGamma::tgammaLower, 170, 165, 1.531729671362682445715419794880088619901822603944331733e304, tolerance);
        // *** Increased from 10 * tolerance ***
        assertClose(BoostGamma::tgamma, 170, 170, 2.090991698081449410761040647015858316167077909285580375e304, 16 * tolerance);
        assertClose(BoostGamma::tgammaLower, 170, 170, 2.178076310923255864178211241883708221901740726771155728e304, 16 * tolerance);
        assertClose(BoostGamma::tgamma, 170, 190, 2.8359275512790301602903689596273175148895758522893941392e303, 10 * tolerance);
        assertClose(BoostGamma::tgammaLower, 170, 190, 3.985475253876802258910214992936834786579861050827796689e304, 10 * tolerance);
        // *** Increased from 10 * tolerance ***
        assertClose(BoostGamma::tgamma, 170, 1000, 6.1067635957780723069200425769800190368662985052038980542e72, 16 * tolerance);

        assertClose(BoostGamma::tgammaLower, 185, 1, 0.001999286058955490074702037576083582139834300307968257924836, tolerance);
        assertClose(BoostGamma::tgamma, 185, 1500, 1.037189524841404054867100938934493979112615962865368623e-67, tolerance * 10);

        assertClose(BoostGamma::tgamma, 36, Math.scalb(1.0, -26), 1.03331479663861449296666513375232000000e40, tolerance * 10);
        assertClose(BoostGamma::tgamma, 50.5, Math.scalb(1.0, -17), 4.2904629123519598109157551960589377e63, tolerance * 10);
        assertClose(BoostGamma::tgamma, 164.5, 0.125, 2.5649307433687542701168405519538910e292, tolerance * 10);
        //
        // Check very large parameters, see: https://github.com/boostorg/math/issues/168
        //
        final double maxVal = Double.MAX_VALUE;
        final double largeVal = maxVal * 0.99f;
        Assertions.assertEquals(BoostGamma.tgamma(22.25, maxVal), 0);
        Assertions.assertEquals(BoostGamma.tgamma(22.25, largeVal), 0);
        Assertions.assertEquals(BoostGamma.tgammaLower(22.25, maxVal), BoostGamma.tgamma(22.25));
        Assertions.assertEquals(BoostGamma.tgammaLower(22.25, largeVal), BoostGamma.tgamma(22.25));
        Assertions.assertEquals(BoostGamma.gammaQ(22.25, maxVal), 0);
        Assertions.assertEquals(BoostGamma.gammaQ(22.25, largeVal), 0);
        Assertions.assertEquals(BoostGamma.gammaP(22.25, maxVal), 1);
        Assertions.assertEquals(BoostGamma.gammaP(22.25, largeVal), 1);
        Assertions.assertEquals(BoostGamma.tgamma(22.25, Double.POSITIVE_INFINITY), 0);
        Assertions.assertEquals(BoostGamma.tgammaLower(22.25, Double.POSITIVE_INFINITY), BoostGamma.tgamma(22.25));
        Assertions.assertEquals(BoostGamma.gammaQ(22.25, Double.POSITIVE_INFINITY), 0);
        Assertions.assertEquals(BoostGamma.gammaP(22.25, Double.POSITIVE_INFINITY), 1);
        //
        // Large arguments and small parameters, see
        // https://github.com/boostorg/math/issues/451:
        //
        Assertions.assertEquals(BoostGamma.gammaQ(1770, 1e-12), 1);
        Assertions.assertEquals(BoostGamma.gammaP(1770, 1e-12), 0);
    }

    @ParameterizedTest
    @EnumSource(value = BiTestCase.class, mode = Mode.MATCH_ANY, names = {"IGAMMA_U.*"})
    void testIGammaUpper(BiTestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @EnumSource(value = BiTestCase.class, mode = Mode.MATCH_ANY, names = {"IGAMMA_L.*"})
    void testIGammaLower(BiTestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @EnumSource(value = BiTestCase.class, mode = Mode.MATCH_ANY, names = {"IGAMMA_Q.*"})
    void testIGammaQ(BiTestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @EnumSource(value = BiTestCase.class, mode = Mode.MATCH_ANY, names = {"IGAMMA_P.*"})
    void testIGammaP(BiTestCase tc) {
        assertFunction(tc);
    }

    /**
     * Test gamma Q with a value {@code a} so small that {@code tgamma(a) == inf}.
     */
    @ParameterizedTest
    @CsvSource({
        // This data needs verifying. Matlab, maxima and Boost all compute different values.
        // The values are the exponent of 2 to ensure tiny values are machine representable.
        // The following are from Boost using long double. This at least verifies that
        // the updated code computes close to the long double equivalent in the source
        // implementation.
        "-1074,-1074,3.67517082493672000135e-321",
        "-1074,-1073,3.67174622284245611609e-321",
        "-1074,-1072,3.66832162074819223109e-321",
        "-1074,-1000,3.4217502699611925034e-321",
        "-1074,-100,3.3960838512369590694e-322",
        "-1074,-10,3.13990203220802065461e-323",
        "-1074,-5,1.44243837956526236902e-323",
        "-1030,-1074,6.46542888972966038541e-308",
        "-1030,-1073,6.45940426601262169248e-308",
        "-1030,-1072,6.45337964229558300003e-308",
        "-1030,-1000,6.01960673466879712923e-308",
        "-1030,-100,5.97445389333973743599e-309",
        "-1030,-10,5.52377407118433787103e-310",
        "-1030,-5,2.53756443309180378013e-310",
    })
    void testGammaQTinyA(int ba, int bx, double q) {
        final double a = Math.scalb(1.0, ba);
        final double x = Math.scalb(1.0, bx);
        final double actualQ = BoostGamma.gammaQ(a, x);

        // Without changes to the Boost code that normalises by tgamma(a)
        // the result is zero. Check this does not occur.
        Assertions.assertNotEquals(0.0, actualQ);

        // Change tolerance for very small sub-normal result
        if (q < 1e-320) {
            // Within 1 ULP
            TestUtils.assertEquals(q, actualQ, 1);
        } else {
            // Sub-normal argument.
            // The worst case here is 260 ULP. The relative error is OK.
            final double relError = 1e-13;
            Assertions.assertEquals(q, actualQ, q * relError);
        }
    }

    /**
     * Gamma function with Lanczos support.
     *
     * <p>This is the original Boost implementation here for reference. For {@code z}
     * in the range [-20, 20] the function is not as accurate as the NSWC
     * Library of Mathematics Subroutines.
     *
     * @param z Argument z
     * @return gamma value
     */
    static double tgammaOriginal(double z) {
        double result = 1;

        if (z <= 0) {
            if (Math.rint(z) == z) {
                // Pole error
                return Double.NaN;
            }
            if (z <= -20) {
                result = BoostGamma.tgamma(-z) * BoostGamma.sinpx(z);
                // Checks for overflow, sub-normal or underflow have been removed.
                return -Math.PI / result;
            }

            // shift z to > 1:
            // Q. Is this comment old? The shift is to > 0.
            // Spot tests in the Boost resources test z -> 0 due to a bug report.
            while (z < 0) {
                result /= z;
                z += 1;
            }
        }
        //
        // z is > 0
        //

        // Updated condition from z < MAX_FACTORIAL
        if ((Math.rint(z) == z) && (z <= MAX_FACTORIAL + 1)) {
            // Gamma(n) = (n-1)!
            result *= FACTORIAL[(int) z - 1];
        } else if (z < ROOT_EPSILON) {
            result *= 1 / z - EULER;
        } else {
            result *= Lanczos.lanczosSum(z);
            final double zgh = z + Lanczos.g() - 0.5;
            final double lzgh = Math.log(zgh);
            if (z * lzgh > LOG_MAX_VALUE) {
                // we're going to overflow unless this is done with care:
                // Check for overflow has been removed:
                // if (lzgh * z / 2 > LOG_MAX_VALUE) ... overflow
                final double hp = Math.pow(zgh, (z / 2) - 0.25);
                result *= hp / Math.exp(zgh);
                // Check for overflow has been removed:
                // if (Double.MAX_VALUE / hp < result) ... overflow
                result *= hp;
            } else {
                result *= Math.pow(zgh, z - 0.5) / Math.exp(zgh);
            }
        }
        return result;
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
     * Assert the function is close to the expected value.
     *
     * @param fun Function
     * @param x Input value
     * @param y Input value
     * @param expected Expected value
     * @param tolerance the tolerance
     */
    private static void assertClose(DoubleBinaryOperator fun, double x, double y, double expected, int tolerance) {
        final double actual = fun.applyAsDouble(x, y);
        TestUtils.assertEquals(expected, actual, tolerance, null, () -> x + ", " + y);
    }

    /**
     * Assert the function using extended precision.
     *
     * @param tc Test case
     */
    private static void assertFunction(TestCase tc) {
        final TestUtils.RMS rms = new TestUtils.RMS();
        try (DataReader in = new DataReader(tc.getFilename())) {
            for (String[] tokens = in.next(); tokens != null; tokens = in.next()) {
                try {
                    final double x = Double.parseDouble(tokens[0]);
                    final BigDecimal expected = new BigDecimal(tokens[tc.getExpectedField()]);
                    final double actual = tc.getFunction().applyAsDouble(x);
                    TestUtils.assertEquals(expected, actual, tc.getTolerance(), rms::add,
                        () -> tc + " x=" + x);
                } catch (final NumberFormatException ex) {
                    Assertions.fail("Failed to load data: " + Arrays.toString(tokens), ex);
                }
            }
        } catch (final IOException ex) {
            Assertions.fail("Failed to load data: " + tc.getFilename(), ex);
        }

        assertRms(tc, rms);
    }

    /**
     * Assert the function using extended precision.
     *
     * @param tc Test case
     */
    private static void assertFunction(BiTestCase tc) {
        final TestUtils.RMS rms = new TestUtils.RMS();
        try (DataReader in = new DataReader(tc.getFilename())) {
            for (String[] tokens = in.next(); tokens != null; tokens = in.next()) {
                try {
                    final double x = Double.parseDouble(tokens[0]);
                    final double y = Double.parseDouble(tokens[1]);
                    final BigDecimal expected = new BigDecimal(tokens[tc.getExpectedField()]);
                    final double actual = tc.getFunction().applyAsDouble(x, y);
                    TestUtils.assertEquals(expected, actual, tc.getTolerance(), rms::add,
                        () -> tc + " x=" + x + ", y=" + y);
                } catch (final NumberFormatException ex) {
                    Assertions.fail("Failed to load data: " + Arrays.toString(tokens), ex);
                }
            }
        } catch (final IOException ex) {
            Assertions.fail("Failed to load data: " + tc.getFilename(), ex);
        }

        assertRms(tc, rms);
    }

    /**
     * Class to read data fields from a test resource file.
     */
    private static class DataReader implements AutoCloseable {
        /** Pattern to split data fields. */
        private static final Pattern FIELD_PATTERN = Pattern.compile("[, ]+");

        /** Input to read. */
        private final BufferedReader in;

        /**
         * @param filename Resource filename to read
         */
        DataReader(String filename) {
            final InputStream resourceAsStream = this.getClass().getResourceAsStream(filename);
            Assertions.assertNotNull(resourceAsStream, () -> "Could not find resource " + filename);
            in = new BufferedReader(new InputStreamReader(resourceAsStream));
        }

        /**
         * Get the next line of data (or null).
         *
         * @return data
         * @throws IOException Signals that an I/O exception has occurred.
         */
        String[] next() throws IOException {
            for (String line = in.readLine(); line != null; line = in.readLine()) {
                if (line.startsWith("#") || line.trim().isEmpty()) {
                    continue;
                }
                return FIELD_PATTERN.split(line);
            }
            return null;
        }

        @Override
        public void close() throws IOException {
            in.close();
        }
    }

    /**
     * Assert the Root Mean Square (RMS) error of the function is below the allowed maximum.
     *
     * @param te Test error
     */
    private static void assertRms(TestError te, TestUtils.RMS data) {
        final double rms = data.getRMS();
        debugRms(te.toString(), data.getMax(), rms);
        Assertions.assertTrue(rms <= te.getRmsTolerance(),
            () -> String.format("%s RMS %s < %s", te, rms, te.getRmsTolerance()));
    }

    /**
     * Output the maximum and RMS ulp for the named test.
     * Used for reporting the errors and setting appropriate test tolerances.
     * This is relevant across different JDK implementations where the java.util.Math
     * functions used in BoostErf may compute to different accuracy.
     *
     * @param name Test name
     * @param maxUlp Maximum ulp
     * @param rmsUlp RMS ulp
     */
    private static void debugRms(String name, double maxUlp, double rmsUlp) {
        // CHECKSTYLE: stop regexp
        // Debug output of max and RMS error.
        System.out.printf("%-25s   max %10.6g   RMS %10.6g%n", name, maxUlp, rmsUlp);
    }
}
