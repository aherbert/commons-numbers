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
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Stream;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.MethodOrderer;
import org.junit.jupiter.api.Order;
import org.junit.jupiter.api.TestMethodOrder;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.EnumSource.Mode;

/**
 * Tests for {@link BoostDigamma} and {@link BoostTrigamma}.
 *
 * <p>Note: Resource data files used in these tests have been extracted
 * from the Boost test files for the digamma and trigamma functions.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class BoostDiTrigammaTest {
    /**
     * Reporting flag. When {@code >= 0} the version 1.3 Di/Trigamma functions
     * are tested and the error is reported to stdout.
     */
    private static int reporting = -1;

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
     */
    private enum TestCase implements TestError {
        // Note: The original Di/Trigamma fixtures are not run by default.
        // They can be enabled for reporting.
        T_DIGAMMA_DATA(BoostDiTrigammaTest.Digamma13::value, "digamma_data.csv", 9, 0.65),
        T_DIGAMMA_INT_DATA(BoostDiTrigammaTest.Digamma13::value, "digamma_int_data.csv", 14, 1.9),
        T_DIGAMMA_HALF_INT_DATA(BoostDiTrigammaTest.Digamma13::value, "digamma_half_int_data.csv", 125, 12.5),
        T_DIGAMMA_ROOT_DATA(BoostDiTrigammaTest.Digamma13::value, "digamma_root_data.csv", 16000, 1600),
        T_DIGAMMA_SMALL_DATA(BoostDiTrigammaTest.Digamma13::value, "digamma_small_data.csv", 1.1e6, 200000),
        T_DIGAMMA_NEG_DATA(BoostDiTrigammaTest.Digamma13::value, "digamma_neg_data.csv", 54000, 4600),
        DIGAMMA_DATA(BoostDigamma::digamma, "digamma_data.csv", 1.1, 0.4),
        DIGAMMA_INT_DATA(BoostDigamma::digamma, "digamma_int_data.csv", 1.0, 0.4),
        DIGAMMA_HALF_INT_DATA(BoostDigamma::digamma, "digamma_half_int_data.csv", 1.0, 0.35),
        DIGAMMA_ROOT_DATA(BoostDigamma::digamma, "digamma_root_data.csv", 1.5, 0.6),
        DIGAMMA_SMALL_DATA(BoostDigamma::digamma, "digamma_small_data.csv", 0.9, 0.4),
        DIGAMMA_NEG_DATA(BoostDigamma::digamma, "digamma_neg_data.csv", 340, 30),
        T_TRIGAMMA_INT_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_int_data.csv", -22000, 4700),
        T_TRIGAMMA_HALF_INT_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_half_int_data.csv", 20000, 11000),
        T_TRIGAMMA_SMALL_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_small_data.csv", 450000, 46000),
        T_TRIGAMMA_MED_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_medium_data.csv", 6000, 3000),
        T_TRIGAMMA_NEG_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_neg_data.csv", 14000, 7000),
        T_TRIGAMMA_NEG_SMALL_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_neg_small_data.csv", 22, 2.5),
        T_TRIGAMMA_NEG_HALF_DATA(BoostDiTrigammaTest.Trigamma13::value, "trigamma_neg_half_int_data.csv", 37, 37),
        TRIGAMMA_INT_DATA(BoostTrigamma::trigamma, "trigamma_int_data.csv", 1.45, 0.5),
        TRIGAMMA_HALF_INT_DATA(BoostTrigamma::trigamma, "trigamma_half_int_data.csv", 1.25, 0.5),
        TRIGAMMA_SMALL_DATA(BoostTrigamma::trigamma, "trigamma_small_data.csv", 0, 0),
        TRIGAMMA_MED_DATA(BoostTrigamma::trigamma, "trigamma_medium_data.csv", 1.5, 0.5),
        TRIGAMMA_NEG_DATA(BoostTrigamma::trigamma, "trigamma_neg_data.csv", 1.0, 0.45),
        TRIGAMMA_NEG_SMALL_DATA(BoostTrigamma::trigamma, "trigamma_neg_small_data.csv", 1.65, 0.40),
        TRIGAMMA_NEG_HALF_DATA(BoostTrigamma::trigamma, "trigamma_neg_half_int_data.csv", 0.9, 0.5);

//        Example output with reporting = 0
//        T_DIGAMMA_DATA                        max      7.93865   RMS   0.601620   mean      0.0349363  n  500
//        T_DIGAMMA_INT_DATA                    max      12.9555   RMS    1.71244   mean       0.334989  n  100
//        T_DIGAMMA_HALF_INT_DATA               max      117.972   RMS    11.8459   mean        1.82825  n  101
//        T_DIGAMMA_ROOT_DATA                   max      14861.0   RMS    1498.85   mean       -100.196  n  200
//        T_DIGAMMA_SMALL_DATA                  max  1.02488e+06   RMS     185338   mean       -42571.8  n   33
//        T_DIGAMMA_NEG_DATA                    max      52911.0   RMS    4503.61   mean        575.677  n  200
//        DIGAMMA_DATA                          max     0.987831   RMS   0.316799   mean     0.00355768  n  500
//        DIGAMMA_INT_DATA                      max     0.916689   RMS   0.358642   mean     -0.0749518  n  100
//        DIGAMMA_HALF_INT_DATA                 max     0.927770   RMS   0.289511   mean     -0.0282425  n  101
//        DIGAMMA_ROOT_DATA                     max      1.44295   RMS   0.539360   mean      -0.314893  n  200
//        DIGAMMA_SMALL_DATA                    max     0.851235   RMS   0.303128   mean     -0.0242290  n   33
//        DIGAMMA_NEG_DATA                      max      80.9962   RMS    8.17010   mean       -1.41228  n  200
//        T_TRIGAMMA_INT_DATA                   max      20231.9   RMS    4540.32   mean       -1496.00  n  500
//        T_TRIGAMMA_HALF_INT_DATA              max      18843.4   RMS    9616.21   mean       -7081.83  n  101
//        T_TRIGAMMA_SMALL_DATA                 max       431205   RMS    44759.0   mean       -5807.85  n   99
//        T_TRIGAMMA_MED_DATA                   max      5059.44   RMS    2793.63   mean       -2203.45  n  200
//        T_TRIGAMMA_NEG_DATA                   max      13354.0   RMS    6035.51   mean        411.373  n   41
//        T_TRIGAMMA_NEG_SMALL_DATA             max      20.0821   RMS    2.07840   mean       0.239588  n   99
//        T_TRIGAMMA_NEG_HALF_DATA              max      36.8791   RMS    36.4224   mean        36.4213  n   20
//        TRIGAMMA_INT_DATA                     max      1.37569   RMS   0.447972   mean     -0.0137117  n  500
//        TRIGAMMA_HALF_INT_DATA                max      1.17836   RMS   0.457438   mean      0.0780938  n  101
//        TRIGAMMA_SMALL_DATA                   max      0.00000   RMS    0.00000   mean        0.00000  n   99
//        TRIGAMMA_MED_DATA                     max      1.41307   RMS   0.460636   mean     -0.0158794  n  200
//        TRIGAMMA_NEG_DATA                     max     0.900578   RMS   0.411192   mean       0.228960  n   41
//        TRIGAMMA_NEG_SMALL_DATA               max      1.58877   RMS   0.344048   mean    -0.00865301  n   99
//        TRIGAMMA_NEG_HALF_DATA                max     0.841433   RMS   0.469806   mean      -0.296449  n   20

        /** The function. */
        private final DoubleUnaryOperator fun;

        /** The filename containing the test data. */
        private final String filename;

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
            this.fun = fun;
            this.filename = filename;
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
     * Copy of Commons Numbers Digamma from version 1.3.
     * This includes fixes applied before the 1.4 release to the validated domain (NUMBERS-219).
     */
    private static final class Digamma13 {
        /** <a href="https://en.wikipedia.org/wiki/Euler-Mascheroni_constant">Euler-Mascheroni constant</a>. */
        private static final double GAMMA = 0.577215664901532860606512090082;

        /** C limit. */
        private static final double C_LIMIT = 49;
        /** S limit. */
        private static final double S_LIMIT = 1e-5;
        /** Fraction. */
        private static final double F_M1_12 = -1d / 12;
        /** Fraction. */
        private static final double F_1_120 = 1d / 120;
        /** Fraction. */
        private static final double F_M1_252 = -1d / 252;

        /** Private constructor. */
        private Digamma13() {
            // intentionally empty.
        }

        /**
         * Computes the digamma function.
         *
         * <p>This is an independently written implementation of the algorithm described in
         * <a href="https://www.uv.es/~bernardo/1976AppStatist.pdf">Jose Bernardo,
         * Algorithm AS 103: Psi (Digamma) Function, Applied Statistics, 1976</a>.
         * A <a href="https://en.wikipedia.org/wiki/Digamma_function#Reflection_formula">
         * reflection formula</a> is incorporated to improve performance on negative values.
         *
         * <p>Some of the constants have been changed to increase accuracy at the moderate
         * expense of run-time.  The result should be accurate to within {@code 1e-8}.
         * relative tolerance for {@code 0 < x < 1e-5}  and within {@code 1e-8} absolute
         * tolerance otherwise.
         *
         * @param x Argument.
         * @return digamma(x) to within {@code 1e-8} relative or absolute error whichever
         * is larger.
         */
        static double value(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // large negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? x : Double.NaN;
            }

            double digamma = 0;
            if (x < 0) {
                // Use reflection formula to fall back into positive values.
                digamma = -Math.PI / Math.tan(Math.PI * x);
                x = 1 - x;
            }

            if (x <= S_LIMIT) {
                // Use method 5 from Bernardo AS103, accurate to O(x).
                return digamma - GAMMA - 1 / x;
            }

            while (x < C_LIMIT) {
                digamma -= 1 / x;
                x += 1;
            }

            // Use method 4, accurate to O(1/x^8)
            final double inv = 1 / (x * x);
            //            1       1        1         1
            // log(x) -  --- - ------ + ------- - -------
            //           2 x   12 x^2   120 x^4   252 x^6
            digamma += Math.log(x) - 0.5 / x + inv * (F_M1_12 + inv * (F_1_120 + F_M1_252 * inv));

            return digamma;
        }
    }

    /**
     * Copy of Commons Numbers Trigamma from version 1.3.
     * This includes fixes applied before the 1.4 release to the validated domain (NUMBERS-219);
     * and to use the reflection formula (NUMBERS-209).
     */
    private static final class Trigamma13 {
        /** C limit. */
        private static final double C_LIMIT = 49;

        /** S limit. */
        private static final double S_LIMIT = 1e-5;
        /** Fraction. */
        private static final double F_1_6 = 1d / 6;
        /** Fraction. */
        private static final double F_1_30 = 1d / 30;
        /** Fraction. */
        private static final double F_1_42 = 1d / 42;

        /** Private constructor. */
        private Trigamma13() {
            // intentionally empty.
        }

        /**
         * Computes the trigamma function.
         *
         * <p>Special cases:
         * <ul>
         * <li>If the argument is infinity, then the result is 0.</li>
         * <li>If the argument is negative infinity, then the result is nan.</li>
         * <li>If the argument is nan, then the result is nan.</li>
         * </ul>
         *
         * @param x Argument.
         * @return trigamma(x) to within {@code 1e-8} relative or absolute error whichever is larger.
         */
        static double value(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // large negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? 0 : Double.NaN;
            }

            if (x <= 0) {
                // Use the reflection formula:
                // trigamma(x) + trigamma(1 - x) = pi^2 / sin^2(pi * x)
                // to fall back into positive values. Without this the increment
                // loop below does not terminate for large negative x: when
                // ulp(x) >= 2 the update x += 1 leaves x unchanged (or advances
                // a single step and then sticks) and x < C_LIMIT holds forever.

                // negative integers are poles
                if (Math.rint(x) == x) {
                    return Double.POSITIVE_INFINITY;
                }

                final double s = Math.sin(Math.PI * x);
                return Math.PI * Math.PI / (s * s) - value(1 - x);
            }

            if (x <= S_LIMIT) {
                return 1 / (x * x);
            }

            double trigamma = 0;
            while (x < C_LIMIT) {
                trigamma += 1 / (x * x);
                x += 1;
            }

            final double inv = 1 / (x * x);
            //  1     1       1       1        1
            //  - + ----- + ----- - ------ + ------
            //  x   2 x^2   6 x^3   30 x^5   42 x^7
            trigamma += 1 / x + 0.5 * inv + inv / x * (F_1_6 - inv * (F_1_30 + F_1_42 * inv));

            return trigamma;
        }
    }

    @ParameterizedTest
    @CsvSource({
        "Infinity, Infinity",
        "NaN, NaN",
        "-Infinity, NaN",
        // Pole errors
        "0, NaN",
        "-1, NaN",
        "-2, NaN",
    })
    void testDigammaEdgeCases(double x, double expected) {
        Assertions.assertEquals(expected, BoostDigamma.digamma(x));
    }

    @ParameterizedTest
    @MethodSource
    void testDigammaSpot(double x, double z, int ulp) {
        assertClose(BoostDigamma::digamma, x, z, ulp);
    }

    static Stream<Arguments> testDigammaSpot() {
        // Boost uses tolerance = 3
        // 5 fail at tolerance = 0
        // 1 fail at tolerance = 1
        final int tolerance = 2;
        return Stream.of(
            // Result computed using mpmath.digamma (version 1.4.1)
            Arguments.of(Double.MAX_VALUE, 709.78271289338399673222338991066, 0),
            // digamma spot tests extracted from
            // boost/libs/math/test/test_digamma.cpp
            Arguments.of(0.125, -8.3884926632958548678027429230863430000514460424495, tolerance),
            Arguments.of(0.5, -1.9635100260214234794409763329987555671931596046604, tolerance),
            Arguments.of(1, -0.57721566490153286060651209008240243104215933593992, tolerance),
            Arguments.of(1.5, 0.036489973978576520559023667001244432806840395339566, tolerance),
            Arguments.of(1.5 - 1.0 / 32, 0.00686541147073577672813890866512415766586241385896200579891429, tolerance),
            Arguments.of(2, 0.42278433509846713939348790991759756895784066406008, tolerance),
            Arguments.of(8, 2.0156414779556099965363450527747404261006978069172, tolerance),
            Arguments.of(12, 2.4426616799758120167383652547949424463027180089374, tolerance),
            Arguments.of(22, 3.0681430398611966699248760264450329818421699570581, tolerance),
            Arguments.of(50, 3.9019896734278921969539597028823666609284424880275, tolerance),
            Arguments.of(500, 6.2136077650889917423827750552855712637776544784569, tolerance),
            // negative values:
            Arguments.of(-0.125, 7.1959829284523046176757814502538535827603450463013, tolerance),
            Arguments.of(-10.125, 9.9480538258660761287008034071425343357982429855241, tolerance),
            Arguments.of(-10.875, -5.1527360383841562620205965901515879492020193154231, tolerance),
            Arguments.of(-1.5, 0.70315664064524318722569033366791109947350706200623, tolerance),
            // Test cases from Rocco Romeo: Values near 0
            Arguments.of(Math.scalb(1.0, -100), -1.26765060022822940149670320537657721566490153286060651209008e30, tolerance),
            Arguments.of(-Math.scalb(1.0, -100), 1.26765060022822940149670320537542278433509846713939348790992e30, tolerance),
            // Values near 1
            Arguments.of(1.0, -0.577215664901532860606512090082402431042159335939923598805767, tolerance),
            Arguments.of(-1 + Math.scalb(1.0, -20), -1.04857557721314249602848739817764518743062133735858753112190e6, tolerance),
            Arguments.of(-1 - Math.scalb(1.0, -20), 1.04857642278181269259522681939281063878220298942888100442172e6, tolerance)
        );
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {".*DIGAMMA_.*"})
    @Order(1)
    void testDigamma(TestCase tc) {
        assertFunction(tc);
    }

    @ParameterizedTest
    @CsvSource({
        "Infinity, 0",
        "NaN, NaN",
        "-Infinity, NaN",
        // Pole errors
        "0, Infinity",
        "-1, Infinity",
        "-2, Infinity",
    })
    void testTrigammaEdgeCases(double x, double expected) {
        Assertions.assertEquals(expected, BoostTrigamma.trigamma(x));
    }

    @ParameterizedTest
    @MethodSource
    void testTrigammaSpot(double x, double z, int ulp) {
        assertClose(BoostTrigamma::trigamma, x, z, ulp);
    }

    static Stream<Arguments> testTrigammaSpot() {
        // Boost has no spot tests.
        // The Boost data covers integer x up to 500.
        // Supplement with some large real x.

        // Created using mpmath (1.14.1) polygamma
        // from mpmath import polygamma, mp
        // mp.pretty = True
        // def tri(x):
        //   print(f'Arguments.of({x}, {polygamma(1, x)}, tolerance),')
        final int tolerance = 1;
        return Stream.of(
            Arguments.of(123.45, 0.00813334272030528333111760362484, tolerance),
            Arguments.of(222.33, 0.00450794890925195733799541918749, tolerance),
            Arguments.of(322.33, 0.00310722802544591037240894139603, tolerance),
            Arguments.of(422.33, 0.00237062212663800648031668822254, tolerance),
            Arguments.of(722.33, 0.00138536752429675238049289841163, tolerance),
            Arguments.of(987.65, 0.00101301718531432031617733361854, tolerance),
            Arguments.of(10101, 9.90049996716169996339917355728e-05, tolerance),
            Arguments.of(125377, 7.97597635916378787167637896026e-06, tolerance),
            Arguments.of(356757, 2.8030323203623346759416706408e-06, tolerance),
            // Large x is very accurate (asymptotes ~ 1/x)
            Arguments.of(6786788.563, 1.47345105036028318762770453093e-07, 0),
            Arguments.of(12703812797.123, 7.87165251888650129685144629529e-11, 0),
            Arguments.of(1.23e+50, 8.13008130081300847073926230118e-51, 0),
            Arguments.of(1.23e+100, 8.13008130081300757704653639729e-101, 0),
            Arguments.of(1.23e+200, 8.13008130081300808406093381094e-201, 0),
            Arguments.of(1.23e+300, 8.13008130081300848952672063138e-301, 0),
            Arguments.of(1.23e+305, 8.13008130081300875278444339342e-306, 0),
            Arguments.of(1.23e+307, 8.13008130081300818593426624807e-308, 0),
            Arguments.of(1.23e+308, 8.13008130081300851573800567809e-309, 0),
            // No sub-normals can approach MIN_VALUE (4.9e-324)
            Arguments.of(Double.MAX_VALUE, 5.56268464626800407530763909489e-309, 0)
        );
    }

    @ParameterizedTest
    @EnumSource(value = TestCase.class, mode = Mode.MATCH_ANY, names = {".*TRIGAMMA_.*"})
    @Order(2)
    void testTrigamma(TestCase tc) {
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
        // Skip testing old functions
        if (reporting < 0 && tc.name().startsWith("T_")) {
            return;
        }
        final TestUtils.ErrorStatistics stats = new TestUtils.ErrorStatistics();
        try (DataReader in = new DataReader(tc.getFilename())) {
            while (in.next()) {
                try {
                    final double x = in.getDouble(0);
                    final BigDecimal expected = in.getBigDecimal(1);
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
        if (reporting < 0) {
            return;
        }
        // CHECKSTYLE: stop regexp
        System.out.printf("%-35s   max %12.6g   RMS %10.6g   mean %14.6g  n %4d%n",
            name, maxAbsUlp, rmsUlp, meanUlp, size);
    }
}
