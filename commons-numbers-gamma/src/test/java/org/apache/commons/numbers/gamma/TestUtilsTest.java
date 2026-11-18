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

import java.math.BigDecimal;
import org.apache.commons.numbers.gamma.TestUtils.ErrorStatistics;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Tests for {@link TestUtils}.
 * This verifies the custom ULP assertions function as expected.
 */
class TestUtilsTest {
    @Test
    void testAssertEquals() {
        final double x = 1.23;
        final double ulp = Math.ulp(x);
        Assertions.assertEquals(0, TestUtils.assertEquals(x, x, 0));
        Assertions.assertEquals(0, TestUtils.assertEquals(x, x, 1));
        Assertions.assertEquals(0, TestUtils.assertEquals(x, x, 30));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(x, x + ulp, 0));
        Assertions.assertEquals(1, TestUtils.assertEquals(x, x + ulp, 1));
        Assertions.assertEquals(-1, TestUtils.assertEquals(x, x - ulp, 1));
        Assertions.assertEquals(1, TestUtils.assertEquals(x, x + ulp, 30));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(x, x + 30 * ulp, 29));
        Assertions.assertEquals(30, TestUtils.assertEquals(x, x + 30 * ulp, 30));
        Assertions.assertEquals(-30, TestUtils.assertEquals(x, x - 30 * ulp, 30));

        // Check order and sign
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(-x - ulp, x, 0));
        Assertions.assertEquals(1, TestUtils.assertEquals(-x - ulp, -x, 1));
        Assertions.assertEquals(1, TestUtils.assertEquals(x - ulp, x, 1));
        Assertions.assertEquals(2, TestUtils.assertEquals(x - ulp, x + ulp, 2));

        // Opposite signs
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(x, -x, 500000));
    }

    @Test
    void testAssertEqualsZero() {
        // These are not binary equal but for numeric purposes they are treated as equal
        final double[] zero = {-0.0, 0.0};
        for (final double a : zero) {
            for (final double b : zero) {
                Assertions.assertEquals(0, TestUtils.assertEquals(a, b, 0));
                Assertions.assertEquals(0, TestUtils.assertEquals(a, b, 1));
                Assertions.assertEquals(0, TestUtils.assertEquals(a, b, 30));
            }
        }

        // Difference from zero
        final double x = Double.MIN_VALUE;
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(0.0, x, 0));
        Assertions.assertEquals(1, TestUtils.assertEquals(0.0, x, 1));
        Assertions.assertEquals(1, TestUtils.assertEquals(-0.0, x, 1));
        Assertions.assertEquals(-1, TestUtils.assertEquals(0.0, -x, 1));
        Assertions.assertEquals(-1, TestUtils.assertEquals(-0.0, -x, 1));
        Assertions.assertEquals(2, TestUtils.assertEquals(0.0, 2 * x, 2));
        Assertions.assertEquals(2, TestUtils.assertEquals(-0.0, 2 * x, 2));
    }

    @Test
    void testAssertEqualsEdgeCases() {
        final double nan = Double.NaN;
        final double inf = Double.POSITIVE_INFINITY;
        final double max = Double.MAX_VALUE;
        Assertions.assertEquals(0, TestUtils.assertEquals(nan, nan, 0));
        Assertions.assertEquals(0, TestUtils.assertEquals(inf, inf, 0));
        Assertions.assertEquals(0, TestUtils.assertEquals(max, max, 0));

        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(nan, inf, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(nan, max, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(inf, max, 0));

        Assertions.assertEquals(-1, TestUtils.assertEquals(inf, max, 1));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(inf, -inf, 0));
    }

    @Test
    void testBigDecimalAssertEquals() {
        final double x = 1.23;
        final double ulp = Math.ulp(x);

        final BigDecimal bdx = new BigDecimal(x);

        Assertions.assertEquals(0, TestUtils.assertEquals(bdx, x, 0));
        Assertions.assertEquals(0, TestUtils.assertEquals(bdx, x, 1));
        Assertions.assertEquals(0, TestUtils.assertEquals(bdx, x, 30));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdx, x + ulp, 0));
        Assertions.assertEquals(1, TestUtils.assertEquals(bdx, x + ulp, 1));
        Assertions.assertEquals(-1, TestUtils.assertEquals(bdx, x - ulp, 1));
        Assertions.assertEquals(1, TestUtils.assertEquals(bdx, x + ulp, 30));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdx, x + 30 * ulp, 29));
        Assertions.assertEquals(30, TestUtils.assertEquals(bdx, x + 30 * ulp, 30));
        Assertions.assertEquals(-30, TestUtils.assertEquals(bdx, x - 30 * ulp, 30));

        // Opposite signs
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdx, -x, 500000));

        // High precision: Add 1/4 ulp to the value:
        //    ---|----------|----------|-----
        //       -1         x  |       +1
        //                     y
        final BigDecimal bdy = bdx.add(new BigDecimal(ulp).divide(BigDecimal.valueOf(4)));
        Assertions.assertEquals(0, TestUtils.assertEquals(bdy, x, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdy, x + ulp, 0));
        Assertions.assertEquals(0.75, TestUtils.assertEquals(bdy, x + ulp, 1));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdy, x - ulp, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdy, x - ulp, 1));
        Assertions.assertEquals(-1.25, TestUtils.assertEquals(bdy, x - ulp, 1.5), 1e-16);
    }

    @Test
    void testBigDecimalAssertEqualsZero() {
        // These are not binary equal but for numeric purposes they are treated as equal
        final double[] zero = {-0.0, 0.0};
        final BigDecimal a = BigDecimal.ZERO;
        for (final double b : zero) {
            Assertions.assertEquals(0, TestUtils.assertEquals(a, b, 0));
            Assertions.assertEquals(0, TestUtils.assertEquals(a, b, 1));
            Assertions.assertEquals(0, TestUtils.assertEquals(a, b, 30));
        }

        // Difference from zero
        final double x = Double.MIN_VALUE;
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(a, x, 0));
        Assertions.assertEquals(1, TestUtils.assertEquals(a, x, 1));
        Assertions.assertEquals(-1, TestUtils.assertEquals(a, -x, 1));
        Assertions.assertEquals(2, TestUtils.assertEquals(a, 2 * x, 2));
        Assertions.assertEquals(-2, TestUtils.assertEquals(a, -2 * x, 2));
    }

    @Test
    void testBigDecimalAssertEqualsEdgeCases() {
        final double nan = Double.NaN;
        final double inf = Double.POSITIVE_INFINITY;
        final double max = Double.MAX_VALUE;

        final BigDecimal bdinf = new BigDecimal("1e310");
        final BigDecimal bdmax = new BigDecimal(max);

        Assertions.assertEquals(0, TestUtils.assertEquals(bdinf, inf, 0));
        Assertions.assertEquals(0, TestUtils.assertEquals(bdmax, max, 0));

        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdinf, max, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdinf, nan, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdmax, inf, 0));
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdmax, nan, 0));

        final double x = Math.nextDown(max);
        Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(bdmax, x, 0));
        Assertions.assertEquals(-1, TestUtils.assertEquals(bdmax, x, 1));
    }

    @Test
    void testBigDecimalAssertEqualsOverflowMargin() {
        final double inf = Double.POSITIVE_INFINITY;
        final double max = Double.MAX_VALUE;
        final BigDecimal ulp = new BigDecimal(Math.ulp(max));
        final BigDecimal bdmax = new BigDecimal(max);
        // 2^1024: treated as the value of infinity
        final BigDecimal two1024 = bdmax.add(ulp);
        final BigDecimal half = new BigDecimal("0.5");
        final BigDecimal quarter = new BigDecimal("0.25");

        for (final int sign : new int[] {1, -1}) {
            final BigDecimal sbdmax = sign > 0 ? bdmax : bdmax.negate();
            final BigDecimal s2 = sign > 0 ? two1024 : two1024.negate();
            final BigDecimal sulp = sign > 0 ? ulp : ulp.negate();
            final double sinf = sign * inf;
            final double smax = sign * max;

            // Expected rounds to MAX_VALUE; actual is infinity
            Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(sbdmax, sinf, 0.5));
            Assertions.assertEquals(sign, TestUtils.assertEquals(sbdmax, sinf, 1));
            // Expected 1/4 ulp above MAX_VALUE (rounds to MAX_VALUE)
            final BigDecimal e1 = sbdmax.add(sulp.multiply(quarter));
            Assertions.assertEquals(sign * 0.75, TestUtils.assertEquals(e1, sinf, 1));

            // Expected rounds to infinity: at the rounding threshold MAX_VALUE + ulp/2
            final BigDecimal e2 = sbdmax.add(sulp.multiply(half));
            Assertions.assertEquals(sinf, e2.doubleValue());
            Assertions.assertEquals(0, TestUtils.assertEquals(e2, sinf, 0));
            Assertions.assertEquals(-sign * 0.5, TestUtils.assertEquals(e2, smax, 0.5));
            // Expected 2^1024 (1 + 2^-53): 1.5 ulp above MAX_VALUE
            final BigDecimal e3 = s2.add(sulp.multiply(half));
            Assertions.assertThrows(AssertionError.class, () -> TestUtils.assertEquals(e3, smax, 1));
            Assertions.assertEquals(-sign * 1.5, TestUtils.assertEquals(e3, smax, 1.5));
            final double below = sign > 0 ? Math.nextDown(max) : Math.nextUp(-max);
            Assertions.assertEquals(-sign * 2.5, TestUtils.assertEquals(e3, below, 2.5));

            // Just below 2^1025: large finite error (2^1025 - ulp/2 - MAX_VALUE = (2^53 + 0.5) ulp)
            final BigDecimal s2x2 = s2.add(s2);
            final BigDecimal below2 = s2x2.subtract(sulp.multiply(half));
            Assertions.assertEquals(-sign * (0x1.0p53 + 0.5), failedError(below2, smax));
            // At or above 2^1025 (definite overflow): infinite error
            Assertions.assertEquals(-sign * inf, failedError(s2x2, smax));
            Assertions.assertEquals(-sign * inf, failedError(new BigDecimal(sign > 0 ? "1e310" : "-1e310"), smax));

            // Finite value below the top binade: no margin, infinite error
            final double small = sign * 0x1.0p1022;
            Assertions.assertEquals(sign * inf, failedError(new BigDecimal(small), sinf));
            Assertions.assertEquals(-sign * inf, failedError(e2, small));
        }

        // Opposite signs: no margin
        Assertions.assertEquals(-inf, failedError(bdmax, -inf));
        Assertions.assertEquals(inf, failedError(bdmax.negate(), inf));
    }

    /**
     * Assert the values are not equal within 1 ulp and return the error.
     *
     * @param expected Expected value
     * @param actual Actual value
     * @return the error
     */
    private static double failedError(BigDecimal expected, double actual) {
        final double[] error = {Double.NaN};
        Assertions.assertThrows(AssertionError.class,
            () -> TestUtils.assertEquals(expected, actual, 1, e -> error[0] = e, null));
        return error[0];
    }

    @Test
    void testErrorStatistics() {
        final ErrorStatistics stats = new ErrorStatistics();
        assertErrorStatistics(stats, 0, Double.NaN, Double.NaN);
        stats.add(0);
        assertErrorStatistics(stats, 0, 0, 0);
        stats.add(0.5);
        assertErrorStatistics(stats, 0.5, Math.sqrt(0.25 / 2), 0.25);
        stats.add(-0.75);
        assertErrorStatistics(stats, 0.75, Math.sqrt((0.25 + 0.75 * 0.75) / 3), -0.25 / 3);
        stats.add(-0.25);
        assertErrorStatistics(stats, 0.75, Math.sqrt((0.25 + 0.75 * 0.75 + 0.25 * 0.25) / 4), -0.5 / 4);
    }

    private static void assertErrorStatistics(ErrorStatistics stats, double max, double rms, double mean) {
        Assertions.assertEquals(max, stats.getMaxAbs(), "max absolute");
        Assertions.assertEquals(rms, stats.getRMS(), "rms");
        Assertions.assertEquals(mean, stats.getMean(), "mean");
    }

    @Test
    void testErrorStatisticsMean() {
        final ErrorStatistics stats = new ErrorStatistics();
        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements
        // Peters' example
        stats.add(1);
        stats.add(1e100);
        stats.add(1);
        stats.add(-1e100);
        Assertions.assertEquals(2.0 / 4, stats.getMean());
    }
}
