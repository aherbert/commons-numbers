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
package org.apache.commons.numbers.examples.jmh.arrays;

import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.FourD;
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.ND;
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.ThreeD;
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.TwoD;
import org.apache.commons.numbers.fraction.BigFraction;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.SplittableRandom;
import java.util.stream.Stream;

/**
 * Test each implementation of the LinearCombination interface.
 */
public class LinearCombinationsTest {
    /**
     * Provide instances of the LinearCombination interface as arguments.
     *
     * @return the stream
     */
    static Stream<Arguments> provideLinearCombination() {
        return Stream.of(
            Arguments.of(LinearCombinations.Dekker.INSTANCE),
            Arguments.of(LinearCombinations.DotK.DOT_2),
            Arguments.of(LinearCombinations.DotK.DOT_3),
            Arguments.of(LinearCombinations.DotK.DOT_4),
            Arguments.of(LinearCombinations.DotK.DOT_5),
            Arguments.of(LinearCombinations.DotK.DOT_6),
            Arguments.of(LinearCombinations.DotK.DOT_7),
            Arguments.of(LinearCombinations.Exact.INSTANCE),
            Arguments.of(LinearCombinations.ExtendedPrecision.INSTANCE)
        );
    }

    @ParameterizedTest
    @MethodSource("provideLinearCombination")
    void testSingleElementArray(ND fun) {
        final double[] a = {1.23456789};
        final double[] b = {98765432.1};

        Assertions.assertEquals(a[0] * b[0], fun.value(a, b));
    }

    @ParameterizedTest
    @MethodSource("provideLinearCombination")
    void testTwoSums(ND fun) {
        final BigFraction[] aF = new BigFraction[] {
            BigFraction.of(-1321008684645961L, 268435456L),
            BigFraction.of(-5774608829631843L, 268435456L),
            BigFraction.of(-7645843051051357L, 8589934592L)
        };
        final BigFraction[] bF = new BigFraction[] {
            BigFraction.of(-5712344449280879L, 2097152L),
            BigFraction.of(-4550117129121957L, 2097152L),
            BigFraction.of(8846951984510141L, 131072L)
        };
        final ThreeD fun3 = (ThreeD) fun;

        final int len = aF.length;
        final double[] a = new double[len];
        final double[] b = new double[len];
        for (int i = 0; i < len; i++) {
            a[i] = aF[i].getNumerator().doubleValue() / aF[i].getDenominator().doubleValue();
            b[i] = bF[i].getNumerator().doubleValue() / bF[i].getDenominator().doubleValue();
        }

        // Ensure "array" and "inline" implementations give the same result.
        final double abSumInline = fun3.value(a[0], b[0],
                                              a[1], b[1],
                                              a[2], b[2]);
        final double abSumArray = fun.value(a, b);
        Assertions.assertEquals(abSumInline, abSumArray);

        // Compare with arbitrary precision computation.
        BigFraction result = BigFraction.ZERO;
        for (int i = 0; i < a.length; i++) {
            result = result.add(aF[i].multiply(bF[i]));
        }
        final double expected = result.doubleValue();
        Assertions.assertEquals(expected, abSumInline, "Expecting exact result");

        final double naive = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        Assertions.assertTrue(Math.abs(naive - abSumInline) > 1.5);
    }

    @ParameterizedTest
    @MethodSource("provideLinearCombination")
    void testHuge(ND fun) {
        final int scale = 971;
        final double[] a = new double[] {
            -1321008684645961.0 / 268435456.0,
            -5774608829631843.0 / 268435456.0,
            -7645843051051357.0 / 8589934592.0
        };
        final double[] b = new double[] {
            -5712344449280879.0 / 2097152.0,
            -4550117129121957.0 / 2097152.0,
            8846951984510141.0 / 131072.0
        };
        final ThreeD fun3 = (ThreeD) fun;

        final int len = a.length;
        final double[] scaledA = new double[len];
        final double[] scaledB = new double[len];
        for (int i = 0; i < len; ++i) {
            scaledA[i] = Math.scalb(a[i], -scale);
            scaledB[i] = Math.scalb(b[i], scale);
        }
        final double abSumInline = fun3.value(scaledA[0], scaledB[0],
                                              scaledA[1], scaledB[1],
                                              scaledA[2], scaledB[2]);
        final double abSumArray = fun.value(scaledA, scaledB);

        Assertions.assertEquals(abSumInline, abSumArray);
        Assertions.assertEquals(-1.8551294182586248737720779899, abSumInline, "Expecting exact result");

        final double naive = scaledA[0] * scaledB[0] + scaledA[1] * scaledB[1] + scaledA[2] * scaledB[2];
        Assertions.assertTrue(Math.abs(naive - abSumInline) > 1.5);
    }

    @ParameterizedTest
    @MethodSource("provideLinearCombination")
    void testArrayVsInline(ND fun) {
        // Assume the instance implements the inline functions
        final TwoD fun2 = (TwoD) fun;
        final ThreeD fun3 = (ThreeD) fun;
        final FourD fun4 = (FourD) fun;

        final SplittableRandom rng = new SplittableRandom();

        double sInline;
        double sArray;
        final double scale = 1e17;
        for (int i = 0; i < 1000; ++i) {
            final double u1 = scale * rng.nextDouble();
            final double u2 = scale * rng.nextDouble();
            final double u3 = scale * rng.nextDouble();
            final double u4 = scale * rng.nextDouble();
            final double v1 = scale * rng.nextDouble();
            final double v2 = scale * rng.nextDouble();
            final double v3 = scale * rng.nextDouble();
            final double v4 = scale * rng.nextDouble();

            // One sum.
            sInline = fun2.value(u1, v1, u2, v2);
            sArray = fun.value(new double[] {u1, u2},
                               new double[] {v1, v2});
            Assertions.assertEquals(sInline, sArray);

            // Two sums.
            sInline = fun3.value(u1, v1, u2, v2, u3, v3);
            sArray = fun.value(new double[] {u1, u2, u3},
                               new double[] {v1, v2, v3});
            Assertions.assertEquals(sInline, sArray);

            // Three sums.
            sInline = fun4.value(u1, v1, u2, v2, u3, v3, u4, v4);
            sArray = fun.value(new double[] {u1, u2, u3, u4},
                               new double[] {v1, v2, v3, v4});
            Assertions.assertEquals(sInline, sArray);
        }
    }

    @ParameterizedTest
    @MethodSource("provideLinearCombination")
    void testNonFinite(ND fun) {
        final double[][] a = new double[][] {
            {1, 2, 3, 4},
            {1, Double.POSITIVE_INFINITY, 3, 4},
            {1, 2, Double.POSITIVE_INFINITY, 4},
            {1, Double.POSITIVE_INFINITY, 3, Double.NEGATIVE_INFINITY},
            {1, 2, 3, 4},
            {1, 2, 3, 4},
            {1, 2, 3, 4},
            {1, 2, 3, 4},
            {1, Double.MAX_VALUE, 3, 4},
            {1, 2, Double.MAX_VALUE, 4},
            {1, Double.MAX_VALUE / 2, 3, -Double.MAX_VALUE / 4},
        };
        final double[][] b = new double[][] {
            {1, -2, 3, 4},
            {1, -2, 3, 4},
            {1, -2, 3, 4},
            {1, -2, 3, 4},
            {1, Double.POSITIVE_INFINITY, 3, 4},
            {1, -2, Double.POSITIVE_INFINITY, 4},
            {1, Double.POSITIVE_INFINITY, 3, Double.NEGATIVE_INFINITY},
            {Double.NaN, -2, 3, 4},
            {1, -2, 3, 4},
            {1, -2, 3, 4},
            {1, -2, 3, 4},
        };

        final TwoD fun2 = (TwoD) fun;
        final ThreeD fun3 = (ThreeD) fun;
        final FourD fun4 = (FourD) fun;

        Assertions.assertEquals(-3, fun2.value(a[0][0], b[0][0],
                                               a[0][1], b[0][1]));
        Assertions.assertEquals(6, fun3.value(a[0][0], b[0][0],
                                              a[0][1], b[0][1],
                                              a[0][2], b[0][2]));
        Assertions.assertEquals(22, fun4.value(a[0][0], b[0][0],
                                               a[0][1], b[0][1],
                                               a[0][2], b[0][2],
                                               a[0][3], b[0][3]));
        Assertions.assertEquals(22, fun.value(a[0], b[0]));

        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun2.value(a[1][0], b[1][0],
                                                                     a[1][1], b[1][1]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun3.value(a[1][0], b[1][0],
                                                                     a[1][1], b[1][1],
                                                                     a[1][2], b[1][2]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun4.value(a[1][0], b[1][0],
                                                                     a[1][1], b[1][1],
                                                                     a[1][2], b[1][2],
                                                                     a[1][3], b[1][3]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun.value(a[1], b[1]));

        Assertions.assertEquals(-3, fun2.value(a[2][0], b[2][0],
                                               a[2][1], b[2][1]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun3.value(a[2][0], b[2][0],
                                                                     a[2][1], b[2][1],
                                                                     a[2][2], b[2][2]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun4.value(a[2][0], b[2][0],
                                                                     a[2][1], b[2][1],
                                                                     a[2][2], b[2][2],
                                                                     a[2][3], b[2][3]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun.value(a[2], b[2]));

        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun2.value(a[3][0], b[3][0],
                                                                     a[3][1], b[3][1]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun3.value(a[3][0], b[3][0],
                                                                     a[3][1], b[3][1],
                                                                     a[3][2], b[3][2]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun4.value(a[3][0], b[3][0],
                                                                     a[3][1], b[3][1],
                                                                     a[3][2], b[3][2],
                                                                     a[3][3], b[3][3]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun.value(a[3], b[3]));

        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun2.value(a[4][0], b[4][0],
                                                                     a[4][1], b[4][1]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun3.value(a[4][0], b[4][0],
                                                                     a[4][1], b[4][1],
                                                                     a[4][2], b[4][2]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun4.value(a[4][0], b[4][0],
                                                                     a[4][1], b[4][1],
                                                                     a[4][2], b[4][2],
                                                                     a[4][3], b[4][3]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun.value(a[4], b[4]));

        Assertions.assertEquals(-3, fun2.value(a[5][0], b[5][0],
                                               a[5][1], b[5][1]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun3.value(a[5][0], b[5][0],
                                                                     a[5][1], b[5][1],
                                                                     a[5][2], b[5][2]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun4.value(a[5][0], b[5][0],
                                                                     a[5][1], b[5][1],
                                                                     a[5][2], b[5][2],
                                                                     a[5][3], b[5][3]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun.value(a[5], b[5]));

        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun2.value(a[6][0], b[6][0],
                                                                     a[6][1], b[6][1]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun3.value(a[6][0], b[6][0],
                                                                     a[6][1], b[6][1],
                                                                     a[6][2], b[6][2]));
        Assertions.assertEquals(Double.NaN, fun4.value(a[6][0], b[6][0],
                                                       a[6][1], b[6][1],
                                                       a[6][2], b[6][2],
                                                       a[6][3], b[6][3]));
        Assertions.assertEquals(Double.NaN, fun.value(a[6], b[6]));

        Assertions.assertEquals(Double.NaN, fun2.value(a[7][0], b[7][0],
                                                       a[7][1], b[7][1]));
        Assertions.assertEquals(Double.NaN, fun3.value(a[7][0], b[7][0],
                                                       a[7][1], b[7][1],
                                                       a[7][2], b[7][2]));
        Assertions.assertEquals(Double.NaN, fun4.value(a[7][0], b[7][0],
                                                       a[7][1], b[7][1],
                                                       a[7][2], b[7][2],
                                                       a[7][3], b[7][3]));
        Assertions.assertEquals(Double.NaN, fun.value(a[7], b[7]));

        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun2.value(a[8][0], b[8][0],
                                                                     a[8][1], b[8][1]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun3.value(a[8][0], b[8][0],
                                                                     a[8][1], b[8][1],
                                                                     a[8][2], b[8][2]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun4.value(a[8][0], b[8][0],
                                                                     a[8][1], b[8][1],
                                                                     a[8][2], b[8][2],
                                                                     a[8][3], b[8][3]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun.value(a[8], b[8]));

        Assertions.assertEquals(-3, fun2.value(a[9][0], b[9][0],
                                               a[9][1], b[9][1]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun3.value(a[9][0], b[9][0],
                                                                     a[9][1], b[9][1],
                                                                     a[9][2], b[9][2]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun4.value(a[9][0], b[9][0],
                                                                     a[9][1], b[9][1],
                                                                     a[9][2], b[9][2],
                                                                     a[9][3], b[9][3]));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, fun.value(a[9], b[9]));

        Assertions.assertEquals(-Double.MAX_VALUE, fun2.value(a[10][0], b[10][0],
                                                              a[10][1], b[10][1]));
        Assertions.assertEquals(-Double.MAX_VALUE, fun3.value(a[10][0], b[10][0],
                                                              a[10][1], b[10][1],
                                                              a[10][2], b[10][2]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun4.value(a[10][0], b[10][0],
                                                                     a[10][1], b[10][1],
                                                                     a[10][2], b[10][2],
                                                                     a[10][3], b[10][3]));
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, fun.value(a[10], b[10]));
    }
}
