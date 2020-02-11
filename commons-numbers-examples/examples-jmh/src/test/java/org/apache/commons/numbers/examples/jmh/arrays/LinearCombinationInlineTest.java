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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.SplittableRandom;
import java.util.function.ToDoubleBiFunction;

/**
 * Test the inline version of each variant of the LinearCombination class.
 * This test ensure the speed test on the inlined versions is testing a valid
 * implementation of the generic array method.
 */
public class LinearCombinationInlineTest {
    /**
     * Define a scalar product between arrays of length 2.
     */
    private interface Scalar2 {
        /**
         * Compute the sum of the products of two sequences of factors.
         * @param a1 First factor of the first term.
         * @param b1 Second factor of the first term.
         * @param a2 First factor of the second term.
         * @param b2 Second factor of the second term.
         * @return \( a_1 b_1 + a_2 b_2 \)
         */
        double apply(double a1, double b1,
                     double a2, double b2);
    }

    /**
     * Define a scalar product between arrays of length 3.
     */
    private interface Scalar3 {
        /**
         * Compute the sum of the products of two sequences of factors.
         * @param a1 First factor of the first term.
         * @param b1 Second factor of the first term.
         * @param a2 First factor of the second term.
         * @param b2 Second factor of the second term.
         * @param a3 First factor of the third term.
         * @param b3 Second factor of the third term.
         * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 \)
         */
        double apply(double a1, double b1,
                     double a2, double b2,
                     double a3, double b3);
    }

    /**
     * Define a scalar product between arrays of length 4.
     */
    private interface Scalar4 {
        /**
         * Compute the sum of the products of two sequences of factors.
         * @param a1 First factor of the first term.
         * @param b1 Second factor of the first term.
         * @param a2 First factor of the second term.
         * @param b2 Second factor of the second term.
         * @param a3 First factor of the third term.
         * @param b3 Second factor of the third term.
         * @param a4 First factor of the fourth term.
         * @param b4 Second factor of the fourth term.
         * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 + a_4 b_4 \)
         */
        double apply(double a1, double b1,
                     double a2, double b2,
                     double a3, double b3,
                     double a4, double b4);
    }

    @Test
    public void testBigDecimal() {
        testArrayVsInline(LinearCombinationBigDecimal::value,
                          LinearCombinationBigDecimal::value,
                          LinearCombinationBigDecimal::value,
                          LinearCombinationBigDecimal::value);
    }

    @Test
    public void testDekker() {
        testArrayVsInline(LinearCombinationDekker::value,
                          LinearCombinationDekker::value,
                          LinearCombinationDekker::value,
                          LinearCombinationDekker::value);
    }

//    @Test
//    public void testDot2s() {
//        testArrayVsInline(LinearCombinationDot2s::value,
//                          LinearCombinationDot2s::value,
//                          LinearCombinationDot2s::value,
//                          LinearCombinationDot2s::value);
//    }

    @Test
    public void testDot2sSplitMask() {
        testArrayVsInline(LinearCombinationDot2sSplitMask::value,
                          LinearCombinationDot2sSplitMask::value,
                          LinearCombinationDot2sSplitMask::value,
                          LinearCombinationDot2sSplitMask::value);
    }

    @Test
    public void testDotK() {
        testArrayVsInline(LinearCombinationDotK::value,
                          LinearCombinationDotK::value,
                          LinearCombinationDotK::value,
                          LinearCombinationDotK::value);
    }

    private static void testArrayVsInline(ToDoubleBiFunction<double[], double[]> fun,
                                          Scalar2 fun2,
                                          Scalar3 fun3,
                                          Scalar4 fun4) {
        final SplittableRandom rng = new SplittableRandom();

        double sInline;
        double sArray;
        final double scale = 1e17;
        for (int i = 0; i < 100; ++i) {
            final double u1 = scale * rng.nextDouble();
            final double u2 = scale * rng.nextDouble();
            final double u3 = scale * rng.nextDouble();
            final double u4 = scale * rng.nextDouble();
            final double v1 = scale * rng.nextDouble();
            final double v2 = scale * rng.nextDouble();
            final double v3 = scale * rng.nextDouble();
            final double v4 = scale * rng.nextDouble();

            // One sum.
            sInline = fun2.apply(u1, v1, u2, v2);
            sArray = fun.applyAsDouble(new double[] {u1, u2},
                                       new double[] {v1, v2});
            Assertions.assertEquals(sInline, sArray);

            // Two sums.
            sInline = fun3.apply(u1, v1, u2, v2, u3, v3);
            sArray = fun.applyAsDouble(new double[] {u1, u2, u3},
                                       new double[] {v1, v2, v3});
            Assertions.assertEquals(sInline, sArray);

            // Three sums.
            sInline = fun4.apply(u1, v1, u2, v2, u3, v3, u4, v4);
            sArray = fun.applyAsDouble(new double[] {u1, u2, u3, u4},
                                       new double[] {v1, v2, v3, v4});
            Assertions.assertEquals(sInline, sArray);
        }
    }
}
