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

import java.math.BigDecimal;

/**
 * Computes linear combinations accurately.
 * This method computes the sum of the products
 * <code>a<sub>i</sub> b<sub>i</sub></code> to high accuracy
 * using BigDecimal.
 */
public final class LinearCombinationBigDecimal {
    /** Private constructor. */
    private LinearCombinationBigDecimal() {
        // intentionally empty.
    }

    /**
     * @param a Factors.
     * @param b Factors.
     * @return \( \sum_i a_i b_i \).
     * @throws IllegalArgumentException if the sizes of the arrays are different.
     */
    public static double value(double[] a,
                               double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("Dimension mismatch: " + a.length + " != " + b.length);
        }

        final int len = a.length;

        if (len == 0) {
            return 0;
        }
        if (len == 1) {
            // Revert to scalar multiplication.
            return a[0] * b[0];
        }

        // BigDecimal cannot handle inf/nan so check the IEEE754 result
        double result = a[0] * b[0];
        for (int i = 1; i < len; i++) {
            result += a[i] * b[i];
        }
        if (!Double.isFinite(result)) {
            return result;
        }

        BigDecimal sum = multiply(a[0], b[0]);
        for (int i = 1; i < len; i++) {
            sum = sum.add(multiply(a[i], b[i]));
        }
        return sum.doubleValue();
    }

    /**
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @return \( a_1 b_1 + a_2 b_2 \)
     *
     * @see #value(double, double, double, double, double, double)
     * @see #value(double, double, double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2) {
        // BigDecimal cannot handle inf/nan so check the IEEE754 result
        final double result = a1 * b1 + a2 * b2;
        if (!Double.isFinite(result)) {
            return result;
        }
        return multiply(a1, b1)
                .add(multiply(a2, b2))
                .doubleValue();
    }

    /**
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param a3 First factor of the third term.
     * @param b3 Second factor of the third term.
     * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 \)
     *
     * @see #value(double, double, double, double)
     * @see #value(double, double, double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2,
                               double a3, double b3) {
        // BigDecimal cannot handle inf/nan so check the IEEE754 result
        final double result = a1 * b1 + a2 * b2 + a3 * b3;
        if (!Double.isFinite(result)) {
            return result;
        }
        return multiply(a1, b1)
                .add(multiply(a2, b2))
                .add(multiply(a3, b3))
                .doubleValue();
    }

    /**
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param a3 First factor of the third term.
     * @param b3 Second factor of the third term.
     * @param a4 First factor of the fourth term.
     * @param b4 Second factor of the fourth term.
     * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 + a_4 b_4 \)
     *
     * @see #value(double, double, double, double)
     * @see #value(double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2,
                               double a3, double b3,
                               double a4, double b4) {
        // BigDecimal cannot handle inf/nan so check the IEEE754 result
        final double result = a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
        if (!Double.isFinite(result)) {
            return result;
        }
        return multiply(a1, b1)
                .add(multiply(a2, b2))
                .add(multiply(a3, b3))
                .add(multiply(a4, b4))
                .doubleValue();
    }

    /**
     * Multiply the factors to a BigDecimal.
     *
     * @param a First factor.
     * @param b Second factor.
     * @return the product
     */
    private static BigDecimal multiply(double a, double b) {
        return new BigDecimal(a).multiply(new BigDecimal(b));
    }
}
