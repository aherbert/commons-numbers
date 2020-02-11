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

/**
 * Computes linear combinations accurately.
 * This class computes the sum of the products of two sequences of numbers
 * <code>a<sub>i</sub> b<sub>i</sub></code> to high accuracy.
 * It does so by using extended precision multiplication and addition algorithms to
 * preserve accuracy and reduce cancellation effects.
 *
 * <p>It is based on the 1971 paper
 * <a href="https://doi.org/10.1007/BF01397083">
 * Dekker (1971) A floating-point technique for extending the available precision</a>.
 *
 * <p>Due to the methods used to increase precision it is possible that the computation
 * overflows when the standard dot product does not. This may occur when an individual
 * product exceeds the magnitude of {@link Double#MAX_VALUE}{@code  / 2}.
 * In this case the standard precision result will be returned.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Dot_product">Dot product</a>
 */
public final class LinearCombinationDekker {
    /** Private constructor. */
    private LinearCombinationDekker() {
        // intentionally empty.
    }

    /**
     * Compute the sum of the products of two sequences of factors to 3-fold precision.
     *
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

        // Dekker's dot product method:
        // Dekker's multiply and then Dekker's add2 to sum in double length precision.

        double x = a[0] * b[0];
        double xx = DoublePrecision.productLow(a[0], b[0], x);

        // Remaining split products added to the current sum and round-off stored
        for (int i = 1; i < len; i++) {
            final double y = a[i] * b[i];
            final double yy = DoublePrecision.productLow(a[i], b[i], y);

            // sum this to the previous number using Dekker's add2 algorithm
            final double r = x + y;
            final double s = DoublePrecision.sumLow(x, xx, y, yy, r);

            x = r + s;
            // final split product does not require the round-off but included here for brevity
            xx = r - x + s;
        }

        if (Double.isNaN(x)) {
            // Either we have split infinite numbers or some coefficients were NaNs,
            // just rely on the naive implementation and let IEEE754 handle this
            x = 0;
            for (int i = 0; i < len; ++i) {
                x += a[i] * b[i];
            }
        }

        return x;
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
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
        final double x = a1 * b1;
        final double xx = DoublePrecision.productLow(a1, b1, x);
        final double y = a2 * b2;
        final double yy = DoublePrecision.productLow(a2, b2, y);
        double r = x + y;
        final double s = DoublePrecision.sumLow(x, xx, y, yy, r);
        r = r + s;

        if (Double.isNaN(r)) {
            // Either we have split infinite numbers or some coefficients were NaNs,
            // just rely on the naive implementation and let IEEE754 handle this
            return x + y;
        }

        return r;
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
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
        double x = a1 * b1;
        double xx = DoublePrecision.productLow(a1, b1, x);
        double y = a2 * b2;
        double yy = DoublePrecision.productLow(a2, b2, y);
        double r = x + y;
        double s = DoublePrecision.sumLow(x, xx, y, yy, r);
        x = r + s;
        xx = r - x + s;
        y = a3 * b3;
        yy = DoublePrecision.productLow(a3, b3, y);
        r = x + y;
        s = DoublePrecision.sumLow(x, xx, y, yy, r);
        x = r + s;

        if (Double.isNaN(x)) {
            // Either we have split infinite numbers or some coefficients were NaNs,
            // just rely on the naive implementation and let IEEE754 handle this
            return a1 * b1 + a2 * b2 + y;
        }

        return x;
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
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
        double x = a1 * b1;
        double xx = DoublePrecision.productLow(a1, b1, x);
        double y = a2 * b2;
        double yy = DoublePrecision.productLow(a2, b2, y);
        double r = x + y;
        double s = DoublePrecision.sumLow(x, xx, y, yy, r);
        x = r + s;
        xx = r - x + s;
        y = a3 * b3;
        yy = DoublePrecision.productLow(a3, b3, y);
        r = x + y;
        s = DoublePrecision.sumLow(x, xx, y, yy, r);
        x = r + s;
        xx = r - x + s;
        y = a4 * b4;
        yy = DoublePrecision.productLow(a4, b4, y);
        r = x + y;
        s = DoublePrecision.sumLow(x, xx, y, yy, r);
        x = r + s;

        if (Double.isNaN(x)) {
            // Either we have split infinite numbers or some coefficients were NaNs,
            // just rely on the naive implementation and let IEEE754 handle this
            return a1 * b1 + a2 * b2 + a3 * b3 + y;
        }

        return x;
    }
}
