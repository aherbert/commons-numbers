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

import org.apache.commons.numbers.arrays.LinearCombination;
import org.apache.commons.numbers.arrays.LinearCombination2;
import org.apache.commons.numbers.arrays.LinearCombinationExact;
import org.apache.commons.numbers.arrays.LinearCombinationExact2;
import org.apache.commons.numbers.arrays.LinearCombinationExact3;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

import java.math.MathContext;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.function.ToDoubleBiFunction;

/**
 * Executes a benchmark to measure the speed of operations in the {@link LinearCombination} class.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgs = {"-server", "-Xms512M", "-Xmx512M"})
public class LinearCombinationPerformance {
    /**
     * The seed to use to create the factors.
     * Using a fixed seed ensures the same factors are created for the variable
     * length arrays as for the small fixed size arrays.
     */
    private static final long SEED = System.currentTimeMillis();

    /**
     * The factors to multiply.
     */
    @State(Scope.Benchmark)
    public static class Factors {
        /**
         * The condition number of the generated data.
         */
        @Param({"1e20"})
        private double c;

        /**
         * The number of factors.
         */
        @Param({"1000"})
        private int size;

        /** Factors a. */
        private double[][] a;

        /** Factors b. */
        private double[][] b;

        /**
         * Gets the length of the array of factors.
         * This exists to be overridden by factors of a specific length.
         * The default is to create factors of length 4 for use in the inlined
         * scalar product methods.
         *
         * @return the length
         */
        public int getLength() {
            return 4;
        }

        /**
         * Gets the number of scalar products to compute.
         *
         * @return the size
         */
        public int getSize() {
            return size;
        }

        /**
         * Gets the a factors.
         *
         * @param index the index
         * @return Factors b.
         */
        public double[] getA(int index) {
            return a[index];
        }

        /**
         * Gets the b factors.
         *
         * @param index the index
         * @return Factors b.
         */
        public double[]  getB(int index) {
            return b[index];
        }

        /**
         * Create the factors.
         */
        @Setup
        public void setup() {
            final UniformRandomProvider rng =
                    RandomSource.create(RandomSource.XO_RO_SHI_RO_1024_PP, SEED);
            // Use the ill conditioned data generation method.
            // This requires an array of at least 6.
            final int n = Math.max(6, getLength());
            final double[] x = new double[n];
            final double[] y = new double[n];
            a = new double[size][];
            b = new double[size][];
            // Limit precision to allow large array lengths to be generated.
            final MathContext mathContext = new MathContext(100);
            for (int i = 0; i < size; i++) {
                LinearCombinationUtils.genDot(c, rng, x, y, null, mathContext);
                a[i] = Arrays.copyOf(x, getLength());
                b[i] = Arrays.copyOf(y, getLength());
            }
        }
    }

    /**
     * The factors to multiply of a specific length.
     */
    @State(Scope.Benchmark)
    public static class LengthFactors extends Factors {
        /**
         * The length of each factors array.
         */
        @Param({"2", "3", "4", "8", "16", "32", "64"})
        private int length;

        /** {@inheritDoc} */
        @Override
        public int getLength() {
            return length;
        }
    }

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

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param fun Scalar product function.
     * @return the scalar products
     */
    private static double[] scalar2Product(Factors factors, Scalar2 fun) {
        final double[] result = new double[factors.getSize()];
        for (int i = 0; i < result.length; i++) {
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            result[i] = fun.apply(a[0], b[0], a[1], b[1]);
        }
        return result;
    }

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param fun Scalar product function.
     * @return the scalar products
     */
    private static double[] scalar3Product(Factors factors, Scalar3 fun) {
        final double[] result = new double[factors.getSize()];
        for (int i = 0; i < result.length; i++) {
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            result[i] = fun.apply(a[0], b[0], a[1], b[1], a[2], b[2]);
        }
        return result;
    }

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param fun Scalar product function.
     * @return the scalar products
     */
    private static double[] scalar4Product(Factors factors, Scalar4 fun) {
        final double[] result = new double[factors.getSize()];
        for (int i = 0; i < result.length; i++) {
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            result[i] = fun.apply(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3]);
        }
        return result;
    }

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param fun Scalar product function.
     * @return the scalar products
     */
    private static double[] scalarProduct(LengthFactors factors, ToDoubleBiFunction<double[], double[]> fun) {
        final double[] result = new double[factors.getSize()];
        for (int i = 0; i < result.length; i++) {
            // These should be pre-computed to the correct length
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            result[i] = fun.applyAsDouble(a, b);
        }
        return result;
    }

    /**
     * Copy the original to a new shorter length.
     * Copy of {@link Arrays#copyOf(double[], int)} without the check for minimum length.
     *
     * @param original the original
     * @param newLength the new length
     * @return the copy
     */
    private static double[] copyOf(double[] original, int newLength) {
        double[] copy = new double[newLength];
        System.arraycopy(original, 0, copy, 0, newLength);
        return copy;
    }

    // Benchmark methods.
    //
    // The methods are partially documented as the names are self-documenting.
    // CHECKSTYLE: stop JavadocMethod
    // CHECKSTYLE: stop DesignForExtension

    /**
     * Baseline the standard precision scalar product of length 2.
     */
    @Benchmark
    public double[] scalar2Baseline(Factors factors) {
        return scalar2Product(factors, (a1, b1, a2, b2) -> a1 * b1 + a2 * b2);
    }

    /**
     * Baseline the standard precision scalar product of length 3.
     */
    @Benchmark
    public double[] scalar3Baseline(Factors factors) {
        return scalar3Product(factors, (a1, b1, a2, b2, a3, b3) -> a1 * b1 + a2 * b2 + a3 * b3);
    }

    /**
     * Baseline the standard precision scalar product of length 4.
     */
    @Benchmark
    public double[] scalar4Baseline(Factors factors) {
        return scalar4Product(factors, (a1, b1, a2, b2, a3, b3, a4, b4) -> a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4);
    }

    /**
     * Baseline the standard precision scalar product on two arrays.
     */
    @Benchmark
    public double[] scalarBaseline(LengthFactors factors) {
        return scalarProduct(factors, (a, b) -> {
            double sum = 0;
            for (int i = 0; i < a.length; i++) {
                sum += a[i] * b[i];
            }
            return sum;
        });
    }

    // Original method with Dekker split

    @Benchmark
    public double[] scalar2Dot2s(Factors factors) {
        return scalar2Product(factors, LinearCombinationDot2s::value);
    }

    @Benchmark
    public double[] scalar3Dot2s(Factors factors) {
        return scalar3Product(factors, LinearCombinationDot2s::value);
    }

    @Benchmark
    public double[] scalar4Dot2s(Factors factors) {
        return scalar4Product(factors, LinearCombinationDot2s::value);
    }

    @Benchmark
    public double[] scalarDot2s(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDot2s::value);
    }

    // Original method with masking split.
    // This is faster as it has no branching during the split to handle overflow.
    // The resulting split multiplication can have inexact results due to the loss of 1-bit.

    @Benchmark
    public double[] scalar2Dot2sMask(Factors factors) {
        return scalar2Product(factors, LinearCombinationDot2sSplitMask::value);
    }

    @Benchmark
    public double[] scalar3Dot2sMask(Factors factors) {
        return scalar3Product(factors, LinearCombinationDot2sSplitMask::value);
    }

    @Benchmark
    public double[] scalar4Dot2sMask(Factors factors) {
        return scalar4Product(factors, LinearCombinationDot2sSplitMask::value);
    }

    @Benchmark
    public double[] scalarDot2sMask(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDot2sSplitMask::value);
    }

    // dot3 using Dekker split

    @Benchmark
    public double[] scalar2Dot3(Factors factors) {
        return scalar2Product(factors, LinearCombination::value);
    }

    @Benchmark
    public double[] scalar3Dot3(Factors factors) {
        return scalar3Product(factors, LinearCombination::value);
    }

    @Benchmark
    public double[] scalar4Dot3(Factors factors) {
        return scalar4Product(factors, LinearCombination::value);
    }

    @Benchmark
    public double[] scalarDot3(LengthFactors factors) {
        return scalarProduct(factors, LinearCombination::value);
    }

    // dotK

    @Benchmark
    public double[] scalarDot4(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDotK::value4);
    }

    @Benchmark
    public double[] scalarDot5(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDotK::value5);
    }

    @Benchmark
    public double[] scalarDot6(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDotK::value6);
    }

    @Benchmark
    public double[] scalarDot7(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDotK::value7);
    }

    // dot3 using masking split
    // TODO - bring this into the JMH project

    @Benchmark
    public double[] scalar2Dot3Mask(Factors factors) {
        return scalar2Product(factors, LinearCombination2::value);
    }

    @Benchmark
    public double[] scalar3Dot3Mask(Factors factors) {
        return scalar3Product(factors, LinearCombination2::value);
    }

    @Benchmark
    public double[] scalar4Dot3Mask(Factors factors) {
        return scalar4Product(factors, LinearCombination2::value);
    }

    @Benchmark
    public double[] scalarDot3Mask(LengthFactors factors) {
        return scalarProduct(factors, LinearCombination2::value);
    }

    // Dekker's dot sum

    @Benchmark
    public double[] scalar2Dekker(Factors factors) {
        return scalar2Product(factors, LinearCombinationDekker::value);
    }

    @Benchmark
    public double[] scalar3Dekker(Factors factors) {
        return scalar3Product(factors, LinearCombinationDekker::value);
    }

    @Benchmark
    public double[] scalar4Dekker(Factors factors) {
        return scalar4Product(factors, LinearCombinationDekker::value);
    }

    @Benchmark
    public double[] scalarDekker(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationDekker::value);
    }

    // Exact expansion sum

    @Benchmark
    public double[] scalar2DotExact(Factors factors) {
        return scalar2Product(factors, LinearCombinationExact::value);
    }

    @Benchmark
    public double[] scalar3DotExact(Factors factors) {
        return scalar3Product(factors, LinearCombinationExact::value);
    }

    @Benchmark
    public double[] scalar4DotExact(Factors factors) {
        return scalar4Product(factors, LinearCombinationExact::value);
    }

    @Benchmark
    public double[] scalarDotExact(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationExact::value);
    }

    // Exact expansion sum (version 2)

    //@Benchmark
    public double[] scalar2DotExact2(Factors factors) {
        return scalar2Product(factors, LinearCombinationExact2::value);
    }

    //@Benchmark
    public double[] scalar3DotExact2(Factors factors) {
        return scalar3Product(factors, LinearCombinationExact2::value);
    }

    //@Benchmark
    public double[] scalar4DotExact2(Factors factors) {
        return scalar4Product(factors, LinearCombinationExact2::value);
    }

    //@Benchmark
    public double[] scalarDotExact2(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationExact2::value);
    }

    // Exact expansion sum (version 3)

    //@Benchmark
    public double[] scalar2DotExact3(Factors factors) {
        return scalar2Product(factors, LinearCombinationExact3::value);
    }

    //@Benchmark
    public double[] scalar3DotExact3(Factors factors) {
        return scalar3Product(factors, LinearCombinationExact3::value);
    }

    //@Benchmark
    public double[] scalar4DotExact3(Factors factors) {
        return scalar4Product(factors, LinearCombinationExact3::value);
    }

    @Benchmark
    public double[] scalarDotExact3(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationExact3::value);
    }

    // BigDecimal

    @Benchmark
    public double[] scalar2BigDecimal(Factors factors) {
        return scalar2Product(factors, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public double[] scalar3BigDecimal(Factors factors) {
        return scalar3Product(factors, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public double[] scalar4BigDecimal(Factors factors) {
        return scalar4Product(factors, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public double[] scalarBigDecimal(LengthFactors factors) {
        return scalarProduct(factors, LinearCombinationBigDecimal::value);
    }
}
