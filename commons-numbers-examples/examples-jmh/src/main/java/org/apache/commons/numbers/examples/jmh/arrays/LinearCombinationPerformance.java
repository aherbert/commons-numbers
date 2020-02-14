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
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.FourD;
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.ND;
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.ThreeD;
import org.apache.commons.numbers.examples.jmh.arrays.LinearCombination.TwoD;
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
import org.openjdk.jmh.infra.Blackhole;

import java.math.MathContext;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

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
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param bh Data sink.
     * @param fun Scalar product function.
     */
    private static void scalar2Product(Factors factors, Blackhole bh, TwoD fun) {
        for (int i = 0; i < factors.getSize(); i++) {
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            bh.consume(fun.value(a[0], b[0], a[1], b[1]));
        }
    }

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param bh Data sink.
     * @param fun Scalar product function.
     */
    private static void scalar3Product(Factors factors, Blackhole bh, ThreeD fun) {
        for (int i = 0; i < factors.getSize(); i++) {
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            bh.consume(fun.value(a[0], b[0], a[1], b[1], a[2], b[2]));
        }
    }

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param bh Data sink.
     * @param fun Scalar product function.
     */
    private static void scalar4Product(Factors factors, Blackhole bh, FourD fun) {
        for (int i = 0; i < factors.getSize(); i++) {
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            bh.consume(fun.value(a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3]));
        }
    }

    /**
     * Compute the scalar product for all the factors.
     *
     * @param factors Factors.
     * @param bh Data sink.
     * @param fun Scalar product function.
     */
    private static void scalarProduct(LengthFactors factors, Blackhole bh, ND fun) {
        for (int i = 0; i < factors.getSize(); i++) {
            // These should be pre-computed to the correct length
            final double[] a = factors.getA(i);
            final double[] b = factors.getB(i);
            bh.consume(fun.value(a, b));
        }
    }

    // Benchmark methods.
    //
    // The methods are partially documented as the names are self-documenting.
    // CHECKSTYLE: stop JavadocMethod
    // CHECKSTYLE: stop DesignForExtension

    /**
     * Baseline the standard precision scalar product of length 2.
     * @param factors Factors
     * @param bh Data sink.
     */
    @Benchmark
    public void scalar2Baseline(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, (a1, b1, a2, b2) -> a1 * b1 + a2 * b2);
    }

    /**
     * Baseline the standard precision scalar product of length 3.
     * @param factors Factors
     * @param bh Data sink.
     */
    @Benchmark
    public void scalar3Baseline(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, (a1, b1, a2, b2, a3, b3) -> a1 * b1 + a2 * b2 + a3 * b3);
    }

    /**
     * Baseline the standard precision scalar product of length 4.
     * @param factors Factors
     * @param bh Data sink.
     */
    @Benchmark
    public void scalar4Baseline(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, (a1, b1, a2, b2, a3, b3, a4, b4) -> a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4);
    }

    /**
     * Baseline the standard precision scalar product on two arrays.
     * @param factors Factors
     * @param bh Data sink.
     */
    @Benchmark
    public void scalarBaseline(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, (a, b) -> {
            double sum = 0;
            for (int i = 0; i < a.length; i++) {
                sum += a[i] * b[i];
            }
            return sum;
        });
    }

    // Original method with Dekker split

    @Benchmark
    public void scalar2Dot2s(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationDot2s::value);
    }

    @Benchmark
    public void scalar3Dot2s(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationDot2s::value);
    }

    @Benchmark
    public void scalar4Dot2s(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationDot2s::value);
    }

    @Benchmark
    public void scalarDot2s(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDot2s::value);
    }

    // Original method with masking split.
    // This is faster as it has no branching during the split to handle overflow.
    // The resulting split multiplication can have inexact results due to the loss of 1-bit.

    @Benchmark
    public void scalar2Dot2sMask(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationDot2sSplitMask::value);
    }

    @Benchmark
    public void scalar3Dot2sMask(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationDot2sSplitMask::value);
    }

    @Benchmark
    public void scalar4Dot2sMask(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationDot2sSplitMask::value);
    }

    @Benchmark
    public void scalarDot2sMask(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDot2sSplitMask::value);
    }

    // dot3 using Dekker split

    @Benchmark
    public void scalar2Dot3(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombination::value);
    }

    @Benchmark
    public void scalar3Dot3(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombination::value);
    }

    @Benchmark
    public void scalar4Dot3(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombination::value);
    }

    @Benchmark
    public void scalarDot3(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombination::value);
    }


    @Benchmark
    public void scalar2Dot3b(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinations.DotK.DOT_3);
    }

    @Benchmark
    public void scalar3Dot3b(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinations.DotK.DOT_3);
    }

    @Benchmark
    public void scalar4Dot3b(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinations.DotK.DOT_3);
    }

    @Benchmark
    public void scalarDot3b(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.DotK.DOT_3);
    }

    // dotK

    @Benchmark
    public void scalarDot4(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDotK::value4);
    }

    @Benchmark
    public void scalarDot5(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDotK::value5);
    }

    @Benchmark
    public void scalarDot6(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDotK::value6);
    }

    @Benchmark
    public void scalarDot7(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDotK::value7);
    }

    @Benchmark
    public void scalarDot4b(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.DotK.DOT_4);
    }

    @Benchmark
    public void scalarDot5b(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.DotK.DOT_5);
    }

    @Benchmark
    public void scalarDot6b(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.DotK.DOT_6);
    }

    @Benchmark
    public void scalarDot7b(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.DotK.DOT_7);
    }

    // Dekker's dot sum

    @Benchmark
    public void scalar2Dekker(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationDekker::value);
    }

    @Benchmark
    public void scalar3Dekker(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationDekker::value);
    }

    @Benchmark
    public void scalar4Dekker(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationDekker::value);
    }

    @Benchmark
    public void scalarDekker(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationDekker::value);
    }

    // Exact expansion sum

    @Benchmark
    public void scalar2DotExact(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationExact::value);
    }

    @Benchmark
    public void scalar3DotExact(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationExact::value);
    }

    @Benchmark
    public void scalar4DotExact(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationExact::value);
    }

    @Benchmark
    public void scalarDotExact(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationExact::value);
    }

    // Exact expansion sum (version 2)

    //@Benchmark
    public void scalar2DotExact2(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationExact2::value);
    }

    //@Benchmark
    public void scalar3DotExact2(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationExact2::value);
    }

    //@Benchmark
    public void scalar4DotExact2(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationExact2::value);
    }

    @Benchmark
    public void scalarDotExact2(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationExact2::value);
    }

    // Exact expansion sum (version 3)

    @Benchmark
    public void scalar2DotExact3(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationExact3::value);
    }

    @Benchmark
    public void scalar3DotExact3(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationExact3::value);
    }

    @Benchmark
    public void scalar4DotExact3(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationExact3::value);
    }

    @Benchmark
    public void scalarDotExact3(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationExact3::value);
    }

    @Benchmark
    public void scalar2DotExact3b(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinations.ExtendedPrecision.INSTANCE);
    }

    @Benchmark
    public void scalar3DotExact3b(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinations.ExtendedPrecision.INSTANCE);
    }

    @Benchmark
    public void scalar4DotExact3b(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinations.ExtendedPrecision.INSTANCE);
    }

    @Benchmark
    public void scalarDotExact3b(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.ExtendedPrecision.INSTANCE);
    }

    // BigDecimal

    @Benchmark
    public void scalar2BigDecimal(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public void scalar3BigDecimal(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public void scalar4BigDecimal(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public void scalarBigDecimal(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinationBigDecimal::value);
    }

    @Benchmark
    public void scalar2BigDecimalb(Factors factors, Blackhole bh) {
        scalar2Product(factors, bh, LinearCombinations.Exact.INSTANCE);
    }

    @Benchmark
    public void scalar3BigDecimalb(Factors factors, Blackhole bh) {
        scalar3Product(factors, bh, LinearCombinations.Exact.INSTANCE);
    }

    @Benchmark
    public void scalar4BigDecimalb(Factors factors, Blackhole bh) {
        scalar4Product(factors, bh, LinearCombinations.Exact.INSTANCE);
    }

    @Benchmark
    public void scalarBigDecimalb(LengthFactors factors, Blackhole bh) {
        scalarProduct(factors, bh, LinearCombinations.Exact.INSTANCE);
    }
}
