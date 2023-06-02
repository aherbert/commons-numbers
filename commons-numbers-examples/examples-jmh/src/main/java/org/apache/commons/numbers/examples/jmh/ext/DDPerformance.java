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

// License for the Boost continued fraction adaptation:

//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

package org.apache.commons.numbers.examples.jmh.ext;

import java.util.concurrent.TimeUnit;
import java.util.function.BinaryOperator;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Executes a benchmark to estimate the speed of double-double extended precision number
 * implementations.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgs = {"-server", "-Xms512M", "-Xmx512M"})
public class DDPerformance {

    /** Static mutable DD implementation. */
    static final String IMP_DD_STATIC_MUTABLE = "static mutable";
    /** Static immutable DD implementation. */
    static final String IMP_DD_STATIC_IMMUTABLE = "static immutable";
    /** OO immutable DD implementation. */
    static final String IMP_DD_OO_IMMUTABLE = "OO immutable";

    /**
     * Interface for an {@code (double, int) -> double} function.
     */
    public interface DoubleIntFunction {
        /**
         * Apply the function.
         *
         * @param x Argument.
         * @param n Argument.
         * @return the result
         */
        double apply(double x, int n);
    }

    /**
     * A {@code (double, int)} tuple.
     */
    public static class DoubleInt {
        /** double value. */
        private final double x;
        /** int value. */
        private final int n;

        /**
         * @param x double value
         * @param n int value
         */
        DoubleInt(double x, int n) {
            this.x = x;
            this.n = n;
        }

        /**
         * @return x
         */
        double getX() {
            return x;
        }

        /**
         * @return n
         */
        int getN() {
            return n;
        }
    }

    /**
     * Contains the function to computes the complementary probability {@code P[D_n^+ >= x]}
     * for the one-sided one-sample Kolmogorov-Smirnov distribution.
     */
    @State(Scope.Benchmark)
    public static class KSMethod {
        /** The implementation of the function. */
        @Param({IMP_DD_STATIC_MUTABLE, IMP_DD_STATIC_IMMUTABLE, IMP_DD_OO_IMMUTABLE})
        private String implementation;

        /** The function. */
        private DoubleIntFunction function;

        /**
         * Gets the function.
         *
         * @return the function
         */
        public DoubleIntFunction getFunction() {
            return function;
        }

        /**
         * Create the function.
         */
        @Setup
        public void setup() {
            function = createFunction(implementation);
        }

        /**
         * Creates the function to evaluate the one-sided one-sample Kolmogorov-Smirnov distribution.
         *
         * @param implementation Function implementation
         * @return the function
         */
        static DoubleIntFunction createFunction(String implementation) {
            if (IMP_DD_STATIC_MUTABLE.equals(implementation)) {
                return KolmogorovSmirnovDistribution.One::sf;
            } else if (IMP_DD_STATIC_IMMUTABLE.equals(implementation)) {
                return KolmogorovSmirnovDistribution.One::sf2;
            } else if (IMP_DD_OO_IMMUTABLE.equals(implementation)) {
                return KolmogorovSmirnovDistribution.One::sf3;
            } else {
                throw new IllegalStateException("unknown KS method: " + implementation);
            }
        }
    }

    /**
     * Contains the data to computes the complementary probability {@code P[D_n^+ >= x]}
     * for the one-sided one-sample Kolmogorov-Smirnov distribution.
     */
    @State(Scope.Benchmark)
    public static class KSData {
        // The parameters should be chosen such that the computation takes less than 1 second
        // thus allowing for repeats in the iteration.
        // Maintain n*x*x < 372.5 and n*x > 3 (see KSSample for details).
        // n=10000, values=50, ux=0.15
        // n=100000, values=5, ux=0.05

        /** The sample size for the KS distribution. This should be below the large N limit of 1000000. */
        @Param({"1000"})
        private int n;
        /** The number of values. */
        @Param({"500"})
        private int values;
        /** The lower limit on x. */
        @Param({"0.001"})
        private double lx;
        /** The upper limit on x. */
        @Param({"0.5"})
        private double ux;

        /** The data. */
        private DoubleInt[] data;

        /**
         * Gets the data.
         *
         * @return the data
         */
        public DoubleInt[] getData() {
            return data;
        }

        /**
         * Create the function.
         */
        @Setup
        public void setup() {
            data = createData(n, values, lx, ux);
        }

        /**
         * Creates the data for the one-sided one-sample Kolmogorov-Smirnov distribution.
         * The value {@code x} should by in the range [0, 1].
         *
         * @param n KS sample size
         * @param values Number of values
         * @param lx Lower limit on x
         * @param ux Upper limit on x
         * @return the data
         */
        static DoubleInt[] createData(int n, int values, double lx, double ux) {
            assert n > 0 : "Invalid n";
            assert lx <= ux : "Invalid range";
            if (values <= 1) {
                // Single value
                return new DoubleInt[] {new DoubleInt((lx + ux) * 0.5, n)};
            }
            // Create values between the lower and upper range
            final double inc = (ux - lx) / (values - 1);
            return IntStream.range(0, values)
                            .mapToObj(i -> new DoubleInt(lx + inc * i, n))
                            .toArray(DoubleInt[]::new);
        }
    }

    /**
     * Contains the data to computes the complementary probability {@code P[D_n^+ >= x]}
     * for the one-sided one-sample Kolmogorov-Smirnov distribution.
     */
    @State(Scope.Benchmark)
    public static class KSSample {
        /** The sample size for the KS distribution. This should be below the large N limit of 1000000. */
        @Param({"10000", "100000"})
        private int n;
        /**
         * The KS value (in the range [0, 1].
         *
         * <p>Note that the use of the full computation depends on the parameters.
         * If {@code n*x*x >= 372.5} then the p-value underflows. So large n limits the
         * usable upper range of x. If {@code n*x <= 3} then a faster computation can be performed
         * (either Smirnov-Dwass or exact when {@code nx <= 1}).
         */
        @Param({"0.01", "0.05"})
        private double x;

        /**
         * @return x
         */
        double getX() {
            return x;
        }

        /**
         * @return n
         */
        int getN() {
            return n;
        }
    }

    /**
     * Contains the data to compute the double-double operations.
     */
    @State(Scope.Benchmark)
    public static class OperatorData {
        /** The sample size. */
        @Param({"1000"})
        private int n;

        /** The data. */
        private DD3[] data;
        /** The second data. */
        private DD3[] data2;

        /**
         * Gets the data.
         *
         * @return the data
         */
        public DD3[] getData() {
            return data;
        }

        /**
         * Gets the second data.
         *
         * @return the second data
         */
        public DD3[] getData2() {
            return data2;
        }

        /**
         * Create the data.
         */
        @Setup(Level.Iteration)
        public void setup() {
            data = createData(n);
            data2 = createData(n);
        }

        /**
         * Creates data where the high part is approximately uniform in the range [-1, 1).
         * The actual value may exceed the range by half a ulp due to the low part of the
         * double-double number.
         *
         * @param n sample size
         * @return the data
         */
        static DD3[] createData(int n) {
            UniformRandomProvider rng = RandomSource.XO_RO_SHI_RO_128_PP.create();
            return IntStream.range(0, n)
                            .mapToObj(i -> {
                                // Uniform in [-1, 1) using increments of 2^-53
                                double x = makeSignedDouble(rng.nextLong());
                                // Uniform in [-2^-53, 2^-53)
                                double y = makeSignedDouble(rng.nextLong()) * 0x1.0p-53;
                                return DD3.twoSum(x, y);
                            }).toArray(DD3[]::new);
        }

        /**
         * Creates a signed double in the range {@code [-1, 1)}. The magnitude is sampled evenly from the
         * 2<sup>54</sup> dyadic rationals in the range.
         *
         * <p>Note: This method will not return samples for both -0.0 and 0.0.
         *
         * @param bits the bits
         * @return the double
         */
        private static double makeSignedDouble(long bits) {
            // Use the upper 54 bits on the assumption they are more random.
            // The sign bit is maintained by the signed shift.
            // The next 53 bits generates a magnitude in the range [0, 2^53) or [-2^53,h 0).
            return (bits >> 10) * 0x1.0p-53d;
        }
    }

    /**
     * Contains the data to compute the double-double operations.
     */
    @State(Scope.Benchmark)
    public static class BinaryOperatorMethod {
        /** The implementation of the function. */
        @Param({"fastAdd", "add", "uncheckedMultiply", "multiply"})
        private String implementation;

        /** The function. */
        private BinaryOperator<DD3> function;

        /**
         * Gets the function.
         *
         * @return the function
         */
        public BinaryOperator<DD3> getFunction() {
            return function;
        }

        /**
         * Create the function.
         */
        @Setup
        public void setup() {
            function = createFunction(implementation);
        }

        /**
         * Creates the function to evaluate the one-sided one-sample Kolmogorov-Smirnov distribution.
         *
         * @param implementation Function implementation
         * @return the function
         */
        static BinaryOperator<DD3> createFunction(String implementation) {
            if ("fastAdd".equals(implementation)) {
                return DD3::fastAdd;
            } else if ("add".equals(implementation)) {
                return DD3::add;
            } else if ("uncheckedMultiply".equals(implementation)) {
                return DD3::uncheckedMultiply;
            } else if ("multiply".equals(implementation)) {
                return DD3::multiply;
            } else {
                throw new IllegalStateException("unknown binary operator: " + implementation);
            }
        }
    }

    /**
     * Gets the double-double implementations.
     *
     * @return the implementations
     */
    static Stream<String> getImplementations() {
        return Stream.of(IMP_DD_STATIC_MUTABLE, IMP_DD_STATIC_IMMUTABLE, IMP_DD_OO_IMMUTABLE);
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param fun Function.
     * @param data Data.
     * @param bh Data sink.
     */
    private static void apply(DoubleIntFunction fun, DoubleInt[] data, Blackhole bh) {
        for (final DoubleInt d : data) {
            bh.consume(fun.apply(d.getX(), d.getN()));
        }
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param fun Function.
     * @param data Data.
     * @param data2 Second data.
     * @param bh Data sink.
     */
    private static void apply(BinaryOperator<DD3> fun, DD3[] data, DD3[] data2, Blackhole bh) {
        for (int i = 0; i < data.length; i++) {
            bh.consume(fun.apply(data[i], data2[i]));
        }
    }

    // Benchmark methods.
    // Benchmarks use function references to perform different operations on the numbers.

    /**
     * Benchmark a range of the KS function.
     *
     * @param method Test method.
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void ksRange(KSMethod method, KSData data, Blackhole bh) {
        apply(method.getFunction(), data.getData(), bh);
    }

    /**
     * Benchmark a sample of the KS function.
     *
     * @param method Test method.
     * @param data Test data.
     * @return the sample value
     */
    @Benchmark
    public double ksSample(KSMethod method, KSSample data) {
        return method.getFunction().apply(data.getX(), data.getN());
    }

    /**
     * Benchmark a sample of the KS function.
     *
     * @param method Test method.
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void binaryOperator(BinaryOperatorMethod method, OperatorData data, Blackhole bh) {
        apply(method.getFunction(), data.getData(), data.getData2(), bh);
    }
}
