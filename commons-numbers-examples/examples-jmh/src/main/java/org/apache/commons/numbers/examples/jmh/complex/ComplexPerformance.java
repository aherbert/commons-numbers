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

package org.apache.commons.numbers.examples.jmh.complex;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.numbers.complex.Complex;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
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

import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.function.BiFunction;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;
import java.util.function.UnaryOperator;
import java.util.stream.Stream;

/**
 * Executes a benchmark to measure the speed of operations in the {@link Complex} class.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgs = {"-server", "-Xms512M", "-Xmx512M"})
public class ComplexPerformance {
    /**
     * An array of edge numbers that will produce edge case results from functions:
     * {@code +/-inf, +/-max, +/-min, +/-0, nan}.
     */
    private static final double[] EDGE_NUMBERS = {
        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.MAX_VALUE,
        -Double.MAX_VALUE, Double.MIN_VALUE, -Double.MIN_VALUE, 0.0, -0.0, Double.NaN};
    /** The range for the uniform random numbers. */
    private static final double RANGE = 1e53;
    /** The range for the log-uniform random numbers: 2^30. */
    private static final double LOG_RANGE = Math.log(0x1.0p30);
    /** The exponent for fixed exponent numbers. */
    private static final long FIXED_EXPONENT = Double.doubleToRawLongBits(1.0);

    /** The multiplier for the split. */
    private static final double MULTIPLIER = 1.34217729E8;

    /** 54 shifted 20-bits to align with the exponent of the upper 32-bits of a double. */
    private static final int EXP_54 = 0x36_00000;
    /** Represents an exponent of 500 in unbiased form shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_500 = 0x5f3_00000;
    /** Represents an exponent of 1024 in unbiased form (infinite or nan)
     * shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_1024 = 0x7ff_00000;
    /** Represents an exponent of -500 in unbiased form shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_NEG_500 = 0x20b_00000;
    /** Represents an exponent of -1022 in unbiased form (sub-normal number)
     * shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_NEG_1022 = 0x001_00000;
    /**
     * Contains the size of numbers.
     */
    @State(Scope.Benchmark)
    public static class ComplexNumberSize {
        /**
         * The size of the data.
         */
        @Param({"10000"})
        private int size;

        /**
         * Gets the size.
         *
         * @return the size
         */
        public int getSize() {
            return size;
        }
    }

    /**
     * Contains an array of complex numbers.
     */
    @State(Scope.Benchmark)
    public static class ComplexNumbers extends ComplexNumberSize {
        /** The numbers. */
        protected Complex[] numbers;

        /**
         * The type of the data.
         */
        @Param({"cis", "gaussian", "log-uniform", "uniform", "fixed-exponent", "edge"})
        private String type;

        /**
         * Gets the numbers.
         *
         * @return the numbers
         */
        public Complex[] getNumbers() {
            return numbers;
        }

        /**
         * Create the complex numbers.
         */
        @Setup
        public void setup() {
            numbers = createNumbers(RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP));
        }

        /**
         * Creates the numbers.
         *
         * @param rng Random number generator.
         * @return the random complex number
         */
        Complex[] createNumbers(UniformRandomProvider rng) {
            Supplier<Complex> generator;
            if ("cis".equals(type)) {
                generator = () -> Complex.ofCis(rng.nextDouble() * 2 * Math.PI);
            } else if ("gaussian".equals(type)) {
                final ZigguratNormalizedGaussianSampler s = ZigguratNormalizedGaussianSampler.of(rng);
                generator = () -> Complex.ofCartesian(s.sample(), s.sample());
            } else if ("log-uniform".equals(type)) {
                generator = () -> Complex.ofCartesian(createLogUniformNumber(rng), createLogUniformNumber(rng));
            } else if ("uniform".equals(type)) {
                generator = () -> Complex.ofCartesian(createUniformNumber(rng), createUniformNumber(rng));
            } else if ("fixed-exponent".equals(type)) {
                generator = () -> Complex.ofCartesian(createFixedExponentNumber(rng), createFixedExponentNumber(rng));
            } else if ("edge".equals(type)) {
                generator = () -> Complex.ofCartesian(createEdgeNumber(rng), createEdgeNumber(rng));
            } else {
                throw new IllegalStateException("Unknown number type: " + type);
            }
            return Stream.generate(generator).limit(getSize()).toArray(Complex[]::new);
        }
    }

    /**
     * Contains two arrays of complex numbers.
     */
    @State(Scope.Benchmark)
    public static class TwoComplexNumbers extends ComplexNumbers {
        /** The numbers. */
        private Complex[] numbers2;

        /**
         * Gets the second set of numbers.
         *
         * @return the numbers
         */
        public Complex[] getNumbers2() {
            return numbers2;
        }

        /**
         * Create the complex numbers.
         */
        @Override
        @Setup
        public void setup() {
            // Do not call super.setup() so we recycle the RNG and avoid duplicates
            final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
            numbers = createNumbers(rng);
            numbers2 = createNumbers(rng);
        }
    }

    /**
     * Contains an array of complex numbers and an array of real numbers.
     */
    @State(Scope.Benchmark)
    public static class ComplexAndRealNumbers extends ComplexNumbers {
        /** The numbers. */
        private double[] numbers2;

        /**
         * Gets the second set of numbers.
         *
         * @return the numbers
         */
        public double[] getNumbers2() {
            return numbers2;
        }

        /**
         * Create the complex numbers.
         */
        @Override
        @Setup
        public void setup() {
            // Do not call super.setup() so we recycle the RNG and avoid duplicates
            final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
            numbers = createNumbers(rng);
            numbers2 = Arrays.stream(createNumbers(rng)).mapToDouble(Complex::real).toArray();
        }
    }

    /**
     * Define a function between a complex and real number.
     */
    private interface ComplexRealFunction {
        /**
         * Applies this function to the given arguments.
         *
         * @param z the complex argument
         * @param x the real argument
         * @return the function result
         */
        Complex apply(Complex z, double x);
    }

    /**
     * Creates a random double number in a log-uniform range.
     *
     * @param rng Random number generator.
     * @return the random number
     */
    private static double createLogUniformNumber(UniformRandomProvider rng) {
        // The log-uniform distribution has an upper and lower bound.
        // Here we use the lower bound of 1 which has a logarithm of zero.
        // It can thus be ignored from: e^(uniform(ln(upper) - ln(lower))
        return Math.exp(rng.nextDouble() * LOG_RANGE);
    }

    /**
     * Creates a random double number with a random sign and a large range.
     * The numbers will be uniform over the range.
     *
     * @param rng Random number generator.
     * @return the random number
     */
    private static double createUniformNumber(UniformRandomProvider rng) {
        // Create [-1, 1) then multiply by a range
        return (rng.nextDouble() - rng.nextInt(1)) * RANGE;
    }

    /**
     * Creates a random double number with a fixed exponent and a random 52-bit mantissa.
     * The numbers will be uniform over the range.
     *
     * @param rng Random number generator.
     * @return the random number
     */
    private static double createFixedExponentNumber(UniformRandomProvider rng) {
        final long mask = ((1L << 52) - 1) | 1L << 63;
        final long bits = rng.nextLong() & mask;
        return bits | FIXED_EXPONENT;
    }

    /**
     * Creates a random double number that will be an edge case:
     * {@code +/-inf, +/-max, +/-min, +/-0, nan}.
     *
     * @param rng Random number generator.
     * @return the random number
     */
    private static double createEdgeNumber(UniformRandomProvider rng) {
        return EDGE_NUMBERS[rng.nextInt(EDGE_NUMBERS.length)];
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param numbers Numbers.
     * @param fun Function.
     * @return the result of the function.
     */
    private static boolean[] apply(Complex[] numbers, Predicate<Complex> fun) {
        final boolean[] result = new boolean[numbers.length];
        for (int i = 0; i < numbers.length; i++) {
            result[i] = fun.test(numbers[i]);
        }
        return result;
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param numbers Numbers.
     * @param fun Function.
     * @return the result of the function.
     */
    private static double[] apply(Complex[] numbers, ToDoubleFunction<Complex> fun) {
        final double[] result = new double[numbers.length];
        for (int i = 0; i < numbers.length; i++) {
            result[i] = fun.applyAsDouble(numbers[i]);
        }
        return result;
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param numbers Numbers.
     * @param fun Function.
     * @return the result of the function.
     */
    private static Complex[] apply(Complex[] numbers, UnaryOperator<Complex> fun) {
        final Complex[] result = new Complex[numbers.length];
        for (int i = 0; i < numbers.length; i++) {
            result[i] = fun.apply(numbers[i]);
        }
        return result;
    }

    /**
     * Apply the function to the paired numbers.
     *
     * @param numbers First numbers of the pairs.
     * @param numbers2 Second numbers of the pairs.
     * @param fun Function.
     * @return the result of the function.
     */
    private static Complex[] apply(Complex[] numbers, Complex[] numbers2,
            BiFunction<Complex, Complex, Complex> fun) {
        final Complex[] result = new Complex[numbers.length];
        for (int i = 0; i < numbers.length; i++) {
            result[i] = fun.apply(numbers[i], numbers2[i]);
        }
        return result;
    }

    /**
     * Apply the function to the paired numbers.
     *
     * @param numbers First numbers of the pairs.
     * @param numbers2 Second numbers of the pairs.
     * @param fun Function.
     * @return the result of the function.
     */
    private static Complex[] apply(Complex[] numbers, double[] numbers2,
            ComplexRealFunction fun) {
        final Complex[] result = new Complex[numbers.length];
        for (int i = 0; i < numbers.length; i++) {
            result[i] = fun.apply(numbers[i], numbers2[i]);
        }
        return result;
    }

    /**
     * Identity function. This can be used to measure overhead of object array creation.
     *
     * @param z Complex number.
     * @return the complex number
     */
    private static Complex identity(Complex z) {
        return z;
    }

    /**
     * Copy function. This can be used to measure overhead of object array creation plus
     * new Complex creation.
     *
     * @param z Complex number.
     * @return a copy of the complex number
     */
    private static Complex copy(Complex z) {
        return Complex.ofCartesian(z.real(), z.imag());
    }

    // Benchmark methods.
    //
    // The methods are partially documented as the names are self-documenting.
    // CHECKSTYLE: stop JavadocMethod
    // CHECKSTYLE: stop DesignForExtension
    //
    // Benchmarks use function references to perform different operations on the complex numbers.
    // Tests show that explicit programming of the same benchmarks run in the same time.
    // For reference examples are provided for the fastest operations: real() and conj().

    /**
     * Explicit benchmark without using a method reference.
     * This should run in the same time as {@link #real(ComplexNumbers)}.
     * This is commented out as it exists for reference purposes.
     */
    //@Benchmark
    public double[] real2(ComplexNumbers numbers) {
        final Complex[] z = numbers.getNumbers();
        final double[] result = new double[z.length];
        for (int i = 0; i < z.length; i++) {
            result[i] = z[i].real();
        }
        return result;
    }

    /**
     * Explicit benchmark without using a method reference.
     * This should run in the same time as {@link #conj(ComplexNumbers)}.
     * This is commented out as it exists for reference purposes.
     */
    //@Benchmark
    public Complex[] conj2(ComplexNumbers numbers) {
        final Complex[] z = numbers.getNumbers();
        final Complex[] result = new Complex[z.length];
        for (int i = 0; i < z.length; i++) {
            result[i] = z[i].conj();
        }
        return result;
    }

    /**
     * Baseline the creation of the new array of numbers.
     * This contains the baseline JMH overhead for all the benchmarks that create complex numbers.
     * All other methods are expected to be slower than this.
     */
    @Benchmark
    public Complex[] baselineNewArray(ComplexNumberSize numberSize) {
        return new Complex[numberSize.getSize()];
    }

    /**
     * Baseline the creation of a copy array of numbers.
     * This is commented out as it provides no information other than to demonstrate that
     * {@link #baselineCopy(ComplexNumbers)} is not being optimised to a single array copy
     * operation.
     */
    //@Benchmark
    public Complex[] baselineCopyArray(ComplexNumbers numbers) {
        return Arrays.copyOf(numbers.getNumbers(), numbers.getNumbers().length);
    }

    /**
     * Baseline the creation of the new array of numbers with the same complex number (an identity).
     *
     * <p>Note: This runs much faster than {@link #baselineCopy(ComplexNumbers)}. This is
     * attributed to the identity function not requiring that the fields of the
     * complex are accessed unlike all other methods that do computations on the real and/or
     * imaginary parts. The method is slower than a creation of a new empty array or a
     * copy array thus contains the loop overhead of the benchmarks that create new numbers.
     *
     * @see #baselineNewArray(ComplexNumberSize)
     * @see #baselineCopyArray(ComplexNumbers)
     */
    @Benchmark
    public Complex[] baselineIdentity(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), ComplexPerformance::identity);
    }

    /**
     * Baseline the creation of the new array of numbers with a copy complex number. This
     * measures the overhead of creation of new complex numbers including field access
     * to the real and imaginary parts.
     */
    @Benchmark
    public Complex[] baselineCopy(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), ComplexPerformance::copy);
    }

    // Unary operations that return a boolean

    @Benchmark
    public boolean[] isNaN(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::isNaN);
    }

    @Benchmark
    public boolean[] isInfinite(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::isInfinite);
    }

    @Benchmark
    public boolean[] isFinite(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::isFinite);
    }

    // Unary operations that return a double

    @Benchmark
    public double[] real(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::real);
    }

    @Benchmark
    public double[] imag(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::imag);
    }

    @Benchmark
    public double[] abs(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::abs);
    }

    @Benchmark
    public double[] arg(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::arg);
    }

    @Benchmark
    public double[] norm(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::norm);
    }

    /**
     * This test demonstrates that the {@link Math#hypot(double, double)} method
     * used in abs() is not as fast as using square root of the norm. Hypot is
     * within 1 ULP of the correct answer. The simple {@code Math.sqrt(x * x + y * y)}
     * may be much worse.
     */
    @Benchmark
    public double[] sqrtNorm(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> Math.sqrt(z.norm()));
    }

    @Benchmark
    public double[] absMathHypot(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> Math.hypot(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absFastMathHypot(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> FastMath.hypot(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absNoop(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotNoop(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absReplica(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotReplica(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absFmaSim(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotFmaSim(z.real(), z.imag()));
    }

    //@Benchmark
    public double[] absFmaElseNoAbs(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotFmaElseNoAbs(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absX2Y2(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotX2Y2(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absDekker(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotDekker(z.real(), z.imag()));
    }

    @Benchmark
    public double[] absFma(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), (ToDoubleFunction<Complex>) z -> hypotFma(z.real(), z.imag()));
    }

    private static double hypotNoop(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fffffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fffffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Second re-order in the rare event the upper 32-bits are the same.
        // This could be done in various locations including inside the function x2y2.
        // It must be done for sub-normals and is placed here for clarity that x2y2 assumes x >= y.
        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return b == Double.POSITIVE_INFINITY ? b : a;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(x2y2Noop(a, b)) * rescale;
    }

    /**
     * Return {@code x^2 + y^2} with high accuracy.
     *
     * <p>It is assumed that {@code 2^500 > x >= y > 2^-500}. Thus there will be no
     * overflow or underflow of the result.
     *
     * <p>Translated from "Freely Distributable Math Library" e_hypot.c
     * to minimize rounding errors. The code has been modified to use the split
     * method from Dekker. See the method {@link #hypotNoop(double, double)} for the
     * original license.
     *
     * <p>Note that fdlibm conditioned {@code x > y} using only the upper 20 bits of the
     * mantissa and the 11-bit exponent and thus x could be slightly smaller than y
     * for cases where the upper 31-bits of the unsigned 64-bit representation are identical.
     * There are cases for {@code 2x > y > x} where the result is different if the
     * arguments were reversed so the equality should be conditioned on the entire magnitude
     * of the two floating point values. This is especially relevant if the numbers are
     * sub-normal where the upper 20-bits of the mantissa are zero and the entire information
     * on the magnitude is in the lower 32-bits.
     *
     * @param x Value x.
     * @param y Value y
     * @return x^2 + y^2
     */
    private static double x2y2Noop(double x, double y) {
        // Below we use a Dekker split to create two 26-bit numbers avoiding use of
        // conversion to raw bits. The original fdlibm version used a split of
        // the upper 32-bits and lower 32-bits. This results in the upper part
        // being a 21-bit precision number (including assumed leading 1) and the
        // The different split can create 1 ULP differences
        // in the result of Math.sqrt(x2y2(x, y)) verses the JDK Math.hypot(x, y)
        // implementation of the fdlibm function. Replacing the Dekker split with
        // Double.longBitsToDouble(Double.doubleToRawLongBits(a) & 0xffffffff00000000L)
        // computes the same result as Math.hypot(x, y).
        // Testing with billions of numbers show that either split achieves the target of
        // 1 ULP from the result of a high precision computation
        // (e.g. see o.a.c.numbers.arrays.LinearCombination) and standard Math.sqrt().
        // Performance comparisons show the Dekker split is faster.
        //
        // This could be replaced by Math.fma(x, x, y * y).
        // Tests with FMA on JDK 9 show precision is not significantly changed from this
        // implementation. Using simulated fused multiply addition (FMA)
        // with split precision summation is slower. A full high precision computation
        // is slower again and is a poor trade-off to gain 1 ULP accuracy.
        final double w = x - y;
        if (w > y) {
            return x;
            //return x2y2b(x, y);
//            final double t1 = splitHigh(x);
//            final double t2 = x - t1;
//            return t1 * t1 - (y * (-y) - t2 * (x + t1));
        }
        // 2y > x > y
        return x * x + y * y;
//        final double t = x + x;
//        final double y1 = splitHigh(y);
//        final double y2 = y - y1;
//        final double t1 = splitHigh(t);
//        final double t2 = t - t1;
//        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));
    }

    private static double hypotReplica(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fffffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fffffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Second re-order in the rare event the upper 32-bits are the same.
        // This could be done in various locations including inside the function x2y2.
        // It must be done for sub-normals and is placed here for clarity that x2y2 assumes x >= y.
        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return b == Double.POSITIVE_INFINITY ? b : a;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(x2y2Replica(a, b)) * rescale;
    }

    /**
     * Return {@code x^2 + y^2} with high accuracy.
     *
     * <p>It is assumed that {@code 2^500 > x >= y > 2^-500}. Thus there will be no
     * overflow or underflow of the result.
     *
     * <p>Translated from "Freely Distributable Math Library" e_hypot.c
     * to minimize rounding errors. The code has been modified to use the split
     * method from Dekker. See the method {@link #hypotNoop(double, double)} for the
     * original license.
     *
     * <p>Note that fdlibm conditioned {@code x > y} using only the upper 20 bits of the
     * mantissa and the 11-bit exponent and thus x could be slightly smaller than y
     * for cases where the upper 31-bits of the unsigned 64-bit representation are identical.
     * There are cases for {@code 2x > y > x} where the result is different if the
     * arguments were reversed so the equality should be conditioned on the entire magnitude
     * of the two floating point values. This is especially relevant if the numbers are
     * sub-normal where the upper 20-bits of the mantissa are zero and the entire information
     * on the magnitude is in the lower 32-bits.
     *
     * @param x Value x.
     * @param y Value y
     * @return x^2 + y^2
     */
    private static double x2y2Replica(double x, double y) {
        // Below we use a Dekker split to create two 26-bit numbers avoiding use of
        // conversion to raw bits. The original fdlibm version used a split of
        // the upper 32-bits and lower 32-bits. This results in the upper part
        // being a 21-bit precision number (including assumed leading 1) and the
        // The different split can create 1 ULP differences
        // in the result of Math.sqrt(x2y2(x, y)) verses the JDK Math.hypot(x, y)
        // implementation of the fdlibm function. Replacing the Dekker split with
        // Double.longBitsToDouble(Double.doubleToRawLongBits(a) & 0xffffffff00000000L)
        // computes the same result as Math.hypot(x, y).
        // Testing with billions of numbers show that either split achieves the target of
        // 1 ULP from the result of a high precision computation
        // (e.g. see o.a.c.numbers.arrays.LinearCombination) and standard Math.sqrt().
        // Performance comparisons show the Dekker split is faster.
        //
        // This could be replaced by Math.fma(x, x, y * y).
        // Tests with FMA on JDK 9 show precision is not significantly changed from this
        // implementation. Using simulated fused multiply addition (FMA)
        // with split precision summation is slower. A full high precision computation
        // is slower again and is a poor trade-off to gain 1 ULP accuracy.
        final double w = x - y;
        if (w > y) {
            final double t1 = Double.longBitsToDouble(Double.doubleToLongBits(x) & 0xffff_ffff_0000_0000L);
            final double t2 = x - t1;
            return t1 * t1 - (y * (-y) - t2 * (x + t1));
        }
        // 2y > x > y
        final double t = x + x;
        final double y1 = Double.longBitsToDouble(Double.doubleToLongBits(y) & 0xffff_ffff_0000_0000L);
        final double y2 = y - y1;
        final double t1 = Double.longBitsToDouble(Double.doubleToLongBits(t) & 0xffff_ffff_0000_0000L);
        final double t2 = t - t1;
        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));
    }

    private static double hypotFmaSim(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fff_ffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fff_ffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        // No use of absolutes as the x^2 + y^2 sum removes the signs.
        // No requirement for a second order unless sub-normal.
        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return b == Double.POSITIVE_INFINITY ? Double.POSITIVE_INFINITY : a;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
            // Disallow parallel condition evaluation.
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        //} else {
        //    return Math.sqrt(x2y2b(a, b));
        }

//        final double yy = b * b;
//        final double xx = a * a;
//        final double xHigh = splitHigh(a);
//        final double xLow  = a - xHigh;
//        final double x2Low = squareLow(xLow, xHigh, xx);
//        // Two sum y^2 into the expansion of x^2
//        final double r = xx + yy;
//        return Math.sqrt(xx - r + yy + x2Low + r) * rescale;
        return Math.sqrt(x2y2Fma(a, b)) * rescale;
    }

    private static double hypotFmaElseNoAbs(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fff_ffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fff_ffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        // No use of absolutes as the x^2 + y^2 sum removes the signs.
        if (ha == hb && Math.abs(a) < Math.abs(b)) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return Math.abs(b) == Double.POSITIVE_INFINITY ? Double.POSITIVE_INFINITY : Math.abs(a);
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
            // Disallow parallel condition evaluation.
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return Math.abs(a);
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(x2y2Fma(a, b)) * rescale;
    }

    private static double hypotX2Y2(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fff_ffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fff_ffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        // No use of absolutes as the x^2 + y^2 sum removes the signs.

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return Math.abs(b) == Double.POSITIVE_INFINITY ? Double.POSITIVE_INFINITY : Math.abs(a);
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
            // Disallow parallel condition evaluation.
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return Math.abs(a);
                }
                if (ha == 0 && a == 0.0) {
                    return Math.abs(b);
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(x * x + y * y) * rescale;
    }

    /**
     * Return {@code x^2 + y^2} with high accuracy.
     *
     * <p>It is assumed that {@code 2^500 > x >= y > 2^-500}. Thus there will be no
     * overflow or underflow of the result.
     *
     * <p>Translated from "Freely Distributable Math Library" e_hypot.c
     * to minimize rounding errors. The code has been modified to use the split
     * method from Dekker. See the method {@link #hypotNoop(double, double)} for the
     * original license.
     *
     * <p>Note that fdlibm conditioned {@code x > y} using only the upper 20 bits of the
     * mantissa and the 11-bit exponent and thus x could be slightly smaller than y
     * for cases where the upper 31-bits of the unsigned 64-bit representation are identical.
     * There are cases for {@code 2x > y > x} where the result is different if the
     * arguments were reversed so the equality should be conditioned on the entire magnitude
     * of the two floating point values. This is especially relevant if the numbers are
     * sub-normal where the upper 20-bits of the mantissa are zero and the entire information
     * on the magnitude is in the lower 32-bits.
     *
     * @param x Value x.
     * @param y Value y
     * @return x^2 + y^2
     */
    private static double x2y2Fma(double x, double y) {
        // Below we use a Dekker split to create two 26-bit numbers avoiding use of
        // conversion to raw bits. The original fdlibm version used a split of
        // the upper 32-bits and lower 32-bits. This results in the upper part
        // being a 21-bit precision number (including assumed leading 1) and the
        // The different split can create 1 ULP differences
        // in the result of Math.sqrt(x2y2(x, y)) verses the JDK Math.hypot(x, y)
        // implementation of the fdlibm function. Replacing the Dekker split with
        // Double.longBitsToDouble(Double.doubleToRawLongBits(a) & 0xffffffff00000000L)
        // computes the same result as Math.hypot(x, y).
        // Testing with billions of numbers show that either split achieves the target of
        // 1 ULP from the result of a high precision computation
        // (e.g. see o.a.c.numbers.arrays.LinearCombination) and standard Math.sqrt().
        // Performance comparisons show the Dekker split is faster.
        //
        // This could be replaced by Math.fma(x, x, y * y).
        // Tests with FMA on JDK 9 show precision is not significantly changed from this
        // implementation. Using simulated fused multiply addition (FMA)
        // with split precision summation is slower. A full high precision computation
        // is slower again and is a poor trade-off to gain 1 ULP accuracy.
//        final double w = x - y;
//        if (w > y) {
//            final double t1 = splitHigh(x);
//            final double t2 = x - t1;
//            return t1 * t1 - (y * (-y) - t2 * (x + t1));
//        }
//        // 2y > x > y
//        final double t = x + x;
//        final double y1 = splitHigh(y);
//        final double y2 = y - y1;
//        final double t1 = splitHigh(t);
//        final double t2 = t - t1;
//        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));

        // Simulate a fused multiply add (FMA)
        //return Math.fma(x, x, y * y);
//        final double yy = y * y;
//        final double xx = x * x;
//        final double xHigh = splitHigh(x);
//        final double xLow  = x - xHigh;
//        final double x2Low = squareLow(xLow, xHigh, xx);
//        // Two sum y^2 into the expansion of x^2
//        double e2 = yy + x2Low;
//        final double e1 = sumLow(yy, x2Low, e2);
//        final double e3 = e2 + xx;
//        e2 = sumLow(e2, xx, e3);
//        return e1 + e2 + e3;

        // FMA does not require the use of abs since we have x^2 and y^2.
        // FMA 2 as a Dekker sum with y round-off of zero.
        final double yy = y * y;
        final double xx = x * x;
        final double xHigh = splitHigh(x);
        final double xLow  = x - xHigh;
        final double x2Low = squareLow(xLow, xHigh, xx);
        // Two sum y^2 into the expansion of x^2
        final double r = xx + yy;
        return xx - r + yy + x2Low + r;

//        // Do a Dekker summation
//        final double xx = x * x;
//        final double yy = y * y;
//        final double xHigh = splitHigh(x);
//        final double xLow = x - xHigh;
//        final double yHigh = splitHigh(y);
//        final double yLow = y - yHigh;
//        final double x2Low = squareLow(xLow, xHigh, xx);
//        final double y2Low = squareLow(yLow, yHigh, yy);
//
//        final double r = xx + yy;
//        return xx - r + yy + y2Low + x2Low + r;
    }

    private static double hypotDekker(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fffffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fffffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Second re-order in the rare event the upper 32-bits are the same.
        // This could be done in various locations including inside the function x2y2.
        // It must be done for sub-normals and is placed here for clarity that x2y2 assumes x >= y.
        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return b == Double.POSITIVE_INFINITY ? b : a;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(x2y2Dekker(a, b)) * rescale;
    }

    private static double x2y2Dekker(double x, double y) {
        // Below we use a Dekker split to create two 26-bit numbers avoiding use of
        // conversion to raw bits. The original fdlibm version used a split of
        // the upper 32-bits and lower 32-bits. This results in the upper part
        // being a 21-bit precision number (including assumed leading 1) and the
        // The different split can create 1 ULP differences
        // in the result of Math.sqrt(x2y2(x, y)) verses the JDK Math.hypot(x, y)
        // implementation of the fdlibm function. Replacing the Dekker split with
        // Double.longBitsToDouble(Double.doubleToRawLongBits(a) & 0xffffffff00000000L)
        // computes the same result as Math.hypot(x, y).
        // Testing with billions of numbers show that either split achieves the target of
        // 1 ULP from the result of a high precision computation
        // (e.g. see o.a.c.numbers.arrays.LinearCombination) and standard Math.sqrt().
        // Performance comparisons show the Dekker split is faster.
        //
        // This could be replaced by Math.fma(x, x, y * y).
        // Tests with FMA on JDK 9 show precision is not significantly changed from this
        // implementation. Using simulated fused multiply addition (FMA)
        // with split precision summation is slower. A full high precision computation
        // is slower again and is a poor trade-off to gain 1 ULP accuracy.
//        final double w = x - y;
//        if (w > y) {
//            final double t1 = splitHigh(x);
//            final double t2 = x - t1;
//            return t1 * t1 - (y * (-y) - t2 * (x + t1));
//        }
//        // 2y > x > y
//        final double t = x + x;
//        final double y1 = splitHigh(y);
//        final double y2 = y - y1;
//        final double t1 = splitHigh(t);
//        final double t2 = t - t1;
//        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));

        // Simulate a fused multiply add (FMA)
        //return Math.fma(x, x, y * y);
//        final double yy = y * y;
//        final double xx = x * x;
//        final double xHigh = splitHigh(x);
//        final double xLow  = x - xHigh;
//        final double x2Low = squareLow(xLow, xHigh, xx);
//        // Two sum y^2 into the expansion of x^2
//        double e2 = yy + x2Low;
//        final double e1 = sumLow(yy, x2Low, e2);
//        final double e3 = e2 + xx;
//        e2 = sumLow(e2, xx, e3);
//        return e1 + e2 + e3;

//        // FMA does not require the use of abs since we have x^2 and y^2.
//        // FMA 2 as a Dekker sum with y round-off of zero.
//        final double yy = y * y;
//        final double xx = x * x;
//        final double xHigh = splitHigh(x);
//        final double xLow  = x - xHigh;
//        final double x2Low = squareLow(xLow, xHigh, xx);
//        // Two sum y^2 into the expansion of x^2
//        final double r = xx + yy;
//        return xx - r + yy + x2Low + r;

        // Do a Dekker summation
        final double xx = x * x;
        final double yy = y * y;
        final double xHigh = splitHigh(x);
        final double xLow = x - xHigh;
        final double yHigh = splitHigh(y);
        final double yLow = y - yHigh;
        final double x2Low = squareLow(xLow, xHigh, xx);
        final double y2Low = squareLow(yLow, yHigh, yy);

        final double r = xx + yy;
        return xx - r + yy + y2Low + x2Low + r;
    }

    private static double hypotFma(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates the x^2+y^2 sum for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // to ensure the function is commutative.
        // The current method of doing the comparison only when the upper 32-bits match
        // is consistently fast in JMH performance tests with alternative implementations.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by directly writing to the exponent bits.
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation computes
        // the split dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases where the implementations compute
        // a different result. Only 1 ULP differences have been measured.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign bit.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fffffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fffffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_54) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Second re-order in the rare event the upper 32-bits are the same.
        // This could be done in various locations including inside the function x2y2.
        // It must be done for sub-normals and is placed here for clarity that x2y2 assumes x >= y.
        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return b == Double.POSITIVE_INFINITY ? b : a;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(Math.fma(a, a, b * b)) * rescale;
    }

    private static double splitHigh(double a) {
        final double c = MULTIPLIER * a;
        return c - (c - a);
    }

    private static double squareLow(double low, double high, double square) {
        return low * low - ((square - high * high) - 2 * low * high);
    }

    // Unary operations that return a complex number

    @Benchmark
    public Complex[] conj(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::conj);
    }

    @Benchmark
    public Complex[] negate(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::negate);
    }

    @Benchmark
    public Complex[] proj(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::proj);
    }

    @Benchmark
    public Complex[] cos(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::cos);
    }

    @Benchmark
    public Complex[] cosh(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::cosh);
    }

    @Benchmark
    public Complex[] exp(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::exp);
    }

    @Benchmark
    public Complex[] log(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::log);
    }

    @Benchmark
    public Complex[] log10(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::log10);
    }

    @Benchmark
    public Complex[] sin(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::sin);
    }

    @Benchmark
    public Complex[] sinh(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::sinh);
    }

    @Benchmark
    public Complex[] sqrt(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::sqrt);
    }

    @Benchmark
    public Complex[] tan(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::tan);
    }

    @Benchmark
    public Complex[] tanh(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::tanh);
    }

    @Benchmark
    public Complex[] acos(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::acos);
    }

    @Benchmark
    public Complex[] acosh(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::acosh);
    }

    @Benchmark
    public Complex[] asin(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::asin);
    }

    @Benchmark
    public Complex[] asinh(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::asinh);
    }

    @Benchmark
    public Complex[] atan(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::atan);
    }

    @Benchmark
    public Complex[] atanh(ComplexNumbers numbers) {
        return apply(numbers.getNumbers(), Complex::atanh);
    }

    // Binary operations on two complex numbers.

    @Benchmark
    public Complex[] pow(TwoComplexNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::pow);
    }

    @Benchmark
    public Complex[] multiply(TwoComplexNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::multiply);
    }

    @Benchmark
    public Complex[] divide(TwoComplexNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::divide);
    }

    @Benchmark
    public Complex[] add(TwoComplexNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::add);
    }

    @Benchmark
    public Complex[] subtract(TwoComplexNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::subtract);
    }

    // Binary operations on a complex and a real number.
    // These only benchmark methods on the real component as the
    // following are expected to be the same speed as the real-only operations
    // given the equivalent primitive operations:
    // - multiplyImaginary
    // - divideImaginary
    // - addImaginary
    // - subtractImaginary
    // - subtractFrom
    // - subtractFromImaginary

    @Benchmark
    public Complex[] powReal(ComplexAndRealNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::pow);
    }

    @Benchmark
    public Complex[] multiplyReal(ComplexAndRealNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::multiply);
    }

    @Benchmark
    public Complex[] divideReal(ComplexAndRealNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::divide);
    }

    @Benchmark
    public Complex[] addReal(ComplexAndRealNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::add);
    }

    @Benchmark
    public Complex[] subtractReal(ComplexAndRealNumbers numbers) {
        return apply(numbers.getNumbers(), numbers.getNumbers2(), Complex::subtract);
    }
}
