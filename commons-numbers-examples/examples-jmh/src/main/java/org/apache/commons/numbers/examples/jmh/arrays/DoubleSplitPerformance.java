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

import java.util.concurrent.TimeUnit;
import java.util.function.DoubleUnaryOperator;

/**
 * Executes a benchmark to measure the speed of operations in the {@link LinearCombination} class.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgs = {"-server", "-Xms512M", "-Xmx512M"})
public class DoubleSplitPerformance {
    /**
     * The multiplier used to split the double value into high and low parts. From
     * Dekker (1971): "The constant should be chosen equal to 2^(p - p/2) + 1,
     * where p is the number of binary digits in the mantissa". Here p is 53
     * and the multiplier is {@code 2^27 + 1}.
     */
    private static final double MULTIPLIER = 1.34217729E8;

    /** The upper limit above which a number may overflow during the split into a high part.
     * Assuming the multiplier is above 2^27 and the maximum exponent is 1023 then a safe
     * limit is a value with an exponent of (1023 - 27) = 2^996. */
    private static final double SAFE_UPPER = 0x1.0p996;

    /** The scale to use when down-scaling during a split into a high part.
     * This must be smaller than the inverse of the multiplier and a power of 2 for exact scaling. */
    private static final double DOWN_SCALE = 0x1.0p-30;

    /** The scale to use when re-scaling during a split into a high part.
     * This is the inverse of {@link #DOWN_SCALE}. */
    private static final double UP_SCALE = 0x1.0p30;

    /** The mask to zero the lower 27-bits of a long . */
    private static final long ZERO_LOWER_27_BITS = 0xffff_ffff_f800_0000L;

    /**
     * The factors to multiply.
     */
    @State(Scope.Benchmark)
    public static class Factors {
        /** The mask for the sign bit and the mantissa. */
        private static final long SIGN_MATISSA_MASK = 0x800f_ffff_ffff_ffffL;
        /** The exponent for small numbers. */
        private static final long EXP_SMALL = Double.doubleToRawLongBits(1.0);
        /** The exponent for big numbers. */
        private static final long EXP_BIG = Double.doubleToRawLongBits(SAFE_UPPER);

        /**
         * The number of factors.
         */
        @Param({"10000"})
        private int size;

        /**
         * The fraction of small factors.
         *
         * <p>Note: The split numbers are used in multiplications.
         * It is unlikely that many numbers will be larger than the upper limit.
         * These numbers are edge cases that would cause overflow in multiplications if
         * the other number is anywhere close to the same magnitude.
         */
        @Param({"1", "0.999", "0.99", "0.9"})
        private double small;

        /** Factors a. */
        private double[] a;

        /**
         * Gets the a factors.
         *
         * @return Factors.
         */
        public double[] getFactors() {
            return a;
        }

        /**
         * Create the factors.
         */
        @Setup
        public void setup() {
            final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_1024_PP);
            a = new double[size];
            for (int i = 0; i < size; i++) {
                long bits = rng.nextLong() & SIGN_MATISSA_MASK;
                // The exponent will either be small or big
                if (rng.nextDouble() < small) {
                    bits |= EXP_SMALL;
                } else {
                    bits |= EXP_BIG;
                }
                a[i] = Double.longBitsToDouble(bits);
            }
        }
    }

    /**
     * Compute an operation on the initial double of each factor.
     *
     * @param factors Factors.
     * @param bh Data sink.
     * @param fun Scalar product function.
     */
    private static void apply(Factors factors, Blackhole bh, DoubleUnaryOperator fun) {
        final double[] a = factors.getFactors();
        for (int i = 0; i < a.length; i++) {
            bh.consume(fun.applyAsDouble(a[i]));
        }
    }

    /**
     * Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) creates
     * a big value from which to derive the two split parts.
     * <pre>
     * c = (2^s + 1) * a
     * a_big = c - a
     * a_hi = c - a_big
     * a_lo = a - a_hi
     * a = a_hi + a_lo
     * </pre>
     *
     * <p>The multiplicand allows a p-bit value to be split into
     * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
     * Combined they have (p􏰔-1) bits of significand but the sign bit of {@code a_lo}
     * contains a bit of information. The constant is chosen so that s is ceil(p/2) where
     * the precision p for a double is 53-bits (1-bit of the mantissa is assumed to be
     * 1 for a non sub-normal number) and s is 27.
     *
     * @param value Value.
     * @return the high part of the value.
     */
    private static double splitDekker(double value) {
        // Avoid overflow
        if (value >= SAFE_UPPER || value <= -SAFE_UPPER) {
            // Do scaling.
            final double x = value * DOWN_SCALE;
            final double c = MULTIPLIER * x;
            final double hi = (c - (c - x)) * UP_SCALE;
            if (Double.isInfinite(hi)) {
                // Number is too large.
                // This occurs if value is infinite or close to Double.MAX_VALUE.
                // Note that Dekker's split creates an approximating 26-bit number which may
                // have an exponent 1 greater than the input value. This will overflow if the
                // exponent is already +1023. Revert to the raw upper 26 bits of the 53-bit
                // mantissa (including the assumed leading 1 bit). This conversion will result in
                // the low part being a 27-bit significand and the potential loss of bits during
                // addition and multiplication. (Contrast to the Dekker split which creates two
                // 26-bit numbers with a bit of information moved to the sign of low.)
                // The conversion will maintain Infinite in the high part where the resulting
                // low part a_lo = a - a_hi = inf - inf = NaN.
                return Double.longBitsToDouble(Double.doubleToRawLongBits(value) & ZERO_LOWER_27_BITS);
            }
            return hi;
        }
        // normal conversion
        final double c = MULTIPLIER * value;
        return c - (c - value);
    }

    /**
     * Implement Dekker's method to split a value into two parts.
     * This is as per {@link #splitDekker(double)} but uses a {@link Math#abs(double)} in the
     * condition statement.
     *
     * @param value Value.
     * @return the high part of the value.
     */
    private static double splitDekkerAbs(double value) {
        if (Math.abs(value) >= SAFE_UPPER) {
            final double x = value * DOWN_SCALE;
            final double c = MULTIPLIER * x;
            final double hi = (c - (c - x)) * UP_SCALE;
            if (Double.isInfinite(hi)) {
                return Double.longBitsToDouble(Double.doubleToRawLongBits(value) & ZERO_LOWER_27_BITS);
            }
            return hi;
        }
        final double c = MULTIPLIER * value;
        return c - (c - value);
    }

    /**
     * Implement Dekker's method to split a value into two parts.
     * This is as per {@link #splitDekker(double)} but has no overflow protection.
     *
     * @param value Value.
     * @return the high part of the value.
     */
    private static double splitDekkerRaw(double value) {
        final double c = MULTIPLIER * value;
        return c - (c - value);
    }

    // Benchmark methods.
    //
    // The methods are partially documented as the names are self-documenting.
    // CHECKSTYLE: stop JavadocMethod
    // CHECKSTYLE: stop DesignForExtension

    /**
     * Baseline using a no-operation on the factors.
     *
     * @param factors Factors.
     * @param bh Data sink.
     */
    @Benchmark
    public void baseline(Factors factors, Blackhole bh) {
        apply(factors, bh, a -> a);
    }

    /**
     * Dekker's split with overflow protection for large magnitudes; uses two signed comparisons
     * for the magnitude check.
     *
     * @param factors Factors.
     * @param bh Data sink.
     */
    @Benchmark
    public void dekker(Factors factors, Blackhole bh) {
        apply(factors, bh, DoubleSplitPerformance::splitDekker);
    }

    /**
     * Dekker's split with overflow protection for large magnitudes; uses a single unsigned
     * comparison (with an absolute value) for the magnitude check.
     *
     * @param factors Factors.
     * @param bh Data sink.
     */
    @Benchmark
    public void dekkerAbs(Factors factors, Blackhole bh) {
        apply(factors, bh, DoubleSplitPerformance::splitDekkerAbs);
    }

    /**
     * Dekker's split without overflow protection.
     *
     * @param factors Factors.
     * @param bh Data sink.
     */
    @Benchmark
    public void dekkerRaw(Factors factors, Blackhole bh) {
        apply(factors, bh, DoubleSplitPerformance::splitDekkerRaw);
    }

    /**
     * Split using the upper and lower raw bits from the double.
     *
     * <p>Note: This method will not work for very small sub-normal numbers
     * ({@code <= 27} bits) as the high part will be zero and the low part will
     * have all the information. Methods that assume {@code hi > lo} will have
     * undefined behaviour.
     *
     * @param factors Factors.
     * @param bh Data sink.
     */
    @Benchmark
    public void rawbits(Factors factors, Blackhole bh) {
        apply(factors, bh, a -> Double.longBitsToDouble(Double.doubleToRawLongBits(a) & ZERO_LOWER_27_BITS));
    }
}
