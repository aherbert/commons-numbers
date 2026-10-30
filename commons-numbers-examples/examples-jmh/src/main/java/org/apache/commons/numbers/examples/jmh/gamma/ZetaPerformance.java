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

package org.apache.commons.numbers.examples.jmh.gamma;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;
import org.apache.commons.numbers.gamma.HurwitzZeta;
import org.apache.commons.numbers.gamma.RiemannZeta;
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
 * Executes a benchmark to estimate the speed of zeta function operations.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgs = {"-server", "-Xms512M", "-Xmx512M"})
public class ZetaPerformance {
    /** Number of terms of the series summation S. */
    private static final int N = 9;
    /** Convergence epsilon for the sum of the tail function. This prevents summation
     * of terms that do not affect the final result. Equal to 2^-53. */
    private static final double EPS = 0x1.0p-53;
    /** Asymptotic threshold for large {@code a}. Used when {@code a + N} is not accurate.
     * Equal to 2^53 - N. */
    private static final double LARGE_A = (1L << 53) - N;
    /** Threshold used when {@code a < 1} where {@code a^-s} == zeta(s, a) with the smallest supported
     * {@code s = 1 + 2^-52}. Set using zeta(s, 1) * 2^54. No possible evaluation of
     * all remaining terms for any s could be added to {@code a^-s} when {@code a < 1}.
     *
     * <p>Approximately equal to 2^106 thus remaining terms (all less than 1) cannot be added in
     * double-double precision. */
    private static final double LARGE_T0 = 8.11296384146066920939820185959e+31;
    /** Downscale to use on the initial rising factorial tail term (always above 1).
     * Equal to Double.MIN_NORMAL: 2^-1022. */
    private static final double DOWNSCALE = 0x1.0p-1022;
    /** Upscale to use on the initial probability tail term (always below 1).
     * Multiplicative inverse of Double.MIN_NORMAL: 2^1022. */
    private static final double UPSCALE = 0x1.0p1022;

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}
     * in double-double (DD) precision.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * Provides 25 terms.
     */
    private static final double[] F = {
        0.08333333333333333, // (1 / 6) / 2!
        -0.001388888888888889, // (-1 / 30) / 4!
        3.306878306878307E-5, // (1 / 42) / 6!
        -8.267195767195768E-7, // (-1 / 30) / 8!
        2.08767569878681E-8, // (5 / 66) / 10!
        -5.284190138687493E-10, // (-691 / 2730) / 12!
        1.3382536530684679E-11, // (7 / 6) / 14!
        -3.3896802963225827E-13, // (-3617 / 510) / 16!
        8.586062056277845E-15, // (43867 / 798) / 18!
        -2.174868698558062E-16, // (-174611 / 330) / 20!
        5.5090028283602295E-18, // (854513 / 138) / 22!
        -1.3954464685812522E-19, // (-236364091 / 2730) / 24!
        3.534707039629467E-21, // (8553103 / 6) / 26!
        -8.953517427037546E-23, // (-23749461029 / 870) / 28!
        2.267952452337683E-24, // (8615841276005 / 14322) / 30!
        -5.744790668872202E-26, // (-7709321041217 / 510) / 32!
        1.455172475614865E-27, // (2577687858367 / 6) / 34!
        -3.6859949406653103E-29, // (-26315271553053477373 / 1919190) / 36!
        9.336734257095045E-31, // (2929993913841559 / 6) / 38!
        -2.36502241570063E-32, // (-261082718496449122051 / 13530) / 40!
        5.990671762482134E-34, // (1520097643918070802691 / 1806) / 42!
        -1.5174548844682903E-35, // (-27833269579301024235023 / 690) / 44!
        3.843758125454189E-37, // (596451111593912163277961 / 282) / 46!
        -9.736353072646691E-39, // (-5609403368997817686249127547 / 46410) / 48!
        2.466247044200681E-40, // (495057205241079648212477525 / 66) / 50!
    };

    /**
     * Create {@code n} double values in {@code [min, max)}.
     *
     * @param rng the source of randomness
     * @param n the number
     * @param minMax the range as "min,max"
     * @return the doubles
     */
    private static double[] doubles(SplittableRandom rng, int n, String minMax) {
        final int index = minMax.indexOf(',');
        final double min = Double.parseDouble(minMax.substring(0, index));
        final double max = Double.parseDouble(minMax.substring(index + 1));
        return rng.doubles(n, min, max).toArray();
    }

    /**
     * Contains two arrays of numbers.
     */
    @State(Scope.Benchmark)
    public static class HurwitzZetaData {
        /** The size of the data. */
        @Param({"1000"})
        private int size;
        /** The range for s. */
        @Param({"1,4"})
        private String rangeS;
        /** The range for a. */
        @Param({
            // When a in [0, 1] the term a^-s dominates and any further terms
            // do not significantly impact the result
            "0,1",
            // When a in [1,16] the term (a+N) is inexact and an extended
            // precision power function can be used
            "1,16",
            // Avoid the asymptotic approximation at ~2^53 = 9007199254740992
            // This uses 1e12 which is 3/4 of the limit
            "16,1000000000000",
        })
        private String rangeA;

        /** The s values. */
        private double[] s;
        /** The a values. */
        private double[] a;

        /**
         * Gets the size of the array of numbers.
         *
         * @return the size
         */
        public int getSize() {
            return size;
        }

        /**
         * Gets the s values.
         *
         * @return the numbers
         */
        public double[] getS() {
            return s;
        }

        /**
         * Gets the a values.
         *
         * @return the numbers
         */
        public double[] getA() {
            return a;
        }

        /**
         * Create the numbers.
         */
        @Setup(Level.Iteration)
        public void setup() {
            final SplittableRandom rng = new SplittableRandom();
            s = doubles(rng, size, rangeS);
            a = doubles(rng, size, rangeA);
        }
    }

    /**
     * Contains an array of numbers.
     */
    @State(Scope.Benchmark)
    public static class RiemannZetaData {
        /** The size of the data. */
        @Param({"1000"})
        private int size;
        /**
         * The range for s.
         *
         * <p>The Boost implementation has rational approximations for
         * ranges: 0-1; 1-2; 2-4; 4-7*; 7-15*; 15-36*
         * (* use an Math.exp call)
         *
         * <p>Above s=36 the function is effectively 1.0 + 2^-s in double precision.
         */
        @Param({
            "0,4",
            "4,36",
            "0,53",
            })
        private String rangeS;

        /** The s values. */
        private double[] s;

        /**
         * Gets the size of the array of numbers.
         *
         * @return the size
         */
        public int getSize() {
            return size;
        }

        /**
         * Gets the s values.
         *
         * @return the numbers
         */
        public double[] getS() {
            return s;
        }

        /**
         * Create the numbers.
         */
        @Setup(Level.Iteration)
        public void setup() {
            final SplittableRandom rng = new SplittableRandom();
            s = doubles(rng, size, rangeS);
        }
    }

    /**
     * Computes the value of \( \zeta(s, a) \) over the domain \( s \gt 1 \) and
     * \( a \ne 0, -1, -2, \ldots, -n \).
     *
     * <p>Negative \( a \) is not supported apart from special cases.
     *
     * <p>Special cases:
     * <ul>
     * <li>If the argument \( s \) is 1, then the result is positive infinity.</li>
     * <li>If the argument \( s \lt 1 \), then the result is nan.</li>
     * <li>If the argument \( s \) is positive infinity and \( a \lt 1 \) then the result is positive infinity.</li>
     * <li>If the argument \( s \) is positive infinity and \( a = 1 \) then the result is 1.</li>
     * <li>If the argument \( s \) is positive infinity and \( a \gt 1 \) then the result is 0.</li>
     * <li>If the argument \( a \le 0 \) and is an integer, then the result is positive infinity.</li>
     * <li>If the argument \( a \lt 0 \) and \( s \) is not an integer, then the result is nan.</li>
     * <li>If the argument \( a \) is negative infinity, then the result is nan.</li>
     * <li>If either argument is nan, then the result is nan.</li>
     * </ul>
     *
     * @param s Argument.
     * @param a Argument.
     * @return \( \zeta(s, a) \)
     */
    static double zeta(double s, double a) {
        if (Double.isNaN(s) || Double.isNaN(a) || s < 1 || a == Double.NEGATIVE_INFINITY) {
            return Double.NaN;
        }
        // s >= 1
        // a > -infinity
        // Check special cases
        if (s == Double.POSITIVE_INFINITY) {
            if (a == 1) {
                // Riemman zeta
                return 1;
            }
            return a < 1 ? Double.POSITIVE_INFINITY : 0;
        }
        if (s == 1) {
            return Double.POSITIVE_INFINITY;
        }
        if (a <= 0) {
            final double ca = Math.ceil(a);
            if (ca == a) {
                // The term 0^-s is infinity
                return Double.POSITIVE_INFINITY;
            }
            if (Math.floor(s) != s) {
                // pow(a, -s) is not defined for negative a when s is non-integer
                return Double.NaN;
            }
            throw new IllegalStateException("Unsupported");
        }
        return zetaImp(s, a);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)} for positive {@code a}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * <p>This method matches the implementation in {@link HurwitzZeta} and is modified
     * to <strong>remove</strong>:
     *
     * <ul>
     * <li>the extended precision power function
     * <li>summation in double-double precision
     * </ul>
     *
     * <p>The function has a max ULP error of twice the HurwitzZeta function at approximately
     * 4 versus 2.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @return zeta(s, a)
     */
    private static double zetaImp(double s, double a) {
        // First term
        final double t0 = Math.pow(a, -s);

        // Can overflow if 0 <= a < 1
        if (t0 >= LARGE_T0) {
            return t0;
        }

        // Asymptotic behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0 (no series sum S): the I term and the first term of T.
        if (a > LARGE_A) {
            return Math.pow(a, 1 - s) / (s - 1) + t0 * 0.5;
        }

        // Now any (a+n)^-s cannot overflow and the sum cannot overflow

        // Methods are used to increase precision
        // - Extended precision summation (low cost)
        // - Extended precision power function using round-off from (a+k) (selectively applied)
        // - Scaling for sub-normal tail terms (no cost)
        // - Adding series terms in magnitude order (no cost)

        double p = Math.pow(a + N, -s);

        // Initialise sum with the first tail term
        double sum = 0.5 * p;
        // S : k in [0, n-1]
        // Note k=0 already computed as t0
        for (int k = N - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum += Math.pow(a + k, -s);
        }

        // I
        final double ti = Math.pow(a + N, 1 - s) / (s - 1);

        // The dominant term may be t0 when a in [0, 1], or ti when a is large and s is small.
        // Adding in double-double precision handles either magnitude order.
        if (t0 > ti) {
            sum += ti;
            sum += t0;
        } else {
            sum += t0;
            sum += ti;
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s) using division by (a+N)^2.
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        //
        // If a or s are large then (a+n)^-(2k-1+s) will underflow.
        // We can downscale the rising factorial (always above 1)
        // and upscale p (always below 1) by the same amount.
        // Scaling could dynamically use the exponent of f and p. Here we use
        // fixed scaling which can be pre-applied to a possible sub-normal
        // power term. This allows a few terms of the tail to be computed
        // when the zeta result or terms are sub-normal and a significant
        // improvement in accuracy for sub-normal results.
        //
        // The loop will exit immediately if p = (a+n)^-(1-s) = 0:
        // Extend support for large s by computing sub-normal (a+n)^-(1+s)
        // in parts with fixed scaling:
        //   scale * (a+n)^-(1+s) = (a+n)^-(1+s)/2 * scale * (a+n)^-(1+s)/2
        // The limit on (1+s):
        //   9^-(1+s) = 2^-(1075+1022) : 1+s = log(2)/log(9) * (1075+1022) = 662
        //
        // Rising factorial term : (s)_{2k-1}
        // The maximum s is < 662.
        // (662)_{2k-1}; k=53 = Pochammer(662, 105) = 3.96e299
        // The rising factorial will not overflow for k <= 53 before (a+n)^-(1+s) is
        // zero. This is within the length of table F.
        // For larger s where p=0 the loop exits immediately.
        double f = s * DOWNSCALE;
        // 2k - 1
        double k2 = 1;

        // Use to divide by (a+n)^2
        final double apn = Math.pow(a + N, -2);

        // The tail requires the most accurate start-point for (a+n)^-(1+s).
        // Compute using the power function not (a+n)^-s / (a+N).
        // Sub-normal p reduces accuracy:
        // p is currently (a+n)^-s. Test if the desired term is sub-normal.
        if (p > (a + N) * Double.MIN_NORMAL) {
            p = UPSCALE * Math.pow(a + N, -s - 1);
        } else {
            // Compute in parts with scaling: x^y = x^(y/2) * x^(y/2)
            p = Math.pow(a + N, (-s - 1) * 0.5);
            p = p * UPSCALE * p;
        }

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Note: an extended precision sum here has no effect on the final result.
        final double stop = sum * EPS;
        double tsum = 0;
        int i;
        for (i = 0; i < F.length; i++) {
            final double t = f * p * F[i];
            tsum += t;
            if (Math.abs(t) <= stop) {
                break;
            }
            // Update (a+n)^-(2k-1+s)
            p *= apn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
        }
        return sum + tsum;
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param x Numbers.
     * @param y Numbers.
     * @param fun Function.
     * @param bh Data sink.
     */
    private static void apply(double[] x, double[] y, DoubleBinaryOperator fun, Blackhole bh) {
        for (int i = 0; i < x.length; i++) {
            bh.consume(fun.applyAsDouble(x[i], y[i]));
        }
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param x Numbers.
     * @param fun Function.
     * @param bh Data sink.
     */
    private static void apply(double[] x, DoubleUnaryOperator fun, Blackhole bh) {
        for (int i = 0; i < x.length; i++) {
            bh.consume(fun.applyAsDouble(x[i]));
        }
    }

    /**
     * Add function. This can be used to measure the JMH overhead of calling a function
     * on an array of numbers.
     *
     * @param x Number.
     * @param y Number.
     * @return the number
     */
    private static double add(double x, double y) {
        return x + y;
    }

    /**
     * Identity function. This can be used to measure the JMH overhead of calling a function
     * on an array of numbers.
     *
     * @param z Number.
     * @return the number
     */
    private static double identity(double z) {
        return z;
    }

    // Benchmark methods.
    // Benchmarks use function references to perform different operations on the numbers.

    /**
     * Baseline the JMH overhead for all the benchmarks that evaluate a function of
     * an array of a pair of numbers. All other methods are expected to be slower than this.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void baseline2(HurwitzZetaData data, Blackhole bh) {
        apply(data.getS(), data.getA(), ZetaPerformance::add, bh);
    }

    /**
     * Benchmark the HurwitzZeta function.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void hurwitzZeta(HurwitzZetaData data, Blackhole bh) {
        apply(data.getS(), data.getA(), HurwitzZeta::value, bh);
    }

    /**
     * Benchmark a standard precision implementation of the Hurwitz zeta function.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void hurwitzZetaStd(HurwitzZetaData data, Blackhole bh) {
        apply(data.getS(), data.getA(), ZetaPerformance::zeta, bh);
    }

    /**
     * Baseline the JMH overhead for all the benchmarks that evaluate a function of
     * an array of numbers. All other methods are expected to be slower than this.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void baseline1(RiemannZetaData data, Blackhole bh) {
        apply(data.getS(), ZetaPerformance::identity, bh);
    }

    /**
     * Benchmark the RiemannZeta function.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void riemannZeta(RiemannZetaData data, Blackhole bh) {
        apply(data.getS(), RiemannZeta::value, bh);
    }

    /**
     * Benchmark a standard precision implementation of the Riemann zeta function.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void riemannZetaStd(RiemannZetaData data, Blackhole bh) {
        apply(data.getS(), s -> ZetaPerformance.zeta(s, 1), bh);
    }

    // TODO:
    //
    // Add a unit test class that compares this imp vs the final
    // and checks they are close.
}
