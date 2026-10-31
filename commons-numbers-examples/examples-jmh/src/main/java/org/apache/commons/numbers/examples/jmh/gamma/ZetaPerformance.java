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

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.DoubleStream;
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
    /** Maximum {@code s} value where {@code zeta(s, 1.5) > 0}. Used in the reflection
     * identity when {@code s} is odd and {@code a} is negative and half-integer. */
    private static final int MAX_REFLECTION_S = 1837;
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
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * This is the same as table {@link #F} in 36-digits of precision.
     * Provides 50 terms.
     */
    private static final BigDecimal[] FBD = {
        new BigDecimal("0.0833333333333333333333333333333333333"),
        new BigDecimal("-0.00138888888888888888888888888888888889"),
        new BigDecimal("0.0000330687830687830687830687830687830688"),
        new BigDecimal("-8.26719576719576719576719576719576720E-7"),
        new BigDecimal("2.08767569878680989792100903212014323E-8"),
        new BigDecimal("-5.28419013868749318484768220217955668E-10"),
        new BigDecimal("1.33825365306846788328269809751291233E-11"),
        new BigDecimal("-3.38968029632258286683019539124944250E-13"),
        new BigDecimal("8.58606205627784456413590545042562713E-15"),
        new BigDecimal("-2.17486869855806187304151642386591790E-16"),
        new BigDecimal("5.50900282836022951520265260890225488E-18"),
        new BigDecimal("-1.39544646858125233407076862640635498E-19"),
        new BigDecimal("3.53470703962946747169322997780379921E-21"),
        new BigDecimal("-8.95351742703754685040261131811274105E-23"),
        new BigDecimal("2.26795245233768306031095073886816606E-24"),
        new BigDecimal("-5.74479066887220244526388198760701840E-26"),
        new BigDecimal("1.45517247561486490186626486727132934E-27"),
        new BigDecimal("-3.68599494066531017818178247990866037E-29"),
        new BigDecimal("9.33673425709504467203255515278562330E-31"),
        new BigDecimal("-2.36502241570062993455963519636983824E-32"),
        new BigDecimal("5.99067176248213430465991239681965783E-34"),
        new BigDecimal("-1.51745488446829026171081313586471893E-35"),
        new BigDecimal("3.84375812545418823222944529099023211E-37"),
        new BigDecimal("-9.73635307264669103526762127925045418E-39"),
        new BigDecimal("2.46624704420068095710640028028884289E-40"),
        new BigDecimal("-6.24707674182074369314875679472336869E-42"),
        new BigDecimal("1.58240302446449142975108170682876394E-43"),
        new BigDecimal("-4.00827368594893596853001219052198266E-45"),
        new BigDecimal("1.01530758555695563116307139453787623E-46"),
        new BigDecimal("-2.57180415824187174992481940976445489E-48"),
        new BigDecimal("6.51445603523381493155843485864185802E-50"),
        new BigDecimal("-1.65013099068965245550609878047932301E-51"),
        new BigDecimal("4.17983062853947589485018723470940703E-53"),
        new BigDecimal("-1.05876346677029087702704202427911729E-54"),
        new BigDecimal("2.68187919126077066614098485884151034E-56"),
        new BigDecimal("-6.79327935110742120952718029953389461E-58"),
        new BigDecimal("1.72075776166814049053634994075823066E-59"),
        new BigDecimal("-4.35873032934889384340019984977316111E-61"),
        new BigDecimal("1.10407929036846667508383959764442732E-62"),
        new BigDecimal("-2.79666551337813450720479375311862655E-64"),
        new BigDecimal("7.08403650167947019850938842238034933E-66"),
        new BigDecimal("-1.79440740828922406660525730933675599E-67"),
        new BigDecimal("4.54528706361109610708507912463944819E-69"),
        new BigDecimal("-1.15133466319820518127300290078624338E-70"),
        new BigDecimal("2.91636477109236135470336898005297074E-72"),
        new BigDecimal("-7.38723826349733756257337539470973670E-74"),
        new BigDecimal("1.87120931176379530622526187140807466E-75"),
        new BigDecimal("-4.73982855776179940549956344121088708E-77"),
        new BigDecimal("1.20061259933545065198171002211412912E-78"),
        new BigDecimal("-3.04118724151429238303712207483889553E-80"),
    };

    /**
     * Create {@code n} double values in {@code [min, max)}. Optionally the values can be
     * integers in {@code [min, max]}, and/or returned as negative.
     *
     * @param rng the source of randomness
     * @param n the number
     * @param minMax the range as "min,max"
     * @param integers the flag to return integer values
     * @param negative the flag to return negative values
     * @return the doubles
     */
    private static double[] doubles(SplittableRandom rng, int n, String minMax,
        boolean integers, boolean negative) {
        final int index = minMax.indexOf(',');
        final double min = Double.parseDouble(minMax.substring(0, index));
        final double max = Double.parseDouble(minMax.substring(index + 1));
        DoubleStream s = integers ?
            rng.ints(n, (int) min, (int) max + 1).asDoubleStream() :
            rng.doubles(n, min, max);
        if (negative) {
            s = s.map(x -> -x);
        }
        return s.toArray();
    }

    /**
     * Contains arrays of numbers for zeta(s, a).
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
            s = doubles(rng, size, rangeS, false, false);
            a = doubles(rng, size, rangeA, false, false);
        }
    }

    /**
     * Contains arrays of numbers for zeta(s, a) with negative {@code a}.
     */
    @State(Scope.Benchmark)
    public static class NegativeHurwitzZetaData {
        /** The size of the data. */
        @Param({"1000"})
        private int size;
        /** The range for s (upper bound inclusive). */
        @Param({"2,5", "6,35"})
        private String rangeS;
        /** The range for a. */
        @Param({
            "0,10",
            "10,1000",
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
            // Require integer s and negative a
            s = doubles(rng, size, rangeS, true, false);
            a = doubles(rng, size, rangeA, false, true);
        }
    }

    /**
     * Contains arrays of numbers for zeta(s) or zeta(s, 1).
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
         * (* use an Math.exp call).
         *
         * <p>Above s=36 the function is effectively 1.0 + 2^-s in double precision;
         * above s=53 the function is ~1.0.
         */
        @Param({
            "0,4",
            "4,36",
            "36,53",
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
            s = doubles(rng, size, rangeS, false, false);
        }
    }

    /** Context for the BigDecimal zeta implementation. */
    private static final class Context {
        /** Context for double precision. */
        static final Context BD_DOUBLE = new Context(9, 15, new MathContext(18));
        /** Context for double-double precision. */
        static final Context BD_DOUBLE_DOUBLE = new Context(17, 30, new MathContext(36));

        /** N. */
        private final int n;
        /** M. */
        private final int m;
        /** Math context for extended precision evaluations. */
        private final MathContext mc;

        /**
         * Create an instance.
         *
         * @param n the number of terms for the series sum
         * @param m the number of terms for the tail sum
         * @param mc the math context for BigDecimal evaluations
         */
        private Context(int n, int m, MathContext mc) {
            this.n = n;
            this.m = m;
            this.mc = mc;
        }

        /**
         * Gets N.
         * @return n
         */
        int getN() {
            return n;
        }

        /**
         * Gets M.
         * @return m
         */
        int getM() {
            return m;
        }

        /**
         * Gets the math context for extended precision evaluations.
         * @return the math context
         */
        MathContext getMathContext() {
            return mc;
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
            return zetaNegative(s, a, ca);
        }
        return zetaImp(s, a);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)} for positive {@code a}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * <p>This method matches the implementation in {@link HurwitzZeta} and is modified for
     * the positive and negative {@code a} cases:
     *
     * <ul>
     * <li>{@code a > 0}: remove the extended precision power function
     * <li>{@code a > 0}: remove summation in double-double precision
     * <li>{@code a < 0}: use BigDecimal in the evaluation instead of double-double
     * </ul>
     *
     * <p>The function has a max ULP error of twice the HurwitzZeta function at approximately
     * 4 versus 2 when {@code a} is positive, and an increased accuracy when {@code a} is negative.
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
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}
     * when {@code a} is negative. Uses {@link BigDecimal} arithmetic.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be negative, and {@code s} integer.
     *
     * @param ds Argument {@code s > 1} and integer, specified as a double
     * @param a Argument {@code a < 0}
     * @param ca Ceil(a)
     * @return zeta(s, a)
     */
    private static double zetaNegative(double ds, double a, double ca) {
        // a < 0 (non-integer) and s is a positive integer.
        // If s is odd then the negative series sum will be negative
        // and the addition of zeta(s, x > 0) has cancellation.
        // This is largest when a is close to half-integer.

        // Check case of total cancellation.
        final boolean odd = ((long) -ds & 0x1) == 1;
        if (ca == a) {
            // The term 0^-s is infinity
            return Double.POSITIVE_INFINITY;
        }
        final double x = a - ca;
        // Intentional float comparison
        if (odd && x == -0.5) {
            // Compute zeta(s, 1 - a)
            if (ds <= MAX_REFLECTION_S) {
                // Use extended precision but evaluated with precision for a double result
                return zeta((int) ds, BigDecimal.ONE.subtract(new BigDecimal(a)), null,
                    Context.BD_DOUBLE).doubleValue();
            }
            // All possible zeta(s, 1-a >= 1.5) underflow
            return 0.0;
        }

        // Handle cancellation as x -> 0.5
        // Note: 0.5^-1024 overflows.
        // Limit of [nextDown(0.5)^-s - nextUp(0.5)^-s] may have terms above 2^1024.
        // 0.5 +/- 2^-54 (requires extended precision as ulp(0.5) is 2^-53)
        // The largest odd s where the difference is finite:
        // var mc = MathContext.DECIMAL128
        // var a = new BigDecimal(Math.nextDown(0.5))
        // var b = BigDecimal.ONE.subtract(a)
        // var s = -1025
        // while (Double.isFinite(a.pow(s, mc).subtract(b.pow(s, mc), mc).doubleValue())) { s -= 2; }
        // s = -1067 : diff = 3.75e308
        if (ds >= 1067) {
            // Use the dominant term using closest to zero
            return odd && x > -0.5 ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        }

        // Safe to cast
        final int s = (int) ds;

        // Odd computation requires twice the precision of a double (17 digits).
        // Even requires some extra.
        final Context c = odd ? Context.BD_DOUBLE_DOUBLE : Context.BD_DOUBLE;

        // Compute the two terms either side of zero:
        // -1 < xn < 0 < xn + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        final BigDecimal xn = new BigDecimal(x);
        final BigDecimal xp = BigDecimal.ONE.add(xn);
        final MathContext mc = c.getMathContext();
        final BigDecimal pn = xn.pow(-s, mc);
        final BigDecimal pp = xp.pow(-s, mc);

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2.
        // Exit early if remaining terms cannot be added.
        final double d = pn.add(pp, mc).doubleValue();
        if (Math.abs(d) > 0x1p106) {
            // Adding to either side will not change a double result
            return d;
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) + (-1)^s * [ zeta(s, -x) - zeta(s, 1 - a) ]
        // z + za + zb

        BigDecimal z = zeta(s, xp, pp, c);
        BigDecimal za;
        BigDecimal zb;
        // A single call to zeta uses many pow operations;
        // use a direct sum when zeta will use more.
        if (x - a > 2 * (c.getN() + 2)) {
            // Note: The difference (za - zb) incurs cancellation.
            // This should not be an issue as x in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 1)  = 1.6449
            // zeta(2, 2)  = 0.6449
            // zeta(2, 30) = 0.0339
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            za = zeta(s, xn.negate(), pn.abs(), c);
            zb = zeta(s, BigDecimal.ONE.subtract(new BigDecimal(a)), null, c);
            // Both terms are positive. Correct the sign for final addition.
            if (odd) {
                za = za.negate();
            } else {
                zb = zb.negate();
            }
        } else {
            // Sum terms in ascending order of magnitude.
            // Use a double to track the iterations, and mirror with a BigDecimal.
            za = pn;
            zb = BigDecimal.ZERO;
            BigDecimal ba = new BigDecimal(a);
            for (double aa = a; aa < x; aa += 1.0) {
                zb = zb.add(ba.pow(-s, mc), mc);
                ba = ba.add(BigDecimal.ONE);
            }
        }

        // Sum in magnitude order. Use the scale for a fast comparison of magnitude
        // as all results have the same precision: base 10 exponent = precision - scale - 1
        // Smaller scale is a bigger value.
        if (za.scale() < z.scale()) {
            final BigDecimal tmp = za;
            za = z;
            z = tmp;
        }
        // za < z
        if (zb.scale() < z.scale()) {
            final BigDecimal tmp = zb;
            zb = z;
            z = tmp;
        }
        // za,zb < z
        return za.add(zb, mc).add(z, mc).doubleValue();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @param a0 {@code a^-s} (can be null)
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    private static BigDecimal zeta(int s, BigDecimal a, BigDecimal a0, Context c) {
        final int n = c.getN();
        final MathContext mc = c.getMathContext();
        BigDecimal apn = a.add(BigDecimal.valueOf(n));
        BigDecimal p = apn.pow(-s, mc);

        // Initialise sum with the first tail term
        BigDecimal sum = p.multiply(new BigDecimal(0.5));
        // S : k in [0, n-1]
        for (int k = n - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum = sum.add(a.add(BigDecimal.valueOf(k)).pow(-s, mc), mc);
        }
        // First term may be provided
        final BigDecimal t0 = a0 == null ? a.pow(-s, mc) : a0;

        // I : (a+p)^(1-s) / (s-1)
        final BigDecimal ti = apn.pow(1 - s, mc).divide(BigDecimal.valueOf(s - 1), mc);

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (t0.compareTo(ti) > 0) {
            sum = sum.add(ti, mc).add(t0, mc);
        } else {
            sum = sum.add(t0, mc).add(ti, mc);
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.

        // *** This must work for all s where (a+n)^-s is not an underflow ***
        // The function testTailConvergence examines the parameter space in DD arithmetic
        // and verifies convergence is possible for the provided n.

        // Rising factorial term : (s)_{2k-1}
        BigDecimal f = BigDecimal.valueOf(s);

        // Alternative implementations for (a+n)^-(2k-1+s).

        // Initialise (a+n)^-(2k-1+s) to (a+n)^-(1+s)
        p = apn.pow(-(1 + s), mc);
        // Divide by (a+n)^2 using multiplication
        apn = apn.pow(-2, mc);

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        BigDecimal tsum = BigDecimal.ZERO;
        final int stop = sum.scale() + mc.getPrecision();
        int i;
        for (i = 0; i < c.getM(); i++) {
            final BigDecimal t = f.multiply(p, mc).multiply(FBD[i], mc);
            tsum = tsum.add(t, mc);
            // Bigger scale is a smaller value
            if (t.scale() >= stop) {
                break;
            }
            // p = (a+n)^-(2k-1+s)
            p = p.multiply(apn, mc);
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            f = f.multiply(BigDecimal.valueOf((s + (2L * i) + 1) * (s + (2L * i) + 2)), mc);
        }
        return sum.add(tsum, mc);
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
     * Benchmark the HurwitzZeta function with negative a.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void negativeHurwitzZeta(NegativeHurwitzZetaData data, Blackhole bh) {
        apply(data.getS(), data.getA(), HurwitzZeta::value, bh);
    }

    /**
     * Benchmark a BigDecimal implementation of the Hurwitz zeta function with negative a.
     *
     * @param data Test data.
     * @param bh Data sink.
     */
    @Benchmark
    public void negativeHurwitzZetaBD(NegativeHurwitzZetaData data, Blackhole bh) {
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
