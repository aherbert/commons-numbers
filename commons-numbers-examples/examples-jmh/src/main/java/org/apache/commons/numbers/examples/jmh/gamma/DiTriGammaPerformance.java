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
import java.util.function.DoubleUnaryOperator;
import org.apache.commons.numbers.gamma.Digamma;
import org.apache.commons.numbers.gamma.Trigamma;
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

/**
 * Executes a benchmark to estimate the speed of digamma and trigamma functions.
 *
 * @since 1.4
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@State(Scope.Benchmark)
@Fork(value = 1, jvmArgs = {"-server", "-Xms512M", "-Xmx512M"})
public class DiTriGammaPerformance {
    /**
     * Contains an array of numbers.
     */
    @State(Scope.Benchmark)
    public static class NumberData {
        /** The size of the data. */
        @Param({"10000"})
        private int size;
        /** The seed of the data. */
        @Param({"12345"})
        private long seed;
        /**
         * The range of the data.
         *
         * <p>Note Digamma asymptotic value uses log(x). This is at x=10 for Boost digamma
         * and x above 49 for {@link Digamma13}.
         *
         * <p>The {@link Digamma13} implementation uses the reflection formula at
         * {@code x < 0}. Boost digamma uses it at {@code x < -1}. Testing for large
         * negative values should reflect into the asymptotic tail where both use a call
         * to the log function. Evaluation of [-1, 0] shows the cost of the reflection for
         * digamma when it is not yet required.
         *
         * <p>Both {@link Trigamma13} and Boost trigamma use reflection when
         * {@code x < 0}.
         */
        @Param({"-100:-50", "-1:0",
                "0:1", "1:2", "2:4", "4:8", "8:16", "16:32", "32:64", "64:128"})
        private String range;

        /** The numbers. */
        private double[] numbers;

        /**
         * Gets the size of the array of numbers.
         *
         * @return the size
         */
        public int getSize() {
            return size;
        }

        /**
         * Gets the numbers.
         *
         * @return the numbers
         */
        public double[] getNumbers() {
            return numbers;
        }

        /**
         * Create the numbers.
         */
        @Setup
        public void setup() {
            final int index = range.indexOf(':');
            if (index < 0) {
                throw new IllegalStateException("Invalid range = " + range);
            }
            final double min = Double.parseDouble(range.substring(0, index));
            final double max = Double.parseDouble(range.substring(index + 1));
            numbers = new SplittableRandom(seed).doubles(size, min, max).toArray();
        }
    }

    /**
     * Contains a unary operator method.
     */
    @State(Scope.Benchmark)
    public static class MethodData {
        /** The method. */
        @Param({"baseline",
                "digamma", "digammaShift", "digammaLog", "digamma1.3",
                "trigamma", "trigamma1.3", "digamma1.3Log"})
        private String method;

        /** The operator. */
        private DoubleUnaryOperator op;

        /**
         * Gets the operator.
         *
         * @return the operator
         */
        public DoubleUnaryOperator getOperator() {
            return op;
        }

        /**
         * Create the operator.
         */
        @Setup
        public void setup() {
            if ("baseline".equals(method)) {
                op = DoubleUnaryOperator.identity();
            } else if ("digamma".equals(method)) {
                op = Digamma::value;
            } else if ("digammaShift".equals(method)) {
                // Copy implementation that shifts x to [1, 2] but does not
                // perform the sum of the 1 / x terms.
                op = BoostDigamma::digamma;
            } else if ("digammaLog".equals(method)) {
                // Copy implementation that computes only the log term for x > 10.
                // For x < 10 it shows the cost without the shift or evaluation of x in [1, 2].
                op = BoostDigamma::log;
            } else if ("digamma1.3".equals(method)) {
                op = Digamma13::value;
            } else if ("digamma1.3Log".equals(method)) {
                op = Digamma13::log;
            } else if ("trigamma".equals(method)) {
                op = Trigamma::value;
            } else if ("trigamma1.3".equals(method)) {
                op = Trigamma13::value;
            // For testing if negative x is an integer (pole)
            } else if ("rint".equals(method)) {
                op = Math::rint;
            } else if ("floor".equals(method)) {
                op = Math::floor;
            } else {
                throw new IllegalStateException("Unknown method = " + method);
            }
        }
    }

    /**
     * Copy of Commons Numbers Digamma from version 1.3.
     * This includes fixes applied before the 1.4 release to the validated domain (NUMBERS-219).
     */
    private static final class Digamma13 {
        /** <a href="https://en.wikipedia.org/wiki/Euler-Mascheroni_constant">Euler-Mascheroni constant</a>. */
        private static final double GAMMA = 0.577215664901532860606512090082;

        /** C limit. */
        private static final double C_LIMIT = 49;
        /** S limit. */
        private static final double S_LIMIT = 1e-5;
        /** Fraction. */
        private static final double F_M1_12 = -1d / 12;
        /** Fraction. */
        private static final double F_1_120 = 1d / 120;
        /** Fraction. */
        private static final double F_M1_252 = -1d / 252;

        /** Private constructor. */
        private Digamma13() {
            // intentionally empty.
        }

        /**
         * Computes the digamma function.
         *
         * <p>This is an independently written implementation of the algorithm described in
         * <a href="https://www.uv.es/~bernardo/1976AppStatist.pdf">Jose Bernardo,
         * Algorithm AS 103: Psi (Digamma) Function, Applied Statistics, 1976</a>.
         * A <a href="https://en.wikipedia.org/wiki/Digamma_function#Reflection_formula">
         * reflection formula</a> is incorporated to improve performance on negative values.
         *
         * <p>Some of the constants have been changed to increase accuracy at the moderate
         * expense of run-time.  The result should be accurate to within {@code 1e-8}.
         * relative tolerance for {@code 0 < x < 1e-5}  and within {@code 1e-8} absolute
         * tolerance otherwise.
         *
         * @param x Argument.
         * @return digamma(x) to within {@code 1e-8} relative or absolute error whichever
         * is larger.
         */
        static double value(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // large negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? x : Double.NaN;
            }

            double digamma = 0;
            if (x < 0) {
                // Use reflection formula to fall back into positive values.
                digamma = -Math.PI / Math.tan(Math.PI * x);
                x = 1 - x;
            }

            if (x <= S_LIMIT) {
                // Use method 5 from Bernardo AS103, accurate to O(x).
                return digamma - GAMMA - 1 / x;
            }

            while (x < C_LIMIT) {
                digamma -= 1 / x;
                x += 1;
            }

            // Use method 4, accurate to O(1/x^8)
            final double inv = 1 / (x * x);
            //            1       1        1         1
            // log(x) -  --- - ------ + ------- - -------
            //           2 x   12 x^2   120 x^4   252 x^6
            digamma += Math.log(x) - 0.5 / x + inv * (F_M1_12 + inv * (F_1_120 + F_M1_252 * inv));

            return digamma;
        }

        /**
         * Computes the log part of the digamma function.
         *
         * <p>This is not a true digamma function.
         *
         * @param x Argument.
         * @return log part of digamma(x)
         */
        static double log(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // large negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? x : Double.NaN;
            }

            double digamma = 0;
            if (x < 0) {
                // Use reflection formula to fall back into positive values.
                digamma = -Math.PI / Math.tan(Math.PI * x);
                x = 1 - x;
            }

            if (x <= S_LIMIT) {
                // Use method 5 from Bernardo AS103, accurate to O(x).
                return digamma - GAMMA - 1 / x;
            }

            // ----
            // From this point the function only compute the log term
            // ----

            // Avoid the while loop
            if (x < C_LIMIT) {
                x += (int) (1 + C_LIMIT - x);
            }

            digamma += Math.log(x);

            return digamma;
        }
    }

    /**
     * Copy of Commons Numbers Trigamma from version 1.3.
     * This includes fixes applied before the 1.4 release to the validated domain (NUMBERS-219);
     * and to use the reflection formula (NUMBERS-209).
     */
    private static final class Trigamma13 {
        /** C limit. */
        private static final double C_LIMIT = 49;

        /** S limit. */
        private static final double S_LIMIT = 1e-5;
        /** Fraction. */
        private static final double F_1_6 = 1d / 6;
        /** Fraction. */
        private static final double F_1_30 = 1d / 30;
        /** Fraction. */
        private static final double F_1_42 = 1d / 42;

        /** Private constructor. */
        private Trigamma13() {
            // intentionally empty.
        }

        /**
         * Computes the trigamma function.
         *
         * <p>Special cases:
         * <ul>
         * <li>If the argument is infinity, then the result is 0.</li>
         * <li>If the argument is negative infinity, then the result is nan.</li>
         * <li>If the argument is nan, then the result is nan.</li>
         * </ul>
         *
         * @param x Argument.
         * @return trigamma(x) to within {@code 1e-8} relative or absolute error whichever is larger.
         */
        static double value(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // large negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? 0 : Double.NaN;
            }

            if (x <= 0) {
                // Use the reflection formula:
                // trigamma(x) + trigamma(1 - x) = pi^2 / sin^2(pi * x)
                // to fall back into positive values. Without this the increment
                // loop below does not terminate for large negative x: when
                // ulp(x) >= 2 the update x += 1 leaves x unchanged (or advances
                // a single step and then sticks) and x < C_LIMIT holds forever.

                // negative integers are poles
                if (Math.rint(x) == x) {
                    return Double.POSITIVE_INFINITY;
                }

                final double s = Math.sin(Math.PI * x);
                return Math.PI * Math.PI / (s * s) - value(1 - x);
            }

            if (x <= S_LIMIT) {
                return 1 / (x * x);
            }

            double trigamma = 0;
            while (x < C_LIMIT) {
                trigamma += 1 / (x * x);
                x += 1;
            }

            final double inv = 1 / (x * x);
            //  1     1       1       1        1
            //  - + ----- + ----- - ------ + ------
            //  x   2 x^2   6 x^3   30 x^5   42 x^7
            trigamma += 1 / x + 0.5 * inv + inv / x * (F_1_6 - inv * (F_1_30 + F_1_42 * inv));

            return trigamma;
        }
    }

    /**
     * Copy implementation of the Boost digamma function.
     * This is here to test variations of the method to shift x into [1, 2].
     */
    private static final class BoostDigamma {
        /** 1/2. */
        private static final double HALF = 0.5;
        /** Limit for asymptotic expansion. */
        private static final double DIGAMMA_LARGE_LIMIT = 10;
        /** Constant for the digamma rational function. */
        private static final double Y = 0.99558162689208984;
        /** Constant for the digamma rational function. Value = 1.4616321446374059.
         * This is close to the root of the digamma function. 1.461632144968362341262660 (25 digits).
         * Note: ROOT1 + ROOT2 + ROOT3 ~ root of digamma in extended precision:
         * 1.46163214496836234126265954232572132846820517356040000000000 (added as BigDecimal). */
        private static final double ROOT1 = 1569415565.0 / 1073741824;
        /** Constant for the digamma rational function. Value = 3.309564688275257e-10. */
        private static final double ROOT2 = (381566830.0 / 1073741824) / 1073741824;
        /** Constant for the digamma rational function. */
        private static final double ROOT3 = 0.9016312093258695918615325266959189453125e-19;

        /** No instances. */
        private BoostDigamma() {}

        /**
         * Returns the digamma function.
         *
         * @param x Argument.
         * @return digamma(x)
         */
        static double digamma(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? x : Double.NaN;
            }
            return digammaImp(x);
        }

        /**
         * Implementation for the digamma function.
         * Adapted from {@code boost::math::detail::digamma_imp}.
         *
         * @param x Argument.
         * @return digamma(x)
         */
        private static double digammaImp(double x) {
            double result = 0;
            //
            // Check for negative arguments and use reflection:
            //
            if (x <= -1) {
                // Reflect:
                x = 1 - x;
                // Argument reduction for tan:
                double remainder = x - Math.floor(x);
                // Shift to negative if > 0.5:
                if (remainder > HALF) {
                    remainder -= 1;
                } else if (remainder == 0) {
                    // Digamma is undefined at negative integers.
                    // As x -> -n: digamma(x) -> +infty from one side and -infty from the other.
                    return Double.NaN;
                }
                result = Math.PI / Math.tan(Math.PI * remainder);
            } else if (x == 0) {
                // Undefined at the pole:
                // As -x -> 0: digamma(x) -> +infty
                // As x -> 0: digamma(x) -> -infty
                return Double.NaN;
            }
            //
            // If we're above the lower-limit for the
            // asymptotic expansion then use it:
            //
            if (x >= DIGAMMA_LARGE_LIMIT) {
                result += digammaImpLarge(x);
            } else {
                // ----
                // From this point the function only computes the shift but
                // does not update the sum using 1/x.
                // ----

                //
                // If x > 2 reduce to the interval [1,2]:
                //
                while (x > 2) {
                    x -= 1;
                    // result += 1 / x;
                }
                //
                // If x < 1 use recurrence to shift to > 1:
                //
                while (x < 1) {
                    // result -= 1 / x;
                    x += 1;
                }
                result += digammaImp1to2(x);
            }
            return result;
        }

        /**
         * Implementation for the digamma function using the asymptotic expansion.
         * Adapted from {@code boost::math::detail::digamma::digamma_imp_large} for
         * 53-bit precision.
         *
         * @param x Argument {@code >= 10}
         * @return digamma(x)
         */
        private static double digammaImpLarge(double x) {
            x -= 1;
            double result = Math.log(x);
            // Change from boost is to subtract the smaller terms first:
            // original: result += 1 / (2 * x)
            // boost::tools::evaluate_polynomial
            final double z = 1 / (x * x);
            double p;
            p =   -0.4432598039215686;
            p =   0.08333333333333333 + p * z;
            p = -0.021092796092796094 + p * z;
            p =  0.007575757575757576 + p * z;
            p = -0.004166666666666667 + p * z;
            p =  0.003968253968253968 + p * z;
            p = -0.008333333333333333 + p * z;
            p =   0.08333333333333333 + p * z;
            // Combine smaller terms first, then add log(x - 1)
            // original: result -= z * p
            result += 1 / (2 * x) - z * p;
            return result;
        }

        /**
         * Implementation for the digamma function with {@code 1 <= x <= 2}.
         * Adapted from {@code boost::math::detail::digamma_imp_1_2} for
         * 53-bit precision.
         *
         * @param x Argument in [1, 2]
         * @return digamma(x)
         */
        private static double digammaImp1to2(double x) {
            // Root 1 - 3 are successively smaller. Subtract sequentially.
            double g = x - ROOT1;
            g -= ROOT2;
            g -= ROOT3;
            final double xm1 = x - 1;
            // boost::tools::evaluate_polynomial P and Q
            double p;
            p = -0.0020713321167745952;
            p =  -0.045251321448739056 + p * xm1;
            p =   -0.28919126444774784 + p * xm1;
            p =   -0.65031853770896507 + p * xm1;
            p =   -0.32555031186804491 + p * xm1;
            p =    0.25479851061131551 + p * xm1;
            double q;
            q = -0.55789841321675513e-6;
            q =   0.0021284987017821144 + q * xm1;
            q =    0.054151797245674225 + q * xm1;
            q =     0.43593529692665969 + q * xm1;
            q =      1.4606242909763515 + q * xm1;
            q =      2.0767117023730469 + q * xm1;
            q =                     1.0 + q * xm1;
            final double r = p / q;
            return g * Y + g * r;
        }

        /**
         * Computes the log part of the digamma function.
         *
         * <p>This is not a true digamma function.
         *
         * @param x Argument.
         * @return log part of digamma(x)
         */
        private static double log(double x) {
            if (!Double.isFinite(x)) {
                // Function oscillates between infinite poles at integer values for
                // negative x; the -infinity limit is undefined
                return x == Double.POSITIVE_INFINITY ? x : Double.NaN;
            }
            double result = 0;
            //
            // Check for negative arguments and use reflection:
            //
            if (x <= -1) {
                // Reflect:
                x = 1 - x;
                // Argument reduction for tan:
                double remainder = x - Math.floor(x);
                // Shift to negative if > 0.5:
                if (remainder > 0.5) {
                    remainder -= 1;
                } else if (remainder == 0) {
                    // Digamma is undefined at negative integers.
                    // As x -> -n: digamma(x) -> +infty from one side and -infty from the other.
                    return Double.NaN;
                }
                result = Math.PI / Math.tan(Math.PI * remainder);
            } else if (x == 0) {
                // Undefined at the pole:
                // As -x -> 0: digamma(x) -> +infty
                // As x -> 0: digamma(x) -> -infty
                return Double.NaN;
            }
            // ----
            // From this point the function only computes the log term.
            // If the log term is not required then no shift is done.
            // This can be used to assess the cost of the shift and evaluation when
            // comparing to the full function, and the shift-only
            // copy implementation that does not update the sum using 1/x terms
            // but still evaluates x in [1, 2].
            // ----
            if (x >= DIGAMMA_LARGE_LIMIT) {
                result += Math.log(x - 1);
            }
            // else shift to x in [1, 2] and compute
            return result;
        }
    }

    /**
     * Apply the function to all the numbers.
     *
     * @param numbers Numbers.
     * @param fun Function.
     * @param bh Data sink.
     */
    private static void apply(double[] numbers, DoubleUnaryOperator fun, Blackhole bh) {
        for (int i = 0; i < numbers.length; i++) {
            bh.consume(fun.applyAsDouble(numbers[i]));
        }
    }

    // Benchmark methods.
    // Benchmarks use function references to perform different operations on the numbers.

    /**
     * Benchmark the function.
     *
     * @param data Test data.
     * @param method Operator method.
     * @param bh Data sink.
     */
    @Benchmark
    public void operator(NumberData data, MethodData method, Blackhole bh) {
        apply(data.getNumbers(), method.getOperator(), bh);
    }
}
