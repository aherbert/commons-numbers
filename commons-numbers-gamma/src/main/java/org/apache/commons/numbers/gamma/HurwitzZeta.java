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
package org.apache.commons.numbers.gamma;

import java.util.Arrays;
import org.apache.commons.numbers.core.DD;
import org.apache.commons.numbers.core.DDMath;

/**
 * <a href="https://en.wikipedia.org/wiki/Hurwitz_zeta_function">
 * Hurwitz zeta</a> function.
 *
 * <p>\[ \zeta(s, a) = \sum_{k=0}^\infty \frac{1}{(k+a)^s} \]
 *
 * <p>The function is formally defined for complex variable \( s \) with \( \mathrm{Re}(s) \gt 1 \)
 * and real \( a \ne 0, -1, -2, \ldots, -n \). This series is absolutely convergent for the given
 * values of \( s \) and \( a \). Note the special case \( zeta(s, 1) \) is the
 * {@link RiemannZeta Riemann zeta} function. This implementation uses real-valued \( s \gt 1 \).
 *
 * <p>The implementation is performed by spitting the integral into two parts and using
 * the Euler-Maclaurin formula to approximate the second integral \( \mathrm{I} + \mathrm{T} + \mathrm{R} \)
 * with a continuous integral \( \mathrm{I} \), a tail \( \mathrm{T} \), and a residual error term
 * \( \mathrm{R} \) (not computed).
 *
 * <p>\[ \begin{aligned}
 * \zeta(s, a) &amp;= \sum_{k=0}^{N-1} \frac{1}{(k+a)^s} + \sum_{k=N}^\infty \frac{1}{(k+a)^s} \\
 *            &amp;= \mathrm{S} + \left[ \mathrm{I} + \mathrm{T} + \mathrm{R} \right] \\
 * \mathrm{I} &amp;= \int_N^\infty \frac{1}{(a+t)^s} dt = \frac{(a+N)^{1-s}}{s-1} \\
 * \mathrm{T} &amp;= \frac{1}{(a+N)^s} \left( \frac{1}{2} + \sum_{k=1}^M \frac{B_{2k}}{2k!} \frac{(s)_{2k-1}}{(a+N)^{2k-1}} \right) \\
 * \mathrm{R} &amp;= - \int_N^\infty \frac{\tilde{B}_{2M}(t)}{2M!} \frac{(s)_{2M}}{(a+t)^{s+2M}} dt \end{aligned} \]
 *
 * <p>\( B_{2k} \) is a Bernoulli number; \( \tilde{B}_{2M}(t) \) is a generalized Bernoulli number;
 * and \( (s)_n \) is the rising factorial Pochammer function:
 *
 * <p>\[ (s)_n = \prod_{i=0}^{n-1} (s+i) \]
 *
 * <p>These formulas for the real-valued \( s \) are provided in Johansson (2015) as
 * equations 5-9. The implementation omits the residual term \( \mathrm{R} \).
 *
 * <p>The integrals are well defined when \( a + N \gt 0 \). Negative \( a \) may require the
 * sum \( S \) of a large number of terms \( N \). The implementation uses a difference
 * of zeta evaluations to compute the sum \( S \) over
 * \( a, a+1, a+2, \cdots, a - \lceil a \rceil  \) to avoid long runtime times.
 *
 * <p>Negative \( a \) with an odd integer power \( s \) requires summation of negative and positive
 * terms and cancellation reduces accuracy as \( a - \lceil a \rceil \to -\frac{1}{2} \).
 * Total cancellation is handled using the identity:
 *
 * <p>\[ \zeta(s, -\frac{1}{2} - n) = \zeta(s, n + \frac{3}{2}) \]
 *
 * <p>for \( s = 3, 5, 7, \ldots \) and \( n = 0, -1, -2, \ldots \)
 *
 * <p>References
 * <ol>
 * <li>Johansson (2015)
 * Rigorous high-precision computation of the Hurwitz zeta function and its derivatives
 * <a href="https://link.springer.com/article/10.1007/s11075-014-9893-1">Numerical Algorithms (69) 253–270</a></li>
 * <li><a href="https://en.wikipedia.org/wiki/Hurwitz_zeta_function">Hurwitz zeta function (Wikipedia)</a></li>
 * <li><a href="https://en.wikipedia.org/wiki/Euler%E2%80%93Maclaurin_formula">Euler–Maclaurin formula (Wikipedia)</a></li>
 * <li><a href="https://en.wikipedia.org/wiki/Bernoulli_number">Bernoulli number (Wikipedia)</a></li>
 * <li><a href="https://en.wikipedia.org/wiki/Falling_and_rising_factorials">Rising and falling factorials (Wikipedia)</a></li>
 * </ol>
 *
 * @see RiemannZeta
 * @since 1.4
 */
public final class HurwitzZeta {
    /** Number of terms of the series summation S. */
    private static final int N = 9;
    /** Convergence epsilon for the sum of the tail function. This prevents summation
     * of terms that do not affect the final result. Equal to 2^-53. */
    private static final double EPS = 0x1.0p-53;
    /** Mutiplicative inversion of the convergence epsilon. Equal to 2^53. */
    private static final double INV_EPS = 0x1.0p53;
    /** Asymptotic threshold for large {@code a}. Used when {@code a + N} is not accurate.
     * Equal to 2^53 - N. */
    private static final double LARGE_A = (1L << 53) - N;
    /** Threshold used when {@code a < 1} where {@code a^-s} == zeta(s, a) with the smallest supported
     * {@code s = 1 + 2^-52}. Set using zeta(s, 1) * 2^54. No possible evaluation of
     * all remaining terms for any s could be added to {@code a^-s} when {@code a < 1}.
     *
     * <p>Approximately equal to 2^106 thus remaining terms (all less than 1) cannot be added in
     * double-double precision. Allows usage in the negative zeta implementation on the initial
     * terms. */
    private static final double LARGE_T0 = 8.11296384146066920939820185959e+31;
    /** Asymptotic threshold for large {@code s} when {@code a < 0}.
     * Used when any possible {@code [0.5+/-ulp(a)]^-s - [0.5-/+ulp(a)]^-s} will overflow.
     * Using extended precision this is finite: {@code (0.5-2^-54)^-1065 - (0.5+2^-54)^-1065}.
     * Note that {@code 0.5^-1024} overflows and at least one term uses {@code x < 0.5}.
     * So 0.5 values separated by more than +/- 2^-54 will also overflow. */
    private static final int LARGE_S = 1067;
    /** Maximum {@code s} value where {@code zeta(s, 1.5) > 0}. Used in the reflection
     * identity when {@code s} is odd and {@code a} is negative and half-integer. */
    private static final int MAX_REFLECTION_S = 1837;
    /** 1/2. */
    private static final double HALF = 0.5;
    /** Maximum difference in exponents to support double-double addition. */
    private static final int MAX_EXP_DELTA = 106;
    /** Downscale to use on the initial rising factorial tail term (always above 1).
     * Equal to Double.MIN_NORMAL: 2^-1022. */
    private static final double DOWNSCALE = 0x1.0p-1022;
    /** Upscale to use on the initial probability tail term (always below 1).
     * Multiplicative inverse of Double.MIN_NORMAL: 2^1022. */
    private static final double UPSCALE = 0x1.0p1022;
    /** Threshold for z in the Taylor series approximation of (1 + t)^z.
     * Set where |z * t| < 2^-26  =>  z < 2^-26 / 2^-52 = 2^26. */
    private static final double Z_LIMIT = 0x1p26;
    /** Minimum normalised exponent for a double.
     * Note 2^-1075 is the minimum non-zero double.
     * Equivalent to Double.MIN_EXPONENT - 53. */
    private static final int MIN_DOUBLE_EXPONENT = -1075;

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}
     * in double-double (DD) precision.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * Provides M=51 terms. The rising factorial overflows after k=51 for s=1065.
     * The test suite uses max 9 before convergence on double data and
     * 20 on double-double data.
     */
    private static final DD[] FDD = {
        DD.ofSum(0.08333333333333333, 4.625929269271485E-18), // (1 / 6) / 2!
        DD.ofSum(-0.001388888888888889, 5.300543954373577E-20), // (-1 / 30) / 4!
        DD.ofSum(3.306878306878307E-5, -2.2300719288557665E-21), // (1 / 42) / 6!
        DD.ofSum(-8.267195767195768E-7, 3.457597454003665E-23), // (-1 / 30) / 8!
        DD.ofSum(2.08767569878681E-8, -1.2073450591132599E-24), // (5 / 66) / 10!
        DD.ofSum(-5.284190138687493E-10, 3.517096671929869E-27), // (-691 / 2730) / 12!
        DD.ofSum(1.3382536530684679E-11, -2.828354019907999E-29), // (7 / 6) / 14!
        DD.ofSum(-3.3896802963225827E-13, -1.4986928409964295E-29), // (-3617 / 510) / 16!
        DD.ofSum(8.586062056277845E-15, -6.05252374381974E-31), // (43867 / 798) / 18!
        DD.ofSum(-2.174868698558062E-16, 4.961617782549996E-33), // (-174611 / 330) / 20!
        DD.ofSum(5.5090028283602295E-18, -1.49827152194499E-35), // (854513 / 138) / 22!
        DD.ofSum(-1.3954464685812522E-19, -1.0350590497256251E-35), // (-236364091 / 2730) / 24!
        DD.ofSum(3.534707039629467E-21, 1.894231142684204E-37), // (8553103 / 6) / 26!
        DD.ofSum(-8.953517427037546E-23, -5.728752743153026E-39), // (-23749461029 / 870) / 28!
        DD.ofSum(2.267952452337683E-24, 1.3043458462619563E-40), // (8615841276005 / 14322) / 30!
        DD.ofSum(-5.744790668872202E-26, 1.663242973708004E-43), // (-7709321041217 / 510) / 32!
        DD.ofSum(1.455172475614865E-27, -5.613265715443096E-44), // (2577687858367 / 6) / 34!
        DD.ofSum(-3.6859949406653103E-29, 1.0778256413554197E-45), // (-26315271553053477373 / 1919190) / 36!
        DD.ofSum(9.336734257095045E-31, -3.9347970210731877E-47), // (2929993913841559 / 6) / 38!
        DD.ofSum(-2.36502241570063E-32, 2.0347170931532494E-49), // (-261082718496449122051 / 13530) / 40!
        DD.ofSum(5.990671762482134E-34, 1.6265467158179092E-50), // (1520097643918070802691 / 1806) / 42!
        DD.ofSum(-1.5174548844682903E-35, 5.493014407946745E-52), // (-27833269579301024235023 / 690) / 44!
        DD.ofSum(3.843758125454189E-37, -3.685053096067968E-53), // (596451111593912163277961 / 282) / 46!
        DD.ofSum(-9.736353072646691E-39, 2.258059165188444E-55), // (-5609403368997817686249127547 / 46410) / 48!
        DD.ofSum(2.466247044200681E-40, -1.505641802268162E-56), // (495057205241079648212477525 / 66) / 50!
        DD.ofSum(-6.247076741820743E-42, -2.7106815859687654E-58), // (-801165718135489957347924991853 / 1590) / 52!
        DD.ofSum(1.5824030244644914E-43, 2.545428531496969E-60), // (29149963634884862421418123812691 / 798) / 54!
        DD.ofSum(-4.008273685948936E-45, -2.2124211668946826E-61), // (-2479392929313226753685415739663229 / 870) / 56!
        DD.ofSum(1.0153075855569557E-46, -9.404269751258486E-63), // (84483613348880041862046775994036021 / 354) / 58!
        DD.ofSum(-2.5718041582418717E-48, -6.537655454012542E-65), // (-1215233140483755572040304994079820246041491 / 56786730) / 60!
        DD.ofSum(6.514456035233815E-50, -2.763626172529861E-66), // (12300585434086858541953039857403386151 / 6) / 62!
        DD.ofSum(-1.6501309906896525E-51, 3.1794529475063687E-68), // (-106783830147866529886385444979142647942017 / 510) / 64!
        DD.ofSum(4.179830628539476E-53, 2.617556823159939E-69), // (1472600022126335654051619428551932342241899101 / 64722) / 66!
        DD.ofSum(-1.058763466770291E-54, 6.6915528436035195E-71), // (-78773130858718728141909149208474606244347001 / 30) / 68!
        DD.ofSum(2.6818791912607708E-56, -8.70695425146146E-73), // (1505381347333367003803076567377857208511438160235 / 4686) / 70!
        DD.ofSum(-6.793279351107421E-58, 2.795667911354165E-74), // (-5827954961669944110438277244641067365282488301844260429 / 140100870) / 72!
        DD.ofSum(1.7207577616681404E-59, 4.65433497191727E-76), // (34152417289221168014330073731472635186688307783087 / 6) / 74!
        DD.ofSum(-4.358730329348894E-61, 2.8840522874209336E-77), // (-24655088825935372707687196040585199904365267828865801 / 30) / 76!
        DD.ofSum(1.1040792903684666E-62, 6.624841731022409E-79), // (414846365575400828295179035549542073492199375372400483487 / 3318) / 78!
        DD.ofSum(-2.7966655133781345E-64, 2.628041826403209E-81), // (-4603784299479457646935574969019046849794257872751288919656867 / 230010) / 80!
        DD.ofSum(7.084036501679471E-66, -5.026235239023924E-82), // (1677014149185145836823154509786269900207736027570253414881613 / 498) / 82!
        DD.ofSum(-1.794407408289224E-67, 1.5372719769275798E-84), // (-2024576195935290360231131160111731009989917391198090877281083932477 / 3404310) / 84!
        DD.ofSum(4.545287063611096E-69, 9.87696151726261E-87), // (660714619417678653573847847426261496277830686653388931761996983 / 6) / 86!
        DD.ofSum(-1.1513346631982051E-70, -7.192856523313341E-87), // (-1311426488674017507995511424019311843345750275572028644296919890574047 / 61410) / 88!
        DD.ofSum(2.9163647710923614E-72, -3.8911087510195904E-89), // (1179057279021082799884123351249215083775254949669647116231545215727922535 / 272118) / 90!
        DD.ofSum(-7.387238263497337E-74, -6.923136687699924E-90), // (-1295585948207537527989427828538576749659341483719435143023316326829946247 / 1410) / 92!
        DD.ofSum(1.8712093117637953E-75, 1.5886680102062367E-92), // (1220813806579744469607301679413201203958508415202696621436215105284649447 / 6) / 94!
        DD.ofSum(-4.739828557761799E-77, -9.517121002177184E-94), // (-211600449597266513097597728109824233673043954389060234150638733420050668349987259 / 4501770) / 96!
        DD.ofSum(1.2006125993354507E-78, -1.109850335891779E-95), // (67908260672905495624051117546403605607342195728504487509073961249992947058239 / 6) / 98!
        DD.ofSum(-3.0411872415142924E-80, 5.125117133572647E-97), // (-94598037819122125295227433069493721872702841533066936133385696204311395415197247711 / 33330) / 100!
        DD.ofSum(7.703417274705106E-82, 3.948211996024456E-99), // (3204019410860907078243020782116241775491817197152717450679002501086861530836678158791 / 4326) / 102!
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * This is the high part of {@link #FDD}.
     */
    private static final double[] F;

    static {
        F = Arrays.stream(FDD).mapToDouble(DD::hi).toArray();
    }

    /**
     * Compute the power function {@code (x+n)^z} in double precision.
     */
    @FunctionalInterface
    interface PowOperator {
        /**
         * Applies this operator to the given operands.
         *
         * @param x Argument x
         * @param n Argument n
         * @param z Exponent z
         * @return the operator result
         */
        double apply(double x, int n, double z);
    }

    /**
     * Compute the power function {@code x^n} in double-double precision.
     * Results are returned as scaled values to prevent under and overflow.
     */
    @FunctionalInterface
    interface ScaledPowOperator {
        /**
         * Compute number {@code x} raised to the power {@code n}.
         *
         * <p>The value is returned as fractional {@code f} and integral
         * {@code 2^exp} components.
         * <pre>
         * (x+xx)^n = (f+ff) * 2^exp
         * </pre>
         *
         * @param x Number.
         * @param n Power.
         * @param exp Result power of two scale factor (integral exponent).
         * @return Fraction part.
         * @see DD#frexp(int[])
         */
        DD apply(DD x, int n, long[] exp);
    }

    /** Context for the double-double zeta implementation. */
    private static class Context {
        /**
         * Context for double precision.
         * <ul>
         * <li>N = 9
         * <li>53-bits to stop the tail series
         * <li>Regular {@link DD} power function for the most important terms
         * </ul>
         */
        static final Context DOUBLE = new Context(9, 53, 0);
        /**
         * Context for double-double precision.
         * <ul>
         * <li>N = 15
         * <li>106-bits to stop the tail series
         * <li>Accurate {@link DDMath} power function for the most important terms
         * </ul>
         * <p>Note that accuracy is limited to approximately 105-bits.
         */
        static final Context DOUBLE_DOUBLE = new Context(15, 106, 1);

        /** N. */
        private final int n;
        /** Precision in bits for convergence of the tail series. */
        private final int precision;
        /** Option for DD power function. */
        private final int powOption;

        /**
         * Create an instance.
         *
         * @param n Number of terms to sum
         * @param precision Convergence precision in bits for the tail series T
         * @param powOption Option for power function
         */
        Context(int n, int precision, int powOption) {
            this.n = n;
            this.precision = precision;
            this.powOption = powOption;
        }

        /**
         * Gets N.
         * @return n
         */
        int getN() {
            return n;
        }

        /**
         * Gets the precision in bits for the tail series T.
         * @return the precision
         */
        int getPrecision() {
            return precision;
        }

        /**
         * Get the function to compute {@code x^y} terms not used in the series sum.
         * May return an extended precision power function.
         * @return function
         */
        ScaledPowOperator getPow() {
            return (powOption & 1) == 0 ?
                DD::pow :
                DDMath::pow;
        }
    }

    /** No instances. */
    private HurwitzZeta() {}

    /**
     * Computes the value of \( \zeta(s, a) \) over the domain \( s \gt 1 \) and
     * \( a \ne 0, -1, -2, \ldots, -n \).
     *
     * <p>Negative \( a \) is supported when \( s \) is an integer due to the support
     * for powers of negative real numbers.
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
     * <p>Note: This function is evaluated for \( s \gt 1 \) where the series is convergent.
     * The special case of the Riemann zeta function with \( a = 1 \) may be evaluated for
     * \( s \lt 1 \) using {@link RiemannZeta#value(double)}.
     *
     * @param s Argument.
     * @param a Argument.
     * @see Math#pow(double, double)
     * @return \( \zeta(s, a) \)
     */
    public static double value(double s, double a) {
        if (Double.isNaN(s) || Double.isNaN(a) || s < 1 || a == Double.NEGATIVE_INFINITY) {
            return Double.NaN;
        }
        // s > 1
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
            return zetaNegativeImp(s, a, ca);
        }
        // Use the faster Riemann zeta function if applicable (accuracy is similar).
        // Creates a consistent output between zeta(s, 1) and zeta(s).
        // Done after domain validation, e.g.
        // this function will return NaN for s < 1 even when a==1.
        if (a == 1) {
            return RiemannZeta.value(s);
        }
        return zetaImp(s, a);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)} for positive {@code a}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @return zeta(s, a)
     */
    // package-private for testing using a == 1 (Riemann zeta function)
    static double zetaImp(double s, double a) {
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

        // Check the extended precision power can be used.
        // If s is large or there is no round-off from (a+n) use standard precision.
        final PowOperator pow = s >= Z_LIMIT || twoSumLow(a, N, a + N) == 0 ?
            HurwitzZeta::stdPowxpn :
            HurwitzZeta::extPowxpn;

        double p = pow.apply(a, N, -s);

        // Initialise sum with the first tail term
        DD sum = DD.of(0.5 * p);
        // S : k in [0, n-1]
        for (int k = N - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum = sum.add(pow.apply(a, k, -s));
        }

        // I
        final double ti = pow.apply(a, N, 1 - s) / (s - 1);

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (t0 > ti) {
            sum = sum.add(ti).add(t0);
        } else {
            sum = sum.add(t0).add(ti);
        }

        // T
        // Rising factorial term : (s)_{2k-1}
        // Note: When s is large the loop exits before the rising factorial overflows.
        // (a+n) >= 9 : 9^-340 = 0 : max (2k-1+s) = 339
        // (339)_{2k-1}; k=50 = Pochammer(339, 99) = 1.5e256
        // The rising factorial will not overflow for k <= 50 before (a+n)^-(2k-1+s) is
        // zero. This is within the length of table F.
        //
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        //
        // If a or s are large then (a+n)^-(2k-1+s) will underflow.
        // We can downscale the rising factorial (always above 1)
        // and upscale p (always below 1) by the same amount.
        // The loop will still exit immediately if p = (a+n)^-s = 0
        // (zeta evaluation is limited by Math.pow).
        // Scaling could dynamically use the exponent of f and p. Here we use
        // fixed scaling which can be pre-applied to a possible sub-normal
        // power term before division by (a+n). This allows a few more terms
        // of the tail to be computed when the zeta result is sub-normal.
        double f = s * DOWNSCALE;
        // 2k - 1
        double k2 = 1;

        // Use to divide by (a+n)^2
        final double apn = pow.apply(a, N, -2);

        // Sub-normal p reduces accuracy:
        // Test if (a+n)^-(1+s) >= 2^-1022
        if (-(1 + s) * Math.getExponent(a + N) >= Double.MIN_EXPONENT) {
            p = UPSCALE * pow.apply(a, N, -s - 1);
        } else {
            // If (a+n)^-s is non-zero this allows division using normal numbers
            p = UPSCALE * p / (a + N);
        }

        double stop = sum.hi() * EPS;
        double rescale = 1;
        if (stop < Double.MIN_NORMAL) {
            // Compute a scaled tail sum so convergence is on a normal number.
            // Upscaling sum cannot overflow due to overflow check on t0.
            stop = sum.hi() * INV_EPS;
            // Either scale p down or scale f up
            f *= INV_EPS * INV_EPS;
            rescale = EPS * EPS;
        }

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Note: an extended precision sum here has no effect on the final result.
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
        return sum.add(tsum * rescale).hi();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)} for negative {@code a}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be negative.
     *
     * @param s Argument {@code s > 1} and integer
     * @param a Argument {@code a < 0}
     * @param ca Ceil(a)
     * @return zeta(s, a)
     */
    private static double zetaNegativeImp(double s, double a, double ca) {
        // a < 0 (non-integer) and s is a positive integer.
        // If s is odd then the negative series sum will be negative
        // and the addition of zeta(s, x > 0) has cancellation.
        // This is largest when a is close to half-integer.

        // Check case of total cancellation.
        final boolean odd = SpecialMath.isOdd(s);
        final double x = a - ca;

        // Exponent for fractional representation
        final int[] exp = {0};

        // Intentional float comparison
        if (odd && x == -HALF) {
            if (s <= MAX_REFLECTION_S) {
                // Use extended precision but evaluated with precision for a double result
                final DD r = zetaImp((int) s, DD.ONE.subtract(a), null, exp, Context.DOUBLE);
                return Math.scalb(r.hi(), exp[0]);
            }
            // All possible zeta(s, 1-a >= 1.5) underflow
            return 0.0;
        }

        // Handle cancellation as x -> 0.5 using double-double precision.
        // Eventually [0.5 +/- ulp(a)]^-s - [0.5 -/+ ulp(a)]^-s will overflow.
        if (s >= LARGE_S) {
            // Use the dominant term closest to zero
            return odd && x > -HALF ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        }

        // Integer s power function supported by double-double (DD)
        final int is = (int) s;

        // Compute the two terms either side of zero:
        // -1 < x < 0 < x + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        final DD xn = DD.of(x);
        final DD xp = DD.ONE.add(x);
        final long[] expn = {0};
        final long[] expp = {0};
        DD pn = DDMath.pow(xn, -is, expn);
        DD pp = DDMath.pow(xp, -is, expp);
        final long delta = expp[0] - expn[0];
        if (Math.abs(delta) > MAX_EXP_DELTA) {
            // Cannot add these terms given the largest is [0.5, 1.0) with ulp 2^-106.
            // All other individual terms are smaller so further DD computation is not possible.
            return delta > 0 ?
                Math.scalb(pp.hi(), intExponent(expp[0])) :
                Math.scalb(pn.hi(), intExponent(expn[0]));
        }

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 1.65; individual terms x^-s < 1.
        // Exit early if remaining terms cannot be added.
        // The result can be finite even if the rescaled terms are infinite. However
        // we do not support further computation from infinite terms so check
        // if anything can be added to these terms.
        if (Math.max(expn[0], expp[0]) > MAX_EXP_DELTA) {
            // Limit of double-double arithmetic
            // Add a scaled term avoiding overflow by re-scaling after the sum
            return Math.scalb(
                pn.add(pp.scalb((int) delta)).hi(),
                intExponent(expn[0])
            );
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) +/- [ zeta(s, -x) - zeta(s, 1 - a) ]
        // z +/- [ za - zb]

        // We have the two largest power terms for each side.
        // The remaining terms are increasingly smaller. Computing with
        // double-double (DD) precision for the zeta evaluations should handle all
        // but the most extreme cancellation close to the root.

        // Evaluate zeta with extra precision when cancellation occurs
        final Context c = odd ? Context.DOUBLE_DOUBLE : Context.DOUBLE;

        exp[0] = intExponent(expp[0]);
        DD z = zetaImp(is, xp, pp, exp, c);
        int ez = exp[0];

        DD za;
        DD zb;
        int eza;
        int ezb;
        // A single call to zeta uses many pow operations;
        // use a direct sum when zeta will use more.
        // Prepare the exponent for the term we have.
        exp[0] = intExponent(expn[0]);
        if (x - a > 2 * (c.getN() + 2)) {
            // Note: The difference (za - zb) incurs cancellation.
            // This should not be an issue as x in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 1)  = 1.6449
            // zeta(2, 2)  = 0.6449
            // zeta(2, 30) = 0.0339
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            za = zetaImp(is, xn.negate(), pn.abs(), exp, c);
            eza = exp[0];
            zb = zetaImp(is, DD.ONE.subtract(a), null, exp, c);
            ezb = exp[0];
            // Both terms are positive. Correct the sign for final addition.
            if (odd) {
                za = za.negate();
            } else {
                zb = zb.negate();
            }
        } else {
            // Sum terms in ascending order of magnitude
            // Using a double to track the iterations is fine as (a+n) is exact until > x.
            za = pn;
            eza = exp[0];
            zb = DD.ZERO;
            ezb = 0;
            final long[] e = {0};
            for (double aa = a; aa < x; aa += 1.0) {
                final DD r = DD.of(aa).pow(-is, e);
                zb = add(zb, ezb, r, intExponent(e[0]), exp);
                ezb = exp[0];
            }
            // Sort for the summation by magnitude
            if (ezb > eza) {
                final DD tmp = za;
                za = zb;
                zb = tmp;
                final int etmp = eza;
                eza = ezb;
                eza = etmp;
            }
        }

        // Sum in magnitude order
        // zb < za : check z against za
        if (ez > eza) {
            // za,zb < z
            za = add(za, eza, zb, ezb, exp);
            eza = exp[0];
        } else {
            // z,zb < za
            z = add(z, ez, zb, ezb, exp);
            ez = exp[0];
        }
        z = add(za, eza, z, ez, exp);
        return Math.scalb(z.hi(), exp[0]);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)} for positive {@code a}
     * using double-double (DD) arithmetic.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code 1 < s <= 1837}
     * @param a Argument {@code a > 0}
     * @param a0 {@code a^-s} (can be null)
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    private static DD zetaImp(int s, DD a, DD a0, int[] exp, Context c) {
        // Note:
        // The power function uses a long exponent to support the full range
        // of a double to an integer power: (2^1024)^(2^31) = 2^(2^10 * 2^31) = 2^(2^41)
        // Here we can clip the long to an integer. The final result is either
        // infinite or zero if the exponent is outside [-1075, 1024) so special
        // routines are used to clip the long exponents to integer and add them
        // without exceeding the range of an integer.
        // Only the first term is above 1; all other terms are < 1 and can be added
        // exactly even if sub-normal.

        // For the exponent of the power result
        final long[] e = {0};

        // Note:
        // Computing with the standard DD pow function has enough accuracy to not accumulate
        // error to the zeta result when a < 1 as the first term (provided) dominates the result.
        // When a >> 1 then the terms are all similar magnitude and we may benefit from
        // DDmath for the zeta result. However this result is then subtracted from a much larger
        // term on the negative a computation; or only used for double precision on the complete
        // cancellation path. In either case the full ~105 bits in the power result are not
        // required from the pow series.

        // Power function based on the context
        final ScaledPowOperator pow = c.getPow();

        // First term may be provided
        final DD t0;
        final int et0;
        if (a0 == null) {
            t0 = pow.apply(a, -s, e);
            et0 = intExponent(e[0]);
        } else {
            t0 = a0;
            et0 = exp[0];
        }

        // Asymptotic Behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0, the I term and the first term of T.
        if (Math.getExponent(a.hi()) >= c.getPrecision()) {
            // a^(1-s) / (s - 1) + a^-s / 2
            // Note: If a^-s is sub-normal then a^(1-s)=(a^-s * a) may not be
            // sub-normal when s is small as a is large.
            DD r = pow.apply(a, 1 - s, e);
            long er = e[0];
            if (s > 2) {
                r = r.divide(s - 1).frexp(exp);
                er += exp[0];
            }
            // Add 0.5 * t0 so adjust the exponent
            return add(r, intExponent(er), t0, addExponent(et0, -1), exp);
        }

        // Note:
        // This function is never called with 0 < a < 1 unless t0 is provided and finite.
        // No check for overflow/large t0 is required.
        // Underflow can occur for large s if a is large. The function is not called
        // when s is extremely large so this case is handled after the
        // major part of the sum has been computed.

        // Now any (a+n)^-s cannot overflow and the sum cannot overflow.
        // However the terms may be sub-normal and we can add them using a fractional
        // representation with a separate base 2 exponent.

        final int n = c.getN();
        DD apn = a.add(n);

        // Initialise sum with the first tail term: 0.5 * (a+n)^-s
        DD sum = apn.pow(-s, e);
        int es = addExponent(intExponent(e[0]), -1);
        // S : k in [0, n-1]
        for (int k = n - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            DD r = a.add(k).pow(-s, e);
            sum = add(sum, es, r, intExponent(e[0]), exp);
            es = exp[0];
        }

        // I : (a+p)^(1-s) / (s-1)
        DD ti = pow.apply(apn, 1 - s, e);
        long eti = e[0];
        if (s > 2) {
            ti = ti.divide(s - 1).frexp(exp);
            eti += exp[0];
        }

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (et0 > eti) {
            sum = add(sum, es, ti, intExponent(eti), exp);
            sum = add(sum, exp[0], t0, et0, exp);
        } else {
            sum = add(sum, es, t0, et0, exp);
            sum = add(sum, exp[0], ti, intExponent(eti), exp);
        }
        es = exp[0];

        // If the sum of terms is far below min double value then stop.
        // Check if the current sum cannot be added to
        // Double.MIN_VALUE in the working precision.
        if (es < MIN_DOUBLE_EXPONENT - c.getPrecision()) {
            exp[0] = 0;
            return DD.ZERO;
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        // Terms cannot overflow using the fractional representation. This also
        // allows the tail computation to continue further when an individual
        // part of each term computation would underflow.
        // Note:
        // This has a drawback when s > (a+n) and the pochammer function (s)_2k grows
        // faster than the term (a+n)^-(2k-1-s) reduces. For large s and/or a the
        // power term would normally underflow and the tail is ignored. Here the
        // convergence must be possible for all expected values of s and a.

        // Rising factorial term : (s)_{2k-1}
        DD f = DD.of(s).frexp(exp);
        long ef = exp[0];

        // Initialise (a+n)^-(2k-1+s) to (a+n)^-(1+s)
        DD p = pow.apply(apn, -1 - s, e);
        long ep = e[0];
        // Divide by (a+n)^2 using multiplication
        apn = pow.apply(apn, -2, e);
        long eapn = e[0];

        // Set the scale for the tail series

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Use the tail epsilon to obtain a minimum exponent for terms.
        final long stop = (long) es - c.getPrecision();

        // s+2k-1
        long s2k = s + 1L;
        // Initialise with the first term.
        // Avoids addition of fractional representation to zero (which has a zero exponent).
        DD tsum = f.multiply(p).multiply(FDD[0]).frexp(exp);
        long et = ef + ep + exp[0];
        int ets = intExponent(et);
        for (int i = 1; i < FDD.length; i++) {
            if (et < stop) {
                break;
            }
            // p = (a+n)^-(2k-1+s)
            p = p.multiply(apn).frexp(exp);
            ep += exp[0] + eapn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            f = f.multiply(s2k * (s2k + 1L)).frexp(exp);
            ef += exp[0];
            s2k += 2L;
            // FDD is finite. f*p in [0, 2).
            final DD t = f.multiply(p).multiply(FDD[i]).frexp(exp);
            et = ef + ep + exp[0];
            tsum = add(tsum, ets, t, intExponent(et), exp);
            ets = exp[0];
        }
        return add(sum, es, tsum, ets, exp);
    }

    /**
     * Helper function to add the fractional representation of two DD numbers.
     *
     * @param a Addend a.
     * @param ae Exponent of a.
     * @param b Addend b.
     * @param be Exponent of b.
     * @param e Exponent of the result.
     * @return the fraction part of the result
     */
    private static DD add(DD a, int ae, DD b, int be, int[] e) {
        // Add if the overlap is within 106-bits
        final long diff = (long) ae - be;
        if (Math.abs(diff) > MAX_EXP_DELTA) {
            // Special case where the exponent of zero is 0.
            // This may be larger than the exponent of the other term.
            if (ae > be && !a.isZero()) {
                e[0] = ae;
                return a;
            }
            e[0] = be;
            return b;
        }
        final DD r = b.add(a.scalb((int) diff)).frexp(e);
        e[0] = addExponent(e[0], be);
        return r;
    }

    /**
     * Helper function to add the two exponents and maintain the [min, max]
     * range of an integer.
     *
     * @param a Exponent.
     * @param b Exponent.
     * @return (a + b) within the range of an integer
     * @see Integer#MIN_VALUE
     * @see Integer#MAX_VALUE
     */
    private static int addExponent(int a, int b) {
        final int r = a + b;
        // Hackers Delight 2-12:
        // Overflow if both arguments have the opposite sign of the result
        if (((a ^ r) & (b ^ r)) < 0) {
            // Use the sign of the overflow to choose the result
            return r < 0 ? Integer.MAX_VALUE : Integer.MIN_VALUE;
        }
        return r;
    }

    /**
     * Clip the long exponent to an integer range.
     *
     * <p>Note that use of a direct cast to {@code int} is a 32-bit truncation. If all
     * lower bits are zero but the high 32-bits contain non-zero bits then a cast is
     * incorrect. This method ensures {@code long} exponents generated from large
     * power computations are clipped to the range {@code [-(2^31), 2^31)}. This is
     * sufficient to detect infinite or zero values, and correctly scale results
     * using {@link Math#scalb(double, int)} or {@link DD#scalb(int)}.
     *
     * @param e Exponent.
     * @return the exponent as the closest integer
     */
    private static int intExponent(long e) {
        // Test if 32-bit truncation is the same number
        final int r = (int) e;
        if (r != e) {
            return e < 0 ? Integer.MIN_VALUE : Integer.MAX_VALUE;
        }
        return r;
    }

    /**
     * Standard precision {@code (x+n)^z}.
     *
     * @param x Argument x
     * @param n Argument n
     * @param z Exponent z
     * @return the result
     */
    private static double stdPowxpn(double x, int n, double z) {
        return Math.pow(x + n, z);
    }

    /**
     * Extended precision power function {@code (x+n)^z}.
     *
     * <p>This method is used to compute the power function terms of
     * the infinite series {@code (a+k)^-s}. It is limited to the expected arguments
     * of the zeta function use case:
     *
     * <ul>
     * <li>{@code x > 0}
     * <li>{@code n >= 1}
     * <li>{@code x+n >= 1} and finite
     * <li>{@code -(2^26) < z < -1} and finite
     * </ul>
     *
     * <p>This method has no advantage over {@code Math.pow(x+n, z)} when there is no
     * round-off component to (x+n). In the event that {@code (x+n)} is exact this method
     * will return the same result as {@code Math.pow(x+n, z)}.
     *
     * <p>Calling with {@code n=0} will consume resources for no benefit. More generally
     * this method should not be called when {@code x + n} does not create a change of
     * exponent between {@code x} and {@code x + n}, and may offer no benefit even if the
     * exponent does change depending on the arguments. For example if {@code |x| < 2^52}
     * and is an integer then addition of any {@code int} value will be exact.
     *
     * @param x Argument x ({@code 0 < x < 2^53})
     * @param n Argument n ({@code > 0})
     * @param z Exponent z ({@code -(2^26) < z < -1})
     * @return the result
     */
    // package-private for testing
    static double extPowxpn(double x, int n, double z) {
        // (s+ss)^z = s^z * (1+ss/s)^z
        //
        // Taylor series:
        // (1 + t)^z = 1 + zt + z(z-1)t^2/2! + z(z-1)(z-2)t^3/3! + ...
        // for |t| < 1 and real z

        // Compute the exact double-double sum.
        // Avoid generation of a DD intermediate by directly using twoSum.
        double s = x + n;
        double ss = twoSumLow(x, n, s);
        double t;

        // Result cases assuming (x+n) >= 1:
        // s = 1 creates 1.0
        // s > 1 can underflow to 0
        // The correction must always be in the same direction, i.e. make the result smaller.
        // Change a negative round-off to positive by moving 1 ulp from s to ss.
        // This is an unpredictable 50:50 branch.
        if (ss < 0) {
            // x > 0, n >= 1 : if ss < 0 then s > 1
            // Do not use Math.ulp. Here we know s is finite > 1. Use raw bits.
            // These methods are intrinsic candidates.
            t = Double.longBitsToDouble(Double.doubleToRawLongBits(s) - 1L);
            // This will be exact as the double+int was 84 bits maximum
            ss += s - t;
            s = t;
        }

        final double r = Math.pow(s, z);

        // 0 < t < 2^-52 (or < 2^-53 if no bit was moved)
        t = ss / s;

        // Compute (1 + t)^z - 1 using a Taylor series of 2 terms with Horner's rule.
        // Equivalent to expm1(z * log1p(t)). Accurate to z*t < 2^-26.
        // Result is approximately [-1.5e-8, 0] correcting r towards zero.
        t = z * t * (1 + (z - 1) * t * 0.5);

        return r + r * t;
    }

    /**
     * Compute the round-off of the sum of two numbers {@code a} and {@code b} using
     * Knuth two-sum algorithm. The values are not required to be ordered by magnitude,
     * i.e. the result is commutative {@code s = a + b == b + a}.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return the sum round-off
     */
    private static double twoSumLow(double a, double b, double x) {
        // (x, xx) = a + b
        // bVirtual = x - a
        // aVirtual = x - bVirtual
        // bRoundoff = b - bVirtual
        // aRoundoff = a - aVirtual
        // xx = aRoundoff + bRoundoff
        final double bVirtual = x - a;
        return (a - (x - bVirtual)) + (b - bVirtual);
    }
}
