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
     * of terms that do not affect the final result. */
    private static final double EPS = 0x1.0p-53;
    /** Asymptotic threshold for large {@code a}. Used when {@code a + N} is not accurate. */
    private static final double LARGE_A = (1L << 53) - N;
    /** Threshold used when {@code a < 1} where {@code a^-s} == zeta(s, a) with the smallest supported
     * {@code s = 1 + 2^-52}. Set using zeta(s, 1) * 2^54. No possible evaluation of
     * all remaining terms for any s could be added to {@code a^-s} when {@code a < 1}.
     * Approximately equal to 2^106 thus remaining terms (all less than 1) cannot be added in
     * double-double precision. */
    private static final double LARGE_T0 = 8.11296384146066920939820185959e+31;
    /** 0.5. */
    private static final double HALF = 0.5;
    /** Maximum exponent above which 0.5^-s will be infinity. */
    private static final double MAX_S = 1024;
    /** double-double NaN. */
    private static final DD DD_NAN = DD.of(Double.NaN);
    /** Maximum number of terms in the negative series summation.
     * Equal to the number of Math.pow calls in two zeta function evaluations. */
    private static final int MAX_TERMS = 2 * (N + 2);

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
     * Compute the power function {@code (x+y)^z}.
     */
    @FunctionalInterface
    interface PowOperator {
        /**
         * Applies this operator to the given operands.
         *
         * @param x the first operand
         * @param y the second operand
         * @param z the third operand
         * @return the operator result
         */
        double apply(int x, double y, double z);
    }

    /** No instances. */
    private HurwitzZeta() {}

    /**
     * Extended precision {@code (x+y)^z}.
     *
     * <p>Warning: This does not check all pow edge cases and
     * assumes {@code (x+y)} is finite.
     *
     * @param x the first operand
     * @param y the second operand
     * @param z the third operand
     * @return the result
     */
    private static double extendedPowNp(int x, double y, double z) {
        // (s+ss)^z = s^z * (1+ss/s)^z
        //          = s^z * exp(z*log1p(ss/s))
        // ss/s < machine epsilon : log1p(ss/s) ~ ss/s
        // exp(x) = 1 when x < machine epsilon
        final DD s = DD.ofSum(x, y);
        double r = Math.pow(s.hi(), z);
        final double t = z * s.lo();
        if (Math.abs(t) > EPS * s.hi()) {
            r *= Math.exp(t / s.hi());
        }
        return r;
    }

    /**
     * Standard precision {@code (x+y)^z}.
     *
     * @param x the first operand
     * @param y the second operand
     * @param z the third operand
     * @return the result
     */
    private static double powNp(int x, double y, double z) {
        return Math.pow(x + y, z);
    }

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
     * <li>If the argument \( a \le 0 \) and is an integer, then the result is positive infinity.</li>
     * <li>If the argument \( a \lt 0 \) and \( s \) is not an integer, then the result is nan.</li>
     * <li>If the argument \( a \) is negative infinity, then the result is nan.</li>
     * <li>If either argument is nan, then the result is nan.</li>
     * </ul>
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
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
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
        double xn = a - ca;
        // Intentional float comparison
        if (odd && xn == -HALF) {
            return zetaImp(s, 1 - a);
        }

        // Compute the two terms either side of zero:
        // -1 < xn < 0 < xn + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        // One term is < 0.5: 0.5^-s overflows when s >= 1024.
        DD sum = DD_NAN;
        if (s < MAX_S) {
            final int n = (int) -s;
            DD pn = DD.of(xn);
            DD pp = DD.ONE.add(xn);
            // Avoid overflow issues using scaling
            final long[] expn = {0};
            final long[] expp = {0};
            if (odd) {
                // Compute accurately in extended precision to handle cancellation.
                // The double-double result is +/- 1 ULP (105 bit precision).
                pn = DDMath.pow(pn, n, expn);
                pp = DDMath.pow(pp, n, expp);
            } else {
                // Power terms accurate to at least double precision.
                pn = pn.pow(n, expn);
                pp = pp.pow(n, expp);
            }
            // If re-scaling and addition create infinity we exit with the IEEE result.
            // Note: if one side overflows then it is unlikely the other side will
            // bring it back to finite and we do not check.
            pn = pn.scalb((int) expn[0]);
            pp = pp.scalb((int) expp[0]);
            sum = pn.add(pp);
        }
        // Check for overflow or return the IEEE result
        if (!sum.isFinite()) {
            return odd && xn > -HALF ?
                Double.NEGATIVE_INFINITY :
                Double.POSITIVE_INFINITY;
        }

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2.
        // Adding to the existing finite sum cannot trigger overflow.

        double xp = 2 + xn;
        if (odd) {
            // If odd compute the terms in extended precision.
            // This is impractical for large |a| so the number of terms
            // is limited and precision will be lost for large |a|.
            // The number of terms depends on how close to
            // half-integer and the size of the exponent.
            // The sum continues until the cancellation in opposing terms is in
            // the low part of the double-double sum. Computing the remaining
            // terms in double precision will have cancellation, and the difference
            // of their sums will overlap the high part of the double-double sum.
            // The degree of cancellation is dependent on the number of remaining
            // negative terms as the positive zeta evaluation to infinity will be
            // larger and more accurate. Ideally the remaining double precision sums
            // have missing bits that do not affect the result. In practice this
            // strategy works to compute a result with many bit of precision and
            // avoid a useless result with catastrophic cancellation.
            //
            // sum          |--------|--------|
            // sp1        |--------|
            // -sn1        |------xx|
            // sp1 + sn1          |----xxxx|
            //
            // When a is close to integer this exits very fast otherwise the
            // number of terms can be large.
            // In the extreme this is limited to 5430 terms when 0.5 +/- 2^-40.
            final int n = (int) -s;
            final double x = xn;
            for (int i = 0; xn > a; i++) {
                xn -= 1.0;
                // Note: DDMath here makes no difference as s is small and the
                // standard pow function is accurate to ~100 bits. If s is large
                // then the terms rapidly reduce in magnitude compared to 0.5^-s
                // and trailing imprecise bits in the term do not change the sum.
                final DD pn = DD.of(xn).pow(n);
                final DD pp = DD.ofSum(2 + i, x).pow(n);
                final DD term = pn.add(pp);
                if (Math.abs(term.hi()) < Math.abs(sum.lo())) {
                    // Switch to a double precision tail.
                    // Reset xn to compute as part of the remaining series.
                    xn += 1.0;
                    break;
                }
                sum = sum.add(term);
            }
            // advance positive x by the number of terms computed (x - xn)
            xp = 2 + (x - xn) + x;
        }

        // Compute the remaining terms
        final double sn1 = negativeSeriesSum(a, xn, s);
        final double sp1 = zetaImp(s, xp);

        return sum.add(DD.ofSum(sp1, sn1)).hi();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
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

        // Two methods are used to increase precision
        // - Extended precision summation (low cost)
        // - Extended precision power function using round-off from (a+k) (selectively applied)
        // Use of either reduces maximum error by 1 ulp.
        // Typical RMS error < 0.4 ulp : max error < 2 ulp until the result is sub-normal.

        // Check the extended precision power will make a difference.
        // If a < 1 then the term a^-s dominates the result and a few extra digits
        // of precision from the power function on remaining terms is lost.
        // This will also be false if (a+n) is exact, or the round-off
        // cannot be used when s is small.
        final PowOperator pow =
            a > 1 && Math.abs(s * DD.ofSum(a, N).lo()) >= EPS ?
            HurwitzZeta::extendedPowNp :
            HurwitzZeta::powNp;

        double p = pow.apply(N, a, -s);

        // Initialise sum with the first tail term
        DD sum = DD.of(0.5 * p);
        // S : k in [0, n-1]
        for (int k = N - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum = sum.add(pow.apply(k, a, -s));
        }

        // I
        final double ti = pow.apply(N, a, 1 - s) / (s - 1);

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (t0 > ti) {
            sum = sum.add(ti).add(t0);
        } else {
            sum = sum.add(t0).add(ti);
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        p = pow.apply(N, a, -s - 1);
        // Use to divide by (a+n)^2
        final double apn = pow.apply(N, a, -2);

        // Rising factorial term : (s)_{2k-1}
        // Note: When s is large the loop exits before the rising factorial overflows.
        // (a+n) >= 9 : 9^-340 = 0 : max (2k-1+s) = 339
        // (339)_{2k-1}; k=50 = Pochammer(339, 99) = 1.5e256
        // The rising factorial will not overflow for k <= 50 before (a+n)^-(2k-1+s) is
        // zero. This is within the length of table F.
        double f = s;
        // 2k - 1
        double k2 = 1;
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Note: an extended precision sum here has no effect on the final result.
        double tsum = 0;
        final double stop = sum.hi() * EPS;
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
        return sum.add(tsum).hi();
    }

    /**
     * Calculates the sum of terms of the power series.
     *
     * <pre>
     *      b-1   1
     *   sum     ---
     *      k=a  k^m
     * </pre>
     *
     * <p>Assumes {@code a} and {@code b} are negative and separated by an integer
     * distance; and {@code exponent >= 2} and integer.
     *
     * <p>Large ranges may be evaluated using a difference of zeta functions.
     *
     * @param a First term in the series to calculate (negative non-integer).
     * @param b Last term in the series to calculate, exclusive (negative non-integer).
     * @param s Exponent (positive integer).
     * @return the sum
     */
    private static double negativeSeriesSum(double a, double b, double s) {
        // This can be computed using a difference of zeta functions.
        // A single call to zeta uses ~10 pow operations; use zeta when the
        // sum will use more.
        // Note: The difference incurs cancellation.
        // When s is even the function is called with b in -[1, 0) and the
        // series is strongly converging and no issue occurs.
        // When s is odd it may be called with large |b|. In this case many
        // terms have been computed to handle most of the cancellation in the
        // result. Worst case is b ~ 5430:
        // zeta(3, 5430) = 1.696e-08
        // zeta(3, 5450) = 1.683e-08
        // Cancellation in lost bits = exponent(max(a, b)) - exponent(a-b) = 7
        // The result is sufficient for reasonable double precision.
        if (b - a > MAX_TERMS) {
            final int sign = SpecialMath.isOdd(s) ? -1 : 1;
            final double zb = zetaImp(s, 1 - b);
            final double za = zetaImp(s, 1 - a);
            return sign * (zb - za);
        }

        // Sum terms in ascending order of magnitude
        double sum = 0;
        double x = a;
        while (x < b) {
            sum += Math.pow(x, -s);
            x += 1.0;
        }
        return sum;
    }
}
