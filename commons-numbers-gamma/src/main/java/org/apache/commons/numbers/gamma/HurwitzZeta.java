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
 * <p>For positive \( a \) the number of terms \( N \) is chosen dynamically. The
 * terms \( \mathrm{I} + \mathrm{T} \) are written as
 * \( (a+N)^{1-s} \left[ \frac{1}{s-1} + \frac{1}{2(a+N)} + \sum_{k=1}^M c_k \right] \)
 * where the bracket does not require a power function. Starting from \( N = 0 \) the
 * bracket is evaluated if the correction terms \( c_k \) converge within \( M \) terms;
 * otherwise the next term of the sum \( \mathrm{S} \) is added. This minimises the
 * number of power function evaluations. The sum \( \mathrm{S} \) also stops if the
 * remaining terms are negligible.
 *
 * <p>The integrals are well defined when \( a + N \gt 0 \). Negative \( a \) may require the
 * sum \( S \) of a large number of terms \( N \). The implementation uses a difference
 * of zeta evaluations to compute the sum \( S \) over negative
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
 * <p>The function is extremely sensitive to small changes in input for large \( s \) and
 * negative \( a \). For odd \( s \) the function changes rapidly from \( -\infty \) to
 * \( \infty \) somewhere between each negative integer \( a \) close to the half-integer midpoint.
 * For even \( s \) the function is \( \infty \) at each negative integer and may be finite
 * in between. For these reasons results should be treated with caution when \( s \)
 * is large and \( a \) negative.
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
    /**
     * Maximum number of terms N of the series summation S.
     *
     * <p>N is chosen dynamically. At this limit the Euler-Maclaurin tail is used without
     * a convergence requirement; with {@link #M} terms this does not occur on the
     * reference data (the maximum observed N is 7).
     */
    private static final int MAX_N = 9;
    /**
     * Maximum number of correction terms M of the Euler-Maclaurin tail.
     *
     * <p>The correction terms do not require a power function. A larger M allows
     * a smaller N. The error does not increase with the number of terms used.
     */
    private static final int M = 16;
    /** Relative convergence threshold for the Euler-Maclaurin correction terms and the
     * remainder of the series sum. Equal to 2^-56. */
    private static final double CONVERGED = 0x1.0p-56;
    /** Threshold for a negligible first correction term relative to the convergence
     * threshold. The bracket then requires no correction terms and supports any
     * {@code a} without overflow. Equal to 2^-50. */
    private static final double NEGLIGIBLE_C1 = 0x1.0p-50;
    /** Threshold used when {@code a < 1} where {@code a^-s} == zeta(s, a) with the smallest supported
     * {@code s = 1 + 2^-52}. Set using zeta(s, 1) * 2^54. No possible evaluation of
     * all remaining terms for any s could be added to {@code a^-s} when {@code a < 1}.
     *
     * <p>Approximately equal to 2^106 thus remaining terms (all less than 1) cannot be added in
     * double-double precision. */
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
    /** 2. */
    private static final int TWO = 2;
    /**
     * Maximum difference in exponents to support double-double addition.
     * Supports addition of smaller terms to an effectively 107-bit DD number.
     *
     * <p>Note a double-double number (x0,x1) has a floating low part (x1) and the effective
     * mantissa is {@code >= 106} bits. It is only 106 bits if the low part is exactly
     * +/-0.5 ulp(x0). If the low-part is smaller then the exponent of the low
     * part changes to expand the mantissa and a smaller DD number can be added
     * even if their exponents differ by more than 106.
     * <pre>
     * (x0,x1)
     *  x0 1-------
     *  x1         10000000
     *  x1         01-------
     *  x1         001-------
     * + (y0,y1)           1------- --------
     * </pre>
     * <p>The probability of adding a smaller number using an exponent difference of
     * {@code -107 - b} is approximatey {@code 2^-b}. Using b = 0 allows ignoring
     * addition of double-double numbers where the probability of changing the least
     * significant bit of the larger term is less than 0.5, and this is only if the
     * least significant bit is from an effective mantissa of over 107-bits. The
     * chance of changing any bits in the high part is zero. Allows computation
     * of two zeta terms that may have total cancellation in the high 53-bits
     * to have precision in the lower bits approaching the limit of DD arithmetic.
     */
    private static final int MAX_EXP_DELTA = 107;
    /** Downscale to use on the initial rising factorial tail term (always above 1).
     * Equal to 2^-1022 * 2^106 = 2^-916. */
    private static final double DOWNSCALE_916 = 0x1.0p-916;
    /** Downscale to use on the final sum.
     * Equal to 2^-106 */
    private static final double DOWNSCALE_106 = 0x1.0p-106;
    /** Upscale to use on the initial probability tail term (always below 1).
     * Multiplicative inverse of Double.MIN_NORMAL: 2^1022. */
    private static final double UPSCALE_1022 = 0x1.0p1022;
    /** Upscale to use on the sum before scaled tail computation. Allows convergence
     * on a normal number using a sum of normal terms.
     * Multiplicative inverse of DOWNSCALE2: 2^106. */
    private static final double UPSCALE_106 = 0x1.0p106;
    /** Threshold below which a double-double product with a double result is computed
     * using scaling. The low part of the product must not be sub-normal: the result must
     * be above approximately 2^-1074 * 2^106 = 2^-968. Equal to 2^-958 (a margin of 2^10). */
    private static final double DD_PRODUCT_THRESHOLD = 0x1.0p-958;
    /** Threshold for large {@code s} where {@code zeta(s, a) = a^-s} in double precision
     * for all {@code a > 0}.
     *
     * <p>The remainder of the series after the first term is bounded:
     * <pre>
     * sum_{k>=1} (a+k)^-s <= (a+1)^-s * (1 + (a+1) / (s-1))
     * </pre>
     * <p>Relative to {@code a^-s} this is largest for the largest {@code a} with a
     * non-zero result, {@code a^-s >= 2^-1074}: {@code a = 2^(1074/s)}. Smaller {@code a}
     * (including {@code a <= 1}) have a smaller relative remainder. Solving for a
     * remainder below {@code 2^-54 * a^-s} (half an ulp) gives {@code s >= 344.17}.
     * A value of 350 has a remainder below {@code 2^-56.8 * a^-s} (maximum
     * {@code a = 8.39}). If {@code a^-s} underflows to zero ({@code a^-s < 2^-1075})
     * the same bound limits the remainder to below {@code 2^-1131} for {@code s >= 350},
     * so zero is also the correctly rounded result. */
    private static final double SINGLE_TERM_S = 350;
    /** Threshold below which a power term is computed using a split power to avoid
     * sub-normal loss of precision: {@code x^y = x^(y/2) * 2^1022 * x^(y/2)}. The
     * zeta sum is then scaled by 2^1022. Equal to 2^-1022 ({@link Double#MIN_NORMAL}).
     *
     * <p>A single {@code Math.pow} is accurate to 1 ulp but in the sub-normal range the
     * ulp is fixed at 2^-1074: relative precision halves for each binade below 2^-1022.
     * The split power has the error of two pow calls and a multiply but no loss of
     * precision. For a single term the crossover is between 2^-1024 and 2^-1022.
     * However this threshold also decides if the entire sum (later terms, I and the tail)
     * is scaled. A lower threshold leaves smaller terms in sub-normal precision.
     * Testing on the reference data (2^-1022, 2^-1023 and 2^-1024) found 2^-1022
     * equal or better in all cases; lower thresholds increased the max error
     * (e.g. 1.26 to 1.68 ulp on results ~2^-1022). */
    private static final double SCALED_SUM_THRESHOLD = Double.MIN_NORMAL;
    /** Minimum normalised exponent for a double.
     * Note 2^-1075 is the minimum non-zero double.
     * Equivalent to Double.MIN_EXPONENT - 53. */
    private static final int MIN_DOUBLE_EXPONENT = -1075;
    /** The minimum deviation from {@code x=1} to use {@code Math.pow(x, y)} when
     * {@code |y| > SAFE_POW_Y}. Below this threshold a helper function is used to
     * compute {@code (1 + (x-1))^y}. The {@code Math.pow} error occurs for
     * {@code ~1.65e-7 < |x-1| < 2^-20}. */
    private static final double SAFE_DEVIATION_FROM_1 = 0x1.0p-20;
    /** The maximum magnitude of {@code y} to use {@code Math.pow(x, y)} for any {@code x}.
     * Above this threshold a helper function is used to compute {@code x^y} when
     * {@code x} is close to 1. This is the threshold for the large {@code |y|}
     * branch of the fdlibm pow algorithm that contains the error. */
    private static final double SAFE_POW_Y = 0x1.0p31;
    /** High part of ln(2) with 32 significant bits (fdlibm ln2_hi).
     * {@code k * LN2_HI} is exact for {@code |k| < 2^21}. */
    private static final double LN2_HI = 0x1.62e42feep-1;
    /** Low part of ln(2) (fdlibm ln2_lo). {@code LN2_HI + LN2_LO} is ln(2) to ~2^-86. */
    private static final double LN2_LO = 0x1.a39ef35793c76p-33;

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
        // Limit the table size to the maximum number of correction terms.
        // The term M+1 is used to check convergence.
        F = Arrays.stream(FDD).limit(M + 1).mapToDouble(DD::hi).toArray();
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
        // Euler-Maclaurin summation with N=0 if the correction terms converge.
        // This includes the asymptotic behaviour as a -> inf, or tiny s:
        // zeta ~ a^(1-s) / (s-1) + a^-s / 2 = a^(1-s) * (1 / (s-1) + 1 / (2a))
        // https://dlmf.nist.gov/25.11#E43
        // Requires a > 2 so a^(1-s) is not close to 1 (Math.pow is accurate).
        final boolean em0 = a > 2;
        if (em0) {
            final double z = zetaEm0(s, a);
            if (!Double.isNaN(z)) {
                return z;
            }
        }

        // First term. safePow is only required here: after the single-term trap
        // (s < SINGLE_TERM_S) the exponent is small (|y| < 2^31) and calls to Math.pow
        // are accurate for a base close to 1.
        double t = safePow(a, -s);

        // The remaining terms cannot be added to a^-s when:
        // - a^-s is very large (can overflow if 0 <= a < 1), see LARGE_T0
        // - s is large, see SINGLE_TERM_S
        if (t >= LARGE_T0 || s >= SINGLE_TERM_S) {
            return t;
        }

        // Now any (a+k)^-s cannot overflow and the sum cannot overflow.
        // Here a < ~2^66 (zetaEm0 converges for larger a) so a + k is exact in DD.

        // Methods are used to increase precision
        // - Extended precision summation (low cost)
        // - Extended precision power function using round-off from (a+k) (selectively applied)
        // - Computation of sub-normal power terms in parts (low cost)
        // - Scaling for sub-normal terms (no cost)

        // Reciprocal used for the convergence bound (not required to be exact).
        // s > 1 + a * 2^-53 so the division is safe.
        final double rs1 = 1 / (s - 1);

        // Check the extended precision power can be used.
        // If there is no round-off from (a+k) use standard precision.
        // Note: s < SINGLE_TERM_S so the extended precision power is limited to |z| < 350.
        final PowOperator pow = twoSumLow(a, MAX_N, a + MAX_N) == 0 ?
            HurwitzZeta::stdPowxpn :
            HurwitzZeta::extPowxpn;

        // Scaling: if active, all terms and the sum are scaled by 2^1022.
        boolean scaled = false;
        if (t < SCALED_SUM_THRESHOLD) {
            t = scaledPow(pow, a, 0, -s);
            if (t == 0) {
                // (a+k)^(1-s) / (s-1) is not 2^1022 larger than a^-s for any s > 1.
                // The max is ~ 2^53 with a~2^53 and s~1.
                return 0;
            }
            scaled = true;
        }

        // Forward summation of S_k = sum_{i<k} (a+i)^-s. At each k with x = a + k:
        // zeta(s, a) = S_k + x^(1-s) * [1 / (s-1) + 1 / (2x) + sum_{j>=1} c_j]
        // The sum is a normalised double-double (s0, s1) with |s1| <= ulp(s0) / 2.
        // The terms are decreasing so each term can be added using a fast-two-sum.
        // Here a < ~100 (zetaEm0 converges for larger a when s < SINGLE_TERM_S):
        // t_k / t_(k-1) = ((a+k-1) / (a+k))^s < ~0.99 so the computed terms are
        // ordered (s0 >= t_(k-1) > t_k) and the fast-two-sum is exact.
        // This computes the same result as DD.add(double).
        double s0 = 0;
        double s1 = 0;
        for (int k = 0;; k++) {
            if (k != 0) {
                // t = (a+k)^-s
                if (scaled) {
                    t = scaledPow(pow, a, k, -s);
                } else {
                    t = pow.apply(a, k, -s);
                    if (t < SCALED_SUM_THRESHOLD) {
                        // First sub-normal term: activate scaling.
                        // The sum is below ~2^-920 (it has not converged with a sub-normal
                        // term for s > 19.3) so scaling by 2^1022 cannot overflow and is exact.
                        // The scaled sum is normalised.
                        s0 *= UPSCALE_1022;
                        s1 *= UPSCALE_1022;
                        t = scaledPow(pow, a, k, -s);
                        scaled = true;
                    }
                }
            }
            // Euler-Maclaurin tail at x = a + k (k = 0 already tested if a > 2)
            if (k != 0 || !em0) {
                final DD x = DD.ofSum(a, k);
                // Estimate of x^(1-s) (scaled if the sum is scaled) used for convergence.
                final double px = x.multiply(t).hi();
                // The convergence threshold is relative to the total S_k + x^(1-s) / (s-1).
                // At the maximum N the tail is used without a convergence requirement.
                final DD f = emBracket(s, x, s0 / px, k >= MAX_N);
                if (f != null) {
                    // Compute x^(1-s) with a power function. This is more accurate than
                    // the estimate x^-s * x when x^(1-s) / (s-1) dominates (s near 1).
                    final double p;
                    if (!scaled) {
                        p = pow.apply(a, k, 1 - s);
                    } else if (px < SCALED_SUM_THRESHOLD * UPSCALE_1022) {
                        // Sub-normal x^(1-s)
                        p = scaledPow(pow, a, k, 1 - s);
                    } else {
                        p = UPSCALE_1022 * pow.apply(a, k, 1 - s);
                    }
                    if (scaled) {
                        // Scaled result is above 2^-53: the double-double product is exact.
                        return DD.ofSum(s0, s1).add(f.multiply(p)).hi() * Double.MIN_NORMAL;
                    }
                    if (s0 + p * f.hi() < DD_PRODUCT_THRESHOLD) {
                        // Avoid a sub-normal low part of the double-double product
                        return DD.ofSum(s0 * UPSCALE_106, s1 * UPSCALE_106)
                            .add(f.multiply(p * UPSCALE_106)).hi() * DOWNSCALE_106;
                    }
                    return DD.ofSum(s0, s1).add(f.multiply(p)).hi();
                }
            }
            // (s0, s1) += t
            final double q = s0 + t;
            final double e = fastTwoSumLow(s0, t, q) + s1;
            s0 = q + e;
            s1 = fastTwoSumLow(q, e, s0);
            // Convergence: the remainder of the series is bounded by the next term and an
            // integral:  sum_{j>k} (a+j)^-s <= t_{k+1} * (1 + (a+k+1) / (s-1))
            //                               <= t_k * (1 + (a+k+1) / (s-1))
            if (t * (1 + (a + k + 1) * rs1) < CONVERGED * s0) {
                return scaled ? s0 * Double.MIN_NORMAL : s0;
            }
        }
    }

    /**
     * Compute the Hurwitz zeta function using the Euler-Maclaurin sum with N=0 (no
     * series terms):
     *
     * <pre>
     * zeta(s, a) = a^(1-s) * [1 / (s-1) + 1 / (2a) + sum_{k>=1} c_k]
     * c_k = B_2k / (2k)! * (s)_(2k-1) * a^-2k
     * </pre>
     *
     * <p>The bracket does not require a power function. It is evaluated first and the
     * single power {@code a^(1-s)} is computed only if the correction series converges
     * within {@link #M} terms. The remainder after a term is bounded by 2 times the first
     * omitted term. If the first correction term {@code c_1} is negligible this is the
     * asymptotic approximation for large {@code a}, or tiny {@code s}.
     *
     * <p>With corrections the terms {@code 1/(s-1)}, {@code 1/(2a)} and {@code c_1} are
     * summed in double-double precision; the remaining correction terms in double
     * precision. The product with the power is rounded once.
     *
     * <p>Overflow and underflow: after the zero test {@code s < 2^53} ({@code s - 1} is
     * exact) and {@code (s+2k-1)(s+2k)} cannot overflow. If {@code a^2} overflows then
     * {@code g = 0} and the series converges at k=1. The convergence threshold is
     * {@code >= 2^-56 / (s-1)} and one term reduces by at most ~2^-111, so the term
     * that ends the series is normal. The double-double bracket is used only if
     * {@code c_1 >= 2^-50} times the threshold, i.e. {@code a < ~2^66}, so the
     * double-double operations on {@code a} cannot overflow. The terms are decreasing
     * (the series stops if a term increases) and bounded by {@code c_1}.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 2}
     * @return zeta(s, a), or NaN if the correction series did not converge
     */
    private static double zetaEm0(double s, double a) {
        // Zero result. Integral bound: zeta <= a^-s + a^(1-s) / (s-1).
        // a^(1-s) <= 2^(-e * (s-1)) with e = exponent(a) (a >= 2^e, e >= 1)
        // 1 / a <= 2^-e; 1 / (s-1) <= 2^-g with g = exponent(s-1)
        // zeta <= 2^(-e * (s-1)) * 2^(1 - min(e, g))
        // The result is below 2^-1075 (rounds to 0) when e * (s-1) > 1076 - min(e, g).
        // A margin of 2 allows for rounding of the product.
        // This includes all s >= 2^53, and all s >= SINGLE_TERM_S when a > 2^53.
        final int e = Math.getExponent(a);
        if ((s - 1) * e > 1078 - Math.min(e, Math.getExponent(s - 1))) {
            return 0;
        }

        final double a2 = a * a;
        // 1 / (s-1) is a lower bound for the bracket: the corrections sum is positive
        final double stop = CONVERGED / (s - 1);
        // g = (s)_(2k-1) * a^-2k
        double g = s / a2;
        double previous = Double.POSITIVE_INFINITY;
        final double c1 = Math.abs(F[0] * g);
        // Correction terms k >= 2 in double precision
        double sum = 0;
        int k = 1;
        for (;; k++) {
            final double c = Math.abs(F[k - 1] * g);
            if (c < stop) {
                break;
            }
            // Correction terms 1 to k are required.
            // Asymptotic series: once the terms increase they cannot converge.
            if (k > M || c >= previous) {
                return Double.NaN;
            }
            if (k >= 2) {
                sum += F[k - 1] * g;
            }
            previous = c;
            g *= (s + 2 * k - 1) * (s + 2 * k) / a2;
        }

        // s - 1 is exact (s < 2^53).
        final DD f;
        // c_1 is included in the bracket unless negligible. This avoids up to ~0.1 ulp
        // error from omitting a term just below the convergence threshold.
        // A negligible c_1 (a > ~2^33) uses no corrections: a may be very large and the
        // double-double operations on a could overflow.
        if (k == 1 && c1 < stop * NEGLIGIBLE_C1) {
            // No corrections: asymptotic approximation. a may be very large.
            f = DD.fromQuotient(1, s - 1).add(HALF / a);
        } else {
            // Bracket: 1/(s-1) + 1/(2a) + s / (12 a^2) in double-double, then the corrections.
            // Here c_1 >= 2^-50 * stop so a < ~2^66 and the double-double operations are safe.
            // Avoid double-double division: a and s - 1 are exact.
            final DD ra = DD.fromQuotient(1, a);
            f = DD.fromQuotient(1, s - 1)
                .add(ra.multiply(HALF))
                .add(ra.square().multiply(DD.fromQuotient(s, 12)))
                .add(sum);
        }
        // a > 2 so a^(1-s) is not close to 1 and Math.pow is accurate.
        final double p = Math.pow(a, 1 - s);
        if (p < SCALED_SUM_THRESHOLD) {
            // Sub-normal power: limited relative precision. The product would round twice
            // (p, then the result) and can be incorrectly rounded close to a tie.
            // Recompute as a split power with full relative precision:
            // a^(1-s) * 2^1022 = a^((1-s)/2) * 2^1022 * a^((1-s)/2).
            // The scaled power is normal and the scaled result is above ~2^-53 (the bracket
            // is above 1 / (s-1) and s < 2^53) so the double-double product is exact.
            // The final scaling rounds once.
            return f.multiply(scaledPow(HurwitzZeta::stdPowxpn, a, 0, 1 - s)).hi() * Double.MIN_NORMAL;
        }
        // The double-double product is exact only if its low part is not sub-normal.
        // If the result is below ~2^-968 the low part loses precision and the high part
        // can be wrong by up to 1 ulp. Compute using p scaled by 2^106 (exact).
        // The final scaling rounds once if the result is sub-normal.
        if (p * f.hi() < DD_PRODUCT_THRESHOLD) {
            return f.multiply(p * UPSCALE_106).hi() * DOWNSCALE_106;
        }
        return f.multiply(p).hi();
    }

    /**
     * Compute the bracket of the Euler-Maclaurin tail at {@code x}:
     *
     * <pre>
     * 1 / (s-1) + 1 / (2x) + sum_{j>=1} c_j
     * c_j = B_2j / (2j)! * (s)_(2j-1) * x^-2j
     * </pre>
     *
     * <p>The series must converge within {@link #M} terms: {@code |c_j| < 2^-56 *
     * (1/(s-1) + extra)}, where {@code extra} is the current sum relative to
     * {@code x^(1-s)}. The terms {@code 1/(s-1)}, {@code 1/(2x)} and {@code c_1} are
     * summed in double-double precision; the remaining terms in double precision.
     *
     * @param s Argument {@code s > 1} with {@code s < 2^53}
     * @param x Argument {@code x > 0} with {@code x < ~2^66}
     * @param extra Magnitude of the current sum relative to {@code x^(1-s)}.
     * @param force If true return the series truncated before the smallest term (or
     * at M terms) if it does not converge.
     * @return the bracket, or null if the series did not converge
     */
    private static DD emBracket(double s, DD x, double extra, boolean force) {
        final double xd = x.hi();
        final double x2 = xd * xd;
        final double stop = CONVERGED * (1 / (s - 1) + extra);
        // g = (s)_(2j-1) * x^-2j
        double g = s / x2;
        final double c1 = Math.abs(F[0] * g);
        double previous = Double.POSITIVE_INFINITY;
        // Correction terms j >= 2 in double precision
        double sum = 0;
        int j = 1;
        for (;; j++) {
            final double c = Math.abs(F[j - 1] * g);
            if (c < stop) {
                break;
            }
            // Asymptotic series: once the terms increase they cannot converge.
            if (j > M || c >= previous) {
                if (!force) {
                    return null;
                }
                // Truncate: terms 1 to j-1
                break;
            }
            if (j >= 2) {
                sum += F[j - 1] * g;
            }
            previous = c;
            g *= (s + 2 * j - 1) * (s + 2 * j) / x2;
        }
        // Avoid double-double division: s - 1 is exact
        final DD rx = reciprocal(x);
        DD f = DD.fromQuotient(1, s - 1).add(rx.multiply(HALF));
        // c_1 is included unless negligible (or truncated: j == 1 and c_1 >= stop)
        if (j > 1 || (c1 < stop && c1 >= stop * NEGLIGIBLE_C1)) {
            f = f.add(rx.square().multiply(DD.fromQuotient(s, 12))).add(sum);
        }
        return f;
    }

    /**
     * Compute the reciprocal of {@code x = xh + xl} where {@code |xl| <= ulp(xh) / 2}.
     *
     * <p>Uses {@code 1 / (xh + xl) = (1 / xh) * (1 - xl / xh + ...)}: the quotient
     * {@code 1 / xh} in double-double precision with a first order correction for the
     * low part. The omitted terms are below {@code 2^-106} relative. This avoids a
     * double-double division.
     *
     * @param x Value
     * @return 1 / x
     */
    private static DD reciprocal(DD x) {
        final double xh = x.hi();
        final double xl = x.lo();
        final DD r = DD.fromQuotient(1, xh);
        if (xl == 0) {
            return r;
        }
        return DD.ofSum(r.hi(), r.lo() - r.hi() * (xl / xh));
    }

    /**
     * Compute {@code (a+k)^z * 2^1022} using the split power
     * {@code (a+k)^(z/2) * 2^1022 * (a+k)^(z/2)}.
     *
     * @param pow Power function
     * @param a Argument a
     * @param k Argument k
     * @param z Exponent z
     * @return the scaled power
     */
    private static double scaledPow(PowOperator pow, double a, int k, double z) {
        final double p = pow.apply(a, k, z * HALF);
        return p * UPSCALE_1022 * p;
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)} for negative {@code a}.
     *
     * <p>The zeta function for negative a can be computed by summing all terms (a+k)^-s
     * until (a+k) is positive and then adding the zeta evaluation for (a+k). This has issues:
     *
     * <ul>
     * <li>As the magnitude of a increases the number of terms to sum becomes impractical.
     * <li>When s is odd the sum of (a+k)^-s for (a+k) below zero is negative. This results
     * in cancellation when added to the zeta evaluation for (a+k) above zero and loss of
     * precision close to roots of the function.
     * <li>When a is a half-integer and s is odd all the terms for (a+k) below zero
     * match corresponding terms for (a+k) above zero leading to total cancellation of all
     * the negative terms.
     * </ul>
     *
     * <p>This implementation computes the sum of (a+k)^-s below zero and above zero
     * using double-double arithmetic. Addition of double-double sums for the two sides handles
     * most of the cancellation when s is odd.
     *
     * <p>The special case when a is a half-integer and s is odd is handled with a call to
     * evaluate only the remaining terms on the positive side.
     * This is named as the reflection identity.
     *
     * <p>When the number of terms for (a+k) below zero is large this can be computed using
     * a difference of zeta evaluations. This has negligible cancellation as one zeta evaluation
     * will include {@code -1 < (a+k) < 0} and is much larger than the other zeta evaluation.
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
        // Exponent for fractional representation
        final int[] exp = {0};

        // Check case of total cancellation
        final boolean odd = SpecialMath.isOdd(s);
        final double x = a - ca;
        // If the term is half integer and s is even the computation can be recycled.
        // Intentional float comparison.
        final boolean halfInteger = x == -HALF;

        if (odd && halfInteger) {
            // Compute zeta(s, 1 - a)
            if (s <= MAX_REFLECTION_S) {
                // Use extended precision but evaluated with precision for a double result
                final DD r = zetaImp((int) s, DD.ofDifference(1, a), null, exp, Context.DOUBLE);
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
        final DD xp = DD.ofSum(1, x);
        final long[] expn = {0};
        final long[] expp = {0};
        final DD pn = DDMath.pow(xn, -is, expn);
        final DD pp;
        final long delta;
        if (halfInteger) {
            // s must be even so copy the existing positive-valued xn result
            pp = pn;
            expp[0] = expn[0];
            delta = 0;
        } else {
            pp = DDMath.pow(xp, -is, expp);
            delta = expp[0] - expn[0];
            if (Math.abs(delta) > MAX_EXP_DELTA) {
                // Cannot add these terms given the largest is [0.5, 1.0) with ulp 2^-106.
                // All other individual terms are smaller so further DD computation is not possible.
                return delta > 0 ?
                    Math.scalb(pp.hi(), intExponent(expp[0])) :
                    Math.scalb(pn.hi(), intExponent(expn[0]));
            }
        }

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1).
        // This is always < 1.65 (exponent=0); individual terms x^-s < 1 (exponent<0).
        // Exit early if remaining terms cannot be added.
        // The result can be finite even if the rescaled terms are infinite. However
        // we do not support further computation from infinite terms so check
        // if anything can be added to these terms.
        if (Math.max(expn[0], expp[0]) > MAX_EXP_DELTA) {
            // Limit of double-double arithmetic.
            // Add a scaled term avoiding overflow by re-scaling after the sum.
            return Math.scalb(
                pn.add(pp.scalb((int) delta)).hi(),
                intExponent(expn[0])
            );
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) + (-1)^s * [ zeta(s, -x) - zeta(s, 1 - a) ]
        // z + za + zb

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
            if (halfInteger) {
                // Already computed
                za = z;
                eza = ez;
            } else {
                // Take care to change the sign of a provided result for the zeta method.
                za = zetaImp(is, xn.negate(), pn.abs(), exp, c);
                eza = exp[0];
            }
            zb = zetaImp(is, DD.ofDifference(1, a), null, exp, c);
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
            // s >= 2. For any -1 < x < 0 : |x|^-s > zeta(s, |x|+1)
            // No sort required for the summation by magnitude
        }

        // Sum in magnitude order
        // zb < za : check z against za
        if (ez > eza) {
            // za,zb < z : za = za + zb
            za = add(za, eza, zb, ezb, exp);
            eza = exp[0];
        } else {
            // z,zb < za : z = z + zb
            z = add(z, ez, zb, ezb, exp);
            ez = exp[0];
        }
        // z = za + z
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
     * <p>Support for sub-normal double values is provided by returning the
     * result using a fractional representation in {@code [0.5, 1)} with
     * a base-2 exponent.
     *
     * <p>The first term can be provided. This must be the fractional representation
     * with the exponent in {@code exp}.
     *
     * @param s Argument {@code 1 < s <= 1837}
     * @param a Argument {@code a > 0}
     * @param a0 {@code a^-s} (can be null)
     * @param exp base-2 exponent of the result.
     * @param c Evaluation context.
     * @return zeta(s, a) fractional part in [0.5, 1)
     * @see DD#frexp(int[])
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
        // Addition of fractional representations is handled by a special routine
        // that can detect addition of zero. Use of scaled addition is not far slower
        // that standard addition as the result from the power function has to be scaled
        // for use in standard addition. So the amount of scaling is the same on each term
        // and the extra expense is re-scaling the running sum.

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

        // Power function for main terms based on the context
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

        // Edge case handling:
        //
        // Ignore asymptotic behavior as a -> inf.
        // When a is large the series cannot use a+k. Occurs when a > 2^106.
        // This function is never called with large a. The reflection case for
        // half-integer a requires a < 2^53.
        //
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
        // Note k=0 already computed as t0
        for (int k = n - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            final DD r = a.add(k).pow(-s, e);
            sum = add(sum, es, r, intExponent(e[0]), exp);
            es = exp[0];
        }

        // I : (a+p)^(1-s) / (s-1)
        DD ti = pow.apply(apn, 1 - s, e);
        long eti = e[0];
        if (s > TWO) {
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
        final long eapn = e[0];

        // Set the scale for the tail series

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Use the tail precision to obtain a minimum exponent for terms.
        final long stop = (long) es - c.getPrecision();

        // s+2k-1
        long s2k = s + 1L;
        // Initialise with the first term.
        // Avoids addition of fractional representation to zero (which has a zero exponent).
        DD tsum = f.multiply(p).multiply(FDD[0]).frexp(exp);
        long et = ef + ep + exp[0];
        int etsum = intExponent(et);
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
            // f*p in [0.25, 2); FDD*2 is finite.
            final DD t = f.multiply(p).multiply(FDD[i]).frexp(exp);
            et = ef + ep + exp[0];
            tsum = add(tsum, etsum, t, intExponent(et), exp);
            etsum = exp[0];
        }
        return add(sum, es, tsum, etsum, exp);
    }

    /**
     * Helper function to add the fractional representation of two DD numbers.
     *
     * <p>Handles the special case of adding {@link DD#ZERO} which has a zero exponent.
     *
     * @param a Addend a.
     * @param ae Exponent of a.
     * @param b Addend b.
     * @param be Exponent of b.
     * @param e Exponent of the result.
     * @return the fraction part of the result
     * @see DD#frexp(int[])
     */
    static DD add(DD a, int ae, DD b, int be, int[] e) {
        // Add if the exponent overlap changes the larger 106-bit term
        final int diff = intExponent((long) ae - be);
        if (diff > MAX_EXP_DELTA) {
            // return a, or b if a is zero
            if (a.isZero()) {
                e[0] = be;
                return b;
            }
            e[0] = ae;
            return a;
        } else if (diff < -MAX_EXP_DELTA) {
            // return b, or a if b is zero
            if (b.isZero()) {
                e[0] = ae;
                return a;
            }
            e[0] = be;
            return b;
        }
        final DD r = b.add(a.scalb(diff)).frexp(e);
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
        return SpecialMath.addAsInt(a, b);
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
        return SpecialMath.toInt(e);
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
     * the infinite series {@code (x+k)^-s}. It is limited to the expected arguments
     * of the zeta function use case:
     *
     * <ul>
     * <li>{@code x > 0}
     * <li>{@code n >= 0}
     * <li>{@code x+n} finite
     * <li>{@code -350 < z < 0}
     * </ul>
     *
     * <p>The exponent is {@code -s}, {@code 1-s}, or half of these for a scaled power.
     * It is limited by the single-term trap for large {@code s} ({@link #SINGLE_TERM_S}).
     *
     * <p>This method has no advantage over {@code Math.pow(x+n, z)} when there is no
     * round-off component to (x+n). In the event that {@code (x+n)} is exact this method
     * will return the same result as {@code Math.pow(x+n, z)}.
     *
     * <p>Calling with {@code n=0} will return {@code Math.pow(x, z)}. More generally
     * this method should not be called when {@code x + n} does not create a change of
     * exponent between {@code x} and {@code x + n}, and may offer no benefit even if the
     * exponent does change depending on the arguments. For example if {@code |x| < 2^52}
     * and is an integer then addition of any {@code int} value will be exact.
     *
     * @param x Argument x ({@code 0 < x < 2^53})
     * @param n Argument n ({@code >= 0})
     * @param z Exponent z ({@code -350 < z < 0})
     * @return the result
     */
    // package-private for testing
    static double extPowxpn(double x, int n, double z) {
        // (s+ss)^z = s^z * (1+t)^z with t = ss/s and |t| <= 2^-53.
        //
        // Taylor series:
        // (1 + t)^z = 1 + zt + z(z-1)t^2/2! + z(z-1)(z-2)t^3/3! + ...
        // for |t| < 1 and real z
        //
        // For |z| < 350 the second term is below (350 * 2^-53)^2 / 2 < 2^-90 relative
        // to the result and is ignored. The sign of ss is not changed: the
        // correction can be in either direction.

        // Compute the exact double-double sum.
        // Avoid generation of a DD intermediate by directly using twoSum.
        final double s = x + n;
        final double ss = twoSumLow(x, n, s);
        // safePow is not required: the Math.pow error close to x=1 (see safePow)
        // only occurs for |z| > 2^31. Here |z| < 350, e.g. (a+1)^-s for a -> 0.
        final double r = Math.pow(s, z);
        return r + r * (z * ss / s);
    }

    /**
     * Compute the round-off of the sum of two numbers {@code a} and {@code b} using
     * Knuth two-sum algorithm. The values are not required to be ordered by magnitude,
     * i.e. the result is commutative {@code s = a + b == b + a}.
     *
     * <p>Copied from {@code o.a.c.numbers.core.DD}.
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

    /**
     * Compute the round-off of the sum of two numbers {@code a} and {@code b} using
     * Dekker's two-sum algorithm. The values are required to be ordered by magnitude:
     * {@code |a| >= |b|}.
     *
     * <p>If {@code a} is zero and {@code b} is non-zero the returned value is zero.
     *
     * <p>Copied from {@code o.a.c.numbers.core.DD}.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return the sum round-off
     */
    private static double fastTwoSumLow(double a, double b, double x) {
        // (x, xx) = a + b
        // bVirtual = x - a
        // xx = b - bVirtual
        return b - (x - a);
    }

    /**
     * Safe single term power function {@code x^y}.
     *
     * <p>This method is used to compute the power function accurately when
     * {@code x -> 1} and {@code |y|} is large. It is limited to the expected arguments
     * of the zeta function use case:
     *
     * <ul>
     * <li>{@code x > 0} and finite
     * <li>{@code y < 0}
     * </ul>
     *
     * <p>This method has no advantage over {@code Math.pow(x, y)} when {@code x} is not
     * close to 1 or {@code |y| <= 2^31}. The method exists to address
     * <a href="https://bugs.openjdk.org/browse/JDK-8326547">JDK-8326547</a> where
     * {@link Math#pow(double, double)} can be outside the 1 ULP specification for
     * {@code |y| > 2^31} and {@code ~1.65e-7 < |1-x| < 2^-20}. The deviation can be
     * over 700 ULP. The error is in the fdlibm pow algorithm used by
     * {@link StrictMath#pow(double, double)} on all platforms, and by
     * {@code Math.pow} on platforms with no pow intrinsic (e.g. AArch64), or on
     * x86_64 in Java 8 for sub-normal results. This bug has been fixed in JDK 27 but no
     * backport exists at the time of development (Nov 2026). This library supports
     * Java 8 and must return an accurate zeta result for the entire domain
     * {@code a > 0}.
     *
     * <p>Use of this helper function will be removed when backports exist on the minimum
     * target Java version.
     *
     * @param x Argument x
     * @param y Argument y
     * @return the result
     * @see <a href="https://bugs.openjdk.org/browse/JDK-8326547">JDK-8326547:
     * StrictMath.pow returns results outside of a 1 ulp error bound</a>
     */
    // package-private for testing
    static double safePow(double x, double y) {
        // Exact for x in [0.5, 2]: |eps| is either 0 or >= 2^-53.
        final double eps = x - 1;
        // Note: x == 1 (eps == 0) uses Math.pow which is exact:
        // 1 for any finite y; NaN for infinite y.
        // NaN y uses Math.pow.
        if (Math.abs(y) > SAFE_POW_Y && Math.abs(eps) < SAFE_DEVIATION_FROM_1 && eps != 0) {
            return pow1p(eps, y);
        }
        return Math.pow(x, y);
    }

    /**
     * Platform independent implementation of {@code (1+x)^y}.
     *
     * <p>This method is accurate as {@code x -> 0} for large {@code |y|}: the error
     * bound is ~1.5 ULP (assuming {@link Math#exp(double)} is accurate to 1 ULP); the
     * measured maximum error is ~1.25 ULP. It is limited to the expected arguments of the
     * zeta function use case:
     *
     * <ul>
     * <li>{@code x = a - 1} for a double {@code a}: {@code x} is exact and
     * {@code 2^-53 <= |x| < 2^-20}. Either sign of {@code x} is supported.
     * <li>{@code y} not NaN. Either sign of {@code y} is supported. Infinite {@code y}
     * returns the limit 0 or infinity.
     * </ul>
     *
     * <p><strong>Warning:</strong> This method does not expect to be called with
     * {@code x=0} where the result for {@code Math.pow(1, y)} is 1 for any finite
     * argument {@code y} and NaN otherwise. This argument is filtered in
     * {@link #safePow(double, double)}.
     *
     * <p>Method: computes {@code z = y * ln(1 + x)} as a double-double {@code (zh + zl)}
     * and returns {@code exp(zh) * (1 + zl)}. The exponent must be carried to ~70 bits
     * because the absolute error in {@code z} becomes the relative error in the result,
     * and {@code |z|} can be up to ~745 (underflow) or ~710 (overflow); a double
     * {@code z} loses log2(745) ~ 9.5 bits. Note: {@code exp(y * log1p(x))} has errors
     * of 100-1000 ULP.
     *
     * @param x Argument x
     * @param y Argument y
     * @return the result (1+x)^y
     */
    private static double pow1p(double x, double y) {
        // z = y * ln(1+x) = p + t with p = y * x.
        // |t| <= |p| * |x| / 2 + ulp(p) < 4e-4 here, so
        // p < -746 gives z < -745.9 (result < 2^-1075 = e^-745.13: rounds to 0) and
        // p > 710 gives z > 709.9 (result > MAX_VALUE = e^709.78: overflows).
        // Every |y| > 2^996 exits here as |x| >= 2^-53 makes |p| >= 2^943.
        final double p = y * x;
        if (p < -746) {
            return 0;
        }
        if (p > 710) {
            return Double.POSITIVE_INFINITY;
        }
        // Here |y| <= 746 / |x| <= 746 * 2^53 < 2^63 so the split of y in the
        // exact product cannot overflow.

        // ln(1+x) = x + lo, lo = -x^2/2 + x^3/3 - x^4/4 + x^5/5
        // Plain Horner evaluation is sufficient. lo only enters z as y*lo, with
        // |y*lo| <= |p|*|x|/2 <= 746 * 2^-21 ~ 3.6e-4, so a few ulp of rounding error
        // in lo is ~8e-20 absolute in z: < 0.001 ulp in the result.
        // Term contributions to z are p * x^(k-1) / k (|p| <= 746, |x| < 2^-20):
        //   x^4/4: up to ~1.6e-16 (~1.5 ulp) - required.
        //   x^5/5: up to ~1.2e-22 (~1e-6 ulp) - kept as a margin.
        //   x^6/6 (truncated): up to ~1e-28 absolute in z, < 2^-102 relative to ln(1+x).
        final double lo = x * x * (-0.5 + x * (1.0 / 3 - x * (0.25 - x * 0.2)));
        // Exact round-off of y*x. Exact unless |p| < ~2^-968 where any error
        // is below 2^-1000 absolute in z and irrelevant to exp(z).
        final double pe = DD.ofProduct(y, x).lo();
        final double t = pe + y * lo;
        double zh = p + t;
        // fast two-sum, |p| >= |t|
        double zl = t - (zh - p);
        if (zh < -746) {
            // result < 2^-1075: underflows to 0
            return 0;
        }
        if (zh > 710) {
            // result > MAX_VALUE: overflows
            return Double.POSITIVE_INFINITY;
        }
        // Scale into [-700, 700] so exp(zh) neither overflows nor is sub-normal;
        // the final scalb then rounds once.
        int k = 0;
        if (zh < -700) {
            k = (int) Math.ceil((-700 - zh) / LN2_HI);
        } else if (zh > 700) {
            k = -(int) Math.ceil((zh - 700) / LN2_HI);
        }
        if (k != 0) {
            // (zh + zl) + k * ln2, renormalised so |zl| <= ulp(zh)/2.
            // zh + k * LN2_HI is exact: k * LN2_HI (|k| <= 68) is a multiple of 2^-32
            // with at most 39 significant bits; zh (|zh| <= 746) is a multiple of 2^-43;
            // and the sum (|sum| <= 700) is representable on that grid.
            final double h = zh + k * LN2_HI;
            final double l = zl + k * LN2_LO;
            zh = h + l;
            zl = l - (zh - h);
        }
        double r = Math.exp(zh);
        // exp(zh + zl) ~ exp(zh) * (1 + zl), |zl| < ~2^-43
        r += r * zl;
        return k == 0 ? r : Math.scalb(r, -k);
    }
}
