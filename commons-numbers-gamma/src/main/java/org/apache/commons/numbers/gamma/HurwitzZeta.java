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

import org.apache.commons.numbers.core.DD;

/**
 * <a href="https://en.wikipedia.org/wiki/Hurwitz_zeta_function">
 * Hurwitz zeta</a> function.
 *
 * <p>\[ \zeta(s, a) = \sum_{k=0}^\infty \frac{1}{(k+a)^s} \]
 *
 * <p>The function is formally defined for complex variable \( s \) with \( \mathrm{Re}(s) \gt 1 \)
 * and real \( a \ne 0, -1, -2, \cdots \). This series is absolutely convergent for the given
 * values of \( s \) and \( a \). Note the special case \( zeta(s, 1) \) is the
 * {@link RiemannZeta Riemann zeta} function. This implementation uses real-valued \( s \gt 1 \).
 *
 * <p>The implementation is performed by spitting the integral into two parts and using
 * the Euler-Maclaurin formula to approximate the second integral \( \mathrm{I} + \mathrm{T} + \mathrm{R} \)
 * with a continuous integral \( \mathrm{I} \), a tail \( \mathrm{T} \), and a residual error term
 * \( \mathrm{R} \) (not computed).
 *
 * <p>\[ \begin{aligned}
 * \zeta(s, a) &amp;= \sum_{k=0}^{N-1} \frac{1}{(k+a)^s} + \sum_{k=N}^\infty \frac{1}{(k+a)^s} = \mathrm{S} + \mathrm{I} + \mathrm{T} + \mathrm{R} \\
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
 * <p>The integrals are well defined when \( a + N \gt 0 \). Negative \( a \) requires the
 * sum \( S \) of a large number of terms \( N \) with potentially very long runtime times.
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
    private static final int N = 8;
    /** Convergence epsilon for the sum of the tail function. This prevents summation
     * of terms that do not affect the final result. */
    private static final double EPS = 0x1.0p-53;
    /** Asymptotic threshold for large {@code a}. Used when {@code a + N} is not exact. */
    private static final double LARGE_A = (1L << 53) - N;

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses {@code 2k!} divided by Bernoulli number {@code B_2k}.
     * Provides M=14 terms. The test suite uses max 9 before convergence.
     */
    private static final double[] F = {
        12.0, // 2! / (1 / 6)
        -720.0, // 4! / (-1 / 30)
        30240.0, // 6! / (1 / 42)
        -1209600.0, // 8! / (-1 / 30)
        4.790016E7, // 10! / (5 / 66)
        -1.8924375803183792E9, // 12! / (-691 / 2730)
        7.47242496E10, // 14! / (7 / 6)
        -2.950130727918164E12, // 16! / (-3617 / 510)
        1.1646782814350067E14, // 18! / (43867 / 798)
        -4.597978722407473E15, // 20! / (-174611 / 330)
        1.81521054019435456E17, // 22! / (854513 / 138)
        -7.1661652561756672E18, // 24! / (-236364091 / 2730)
        2.82908877253043E20, // 26! / (8553103 / 6)
        -1.1168794925000445E22, // 28! / (-23749461029 / 870)
    };

    /** No instances. */
    private HurwitzZeta() {}

    /**
     * Computes the value of \( \zeta(s, a) \).
     *
     * <p>Special cases:
     * <ul>
     * <li>If the argument \( s \) is 1, then the result is positive infinity.</li>
     * <li>If the argument \( s \lt 1 \), then the result is nan.</li>
     * <li>If the argument \( a \le 0 \) and is an integer, then the result is positive infinity.</li>
     * <li>If the argument \( a \le 0 \) and \( s \) is not an integer, then the result is nan.</li>
     * <li>If the argument \( a \) is negative infinity, then the result is nan.</li>
     * <li>If either argument is nan, then the result is nan.</li>
     * </ul>
     *
     * <p><strong>Warning</strong>
     *
     * <p>Negative \( a \) will have increasing runtime as the magnitude of \( a \) increases,
     * with potentially very long runtimes.
     *
     * @param s Argument.
     * @param a Argument.
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
            if (Math.floor(a) == a) {
                // The term 0^-s is infinity
                return Double.POSITIVE_INFINITY;
            }
            if (Math.floor(s) != s) {
                // pow(a, -s) is not defined for negative a when s is non-integer
                return Double.NaN;
            }
            // a < 0 (non-integer) and s is a positive integer.

            // TODO: Test if this is best way to handle negative a.

            // Sum the series until a is positive.
            // Warning: Max terms ~ 2^53.
            // Terms are ascending magnitude.
            // If s is odd then the sum will be negative and the
            // addition of zeta(s, x) has cancellation.
            // Sum in extended precision. Use a standard sum to catch infinity.
            double sum = 0;
            DD ss = DD.ZERO;
            double x = a;
            double t;
            while (x < 0) {
                t = Math.pow(x, -s);
                sum += t;
                ss = ss.add(t);
                x += 1.0;
            }
            // When a is very close to integer the last term can create overflow.
            // Check the extended precision sum is valid or return the IEEE result.
            if (!ss.isFinite()) {
                return sum;
            }
            // Add the remaining series zeta(s, x) for x > 0
            final double z = zetaImp(s, x);
            ss = ss.add(z);
            return ss.isFinite() ? ss.hi() : sum + z;
        }
        // Use the more accurate Riemann zeta function.
        // This is done after domain validation.
        if (a == 1) {
            return RiemannZeta.value(s);
        }
        return zetaImp(s, a);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a >= 1}
     * @return zeta(s, a)
     */
    static double zetaImp(double s, double a) {
        // Asymptotic Behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0, the I term and the first term of T.
        if (a > LARGE_A) {
            return Math.pow(a, 1 - s) / (s - 1) + Math.pow(a, -s) * 0.5;
        }

        final double apn = a + N;
        double p = Math.pow(apn, -s);

        // Initialise sum with the first tail term
        double sum = 0.5 * p;
        // S : k in [0, n-1]
        for (int k = N - 1; k >= 0; k--) {
            // Descending k sums in order of magnitude for increased precision.
            // Prevents early exit for large s when the term (a+k)^-s is below
            // machine epsilon of the ascending series sum.
            sum += Math.pow(a + k, -s);
        }

        // I
        sum += Math.pow(apn, 1 - s) / (s - 1);

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        // When s is large the loop exits before the rising factorial overflows.

        // Rising factorial term : (s)_{2k-1}
        double f = s;
        // 2k - 1
        double k2 = 1;
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        double tsum = 0;
        final double stop = sum * EPS;
        int i;
        for (i = 0; i < F.length; i++) {
            // p = (a+n)^-(2k-1+s)
            p /= apn;
            final double t = f * p / F[i];
            tsum += t;
            if (Math.abs(t) <= stop) {
                break;
            }
            p /= apn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
        }
        return sum + tsum;
    }
}
