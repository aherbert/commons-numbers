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

/**
 * <a href="https://en.wikipedia.org/wiki/Digamma_function">Digamma function</a>.
 *
 * <p>It is defined as the logarithmic derivative of the \( \Gamma \)
 * ({@link Gamma}) function:
 *
 * <p>\[ \psi(x) = \frac{d}{dx}(\ln \Gamma(x)) = \frac{\Gamma^\prime(x)}{\Gamma(x)} \]
 *
 * <p>Negative arguments are handled using a reflection to make \( x \) positive:
 *
 * <p>\[ \psi(1-x) = \psi(x) + \frac{\pi}{\tan(\pi x)} \]
 *
 * <p>For large \( x \) the function uses the asymptotic expansion:
 *
 * <p>\[ \psi(x) = \ln(x) + \frac{1}{2x} - \sum_{n=1}^\infty \frac{B_{2n}}{2nx^{2n}} \]
 *
 * <p>where \( B_{2n} \) is a Bernoulli number. Note that the expansion is divergent after
 * a small number of terms, and the exact point of divergence depends on \( x \). The
 * expansion can be truncated at a term too small to effect the result. Here the limit
 * is chosen as \( x = 10 \) and the summation evaluated as a polynomial in \( 1 / x^2 \).
 *
 * <p>For \( x \lt 10 \) the argument is shifted to \( [1, 2] \)
 * and a rational approximation is used for the region that contains the root.
 * Arguments are shifted using the recurrence relation:
 *
 * <p>\[ \psi(x + 1) = \psi(x) + \frac{1}{x} \]
 *
 * <p>For negative \( x \) the digamma function transitions between positive and negative
 * infinity between each integer value. Poles at integer \( -x \) are undefined.
 * Negative \( x \) will be increasingly inaccurate close to each root between neighbouring
 * integers due to cancellation of positive and negative terms.
 *
 * <p>This code has been adapted from:
 * <ul>
 *  <li>The <a href="https://www.boost.org/">Boost</a>
 *      {@code c++} implementation {@code <boost/math/special_functions/digamma.hpp>}</li>
 * </ul>
 *
 * @see Gamma
 */
public final class Digamma {
    /** Private constructor. */
    private Digamma() {
        // intentionally empty.
    }

    /**
     * Computes the digamma function.
     *
     * <p>Special cases:
     * <ul>
     * <li>If the argument is infinity, then the result is infinity.</li>
     * <li>If the argument is zero or a negative integer, then the result is nan.</li>
     * <li>If the argument is negative infinity, then the result is nan.</li>
     * <li>If the argument is nan, then the result is nan.</li>
     * </ul>
     *
     * <p>Note: Integer arguments {@code x <= 0} are poles. As {@code x} approaches a pole
     * from below or above the value tends to infinity or -infinity. There is a discontinuity
     * at the pole and this function returns nan.
     *
     * @param x Argument.
     * @return digamma(x)
     */
    public static double value(double x) {
        return BoostDigamma.digamma(x);
    }
}
