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
 * <a href="https://en.wikipedia.org/wiki/Trigamma_function">Trigamma function</a>.
 *
 * <p>It is the derivative of the {@link Digamma digamma function}:
 *
 * <p>\[ \psi_1(x) = \frac{d^2}{dx^2} (\ln \Gamma(x)) = \frac{d}{dx} \psi(x) \]
 *
 * <p>It may defined as the series sum:
 *
 * <p>\[ \psi_1(x) = \sum_{k=0}^\infty \frac{1}{(x+k)^2} \]
 *
 * <p>This is a special case of the Hurwitz zeta function:
 *
 * <p>\[ \psi_1(x) = \zeta(2, x) \]
 *
 * <p>Negative arguments are handled using a reflection:
 *
 * <p>\[ \psi_1(1-x) = \frac{\pi^2}{\sin^2(\pi x)} - \psi_1(x)\]
 *
 * <p>Arguments in \( (0, 1) \) are shifted using:
 *
 * <p>\[ \psi_1(x) = \psi_1(x + 1) + \frac{1}{x^2} \]
 *
 * <p>For \( x \ge 1 \) the function uses rational approximations for the result.
 *
 * <p>The trigamma function has poles at integer arguments \( x = 0, -1, -2, \dots, -n \).
 * As \( x \) approaches a pole from below or above the value tends to infinity. This
 * function returns the natural limit \( +\infty \) for the pole.
 *
 * <p>This code has been adapted from:
 * <ul>
 *  <li>The <a href="https://www.boost.org/">Boost</a>
 *      {@code c++} implementation {@code <boost/math/special_functions/trigamma.hpp>}</li>
 * </ul>
 *
 * @see Digamma
 * @see
 * <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/sf_gamma/trigamma.html">
 * Boost C++ Trigamma function</a>
 */
public final class Trigamma {
    /** Private constructor. */
    private Trigamma() {
        // intentionally empty.
    }

    /**
     * Computes the trigamma function.
     *
     * <p>Special cases:
     * <ul>
     * <li>If the argument is infinity, then the result is 0.</li>
     * <li>If the argument is zero or a negative integer, then the result is infinity.</li>
     * <li>If the argument is negative infinity, then the result is nan.</li>
     * <li>If the argument is nan, then the result is nan.</li>
     * </ul>
     *
     * <p>Note: Integer arguments {@code x <= 0} are poles. As {@code x} approaches a pole
     * from below or above the value tends to infinity. This function returns the natural
     * limit for the pole.
     *
     * @param x Argument.
     * @return trigamma(x)
     */
    public static double value(double x) {
        return BoostTrigamma.trigamma(x);
    }
}
