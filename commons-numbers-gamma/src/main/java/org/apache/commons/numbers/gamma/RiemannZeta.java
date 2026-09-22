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
 * <a href="https://mathworld.wolfram.com/RiemannZetaFunction.html">
 * Riemann zeta</a> function.
 *
 * <p>\[ \zeta(s) = \sum_{k=1}^\infty \frac{1}{k^s} = \frac{1}{1^s} + \frac{1}{2^s} + \frac{1}{3^s} + \cdots \]
 *
 * <p>The function is formally defined for complex variable \( s \) with \( \mathrm{Re}(s) \gt 1 \),
 * and its analytic continuation elsewhere. This implementation uses real-valued \( s \).
 *
 * <p>The reflection formula is used to map \( s \) to the positive domain:
 *
 * <p>\[ \zeta(1 - s) = 2 \sin(\pi \frac{1-s}{2}) (2\pi)^{-s} \Gamma(s) \zeta(s) \]
 *
 * <p>where \( \Gamma(s) \) is the {@link Gamma} function. Negative arguments may be
 * increasingly inaccurate as the magnitude of \( -s \) becomes large.
 *
 * <p>This code has been adapted from:
 * <ul>
 *  <li>The <a href="https://www.boost.org/">Boost</a>
 *      {@code c++} implementation {@code <boost/math/special_functions/zeta.hpp>}.</li>
 * </ul>
 *
 * @see
 * <a href="https://www.boost.org/doc/libs/1_92_0/libs/math/doc/html/math_toolkit/zetas/zeta.html">
 * Boost C++ Riemann Zeta Function</a>
 * @since 1.4
 */
public final class RiemannZeta {

    /** No instances. */
    private RiemannZeta() {}

    /**
     * Computes the value of \( \zeta(s) \).
     *
     * <p>Special cases:
     * <ul>
     * <li>If the argument is 1, then the result is positive infinity.</li>
     * <li>If the argument is a negative even integer, then the result is 0.</li>
     * <li>If the argument is positive infinity, then the result is 1.</li>
     * <li>If the argument is negative infinity, then the result is nan.</li>
     * <li>If the argument is nan, then the result is nan.</li>
     * </ul>
     *
     * @param s Argument.
     * @return \( \zeta(s) \)
     */
    public static double value(double s) {
        return BoostZeta.zeta(s);
    }
}
