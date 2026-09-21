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
 * <a href="https://en.wikipedia.org/wiki/Riemann_zeta_function">
 * Riemann zeta</a> function.
 *
 * <p>\[ \zeta(s) = \sum_{k=1}^\infty \frac{1}{k^s} = \frac{1}{1^s} + \frac{1}{2^s} + \frac{1}{3^s} + \cdots \]
 *
 * <p>The function is formally defined for complex variable {@code s} with {@code Re(s) > 1},
 * and its analytic continuation elsewhere.
 *
 * <p>This implementation uses real-valued {@code s != 1}.
 *
 * <p>The implementation uses the reflection formula when {@code s < 0}.
 *
 * <p>TODO - add reflection formula
 *
 * <p>References
 * <ol>
 * <li><a href="https://en.wikipedia.org/wiki/Riemann_zeta_function">Riemann zeta function (Wikipedia)</a></li>
 * </ol>
 *
 * @since 1.4
 */
public final class RiemannZeta {

    /** No instances. */
    private RiemannZeta() {}

    /**
     * Computes the value of \( \zeta(s) \).
     *
     * <p>Returns positive infinity if {@code s == 1}.
     *
     * @param s Argument.
     * @return \( \zeta(s) \)
     */
    public static double value(double s) {
        return BoostZeta.zeta(s);
    }
}
