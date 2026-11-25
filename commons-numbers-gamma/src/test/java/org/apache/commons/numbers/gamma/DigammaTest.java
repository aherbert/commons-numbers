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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Tests for {@link Digamma}.
 */
class DigammaTest {
    @Test
    void testDigammaLargeArgs() {
        final int ulps = 1;
        // Reference values from WolframAlpha.
        // This used the exact decimal representation.
        // Where the spot values are not representable doubles the
        // reference values are computed using mpmath digamma (1.4.1).
        TestUtils.assertEquals(4.6001618527380874002, Digamma.value(100), ulps);
        TestUtils.assertEquals(3.9019896734278921970, Digamma.value(50), ulps);
        TestUtils.assertEquals(2.9705239922421490509, Digamma.value(20), ulps);
        TestUtils.assertEquals(2.9958363947076465821, Digamma.value(20.5), ulps);
        TestUtils.assertEquals(2.26221435709414812355615936422, Digamma.value(10.1), ulps);
        TestUtils.assertEquals(2.11685881890043800881857132362, Digamma.value(8.8), ulps);
        TestUtils.assertEquals(1.8727843350984671394, Digamma.value(7), ulps);
        TestUtils.assertEquals(0.42278433509846713939, Digamma.value(2), ulps);
        TestUtils.assertEquals(-100.560885457868672415475307272, Digamma.value(0.01), ulps);
        TestUtils.assertEquals(-4.03903989659218953144202580054, Digamma.value(-0.8), ulps);
        TestUtils.assertEquals(4.2003210041401871251493214984, Digamma.value(-6.3), ulps);
        TestUtils.assertEquals(-0.0000311062512303416496575636842388, Digamma.value(1.4616), ulps);
    }

    @Test
    void testDigammaSmallArgs() {
        // Values for negative powers of 10 from 1 to 30 as computed by WebMathematica with 20 digits.
        // For example, to compute trigamma($i) = polygamma(1, $i), use
        //
        // https://www.wolframalpha.com/input?i=polygamma%281%2C+0.1%29
        double[] expected = {-10.423754940411076795, -100.56088545786867450, -1000.5755719318103005,
            -10000.577051183514335, -100000.57719921568107, -1.0000005772140199687e6, -1.0000000577215500408e7,
            -1.0000000057721564845e8, -1.0000000005772156633e9, -1.0000000000577215665e10, -1.0000000000057721566e11,
            -1.0000000000005772157e12, -1.0000000000000577216e13, -1.0000000000000057722e14, -1.0000000000000005772e15, -1e+16,
            -1e+17, -1e+18, -1e+19, -1e+20, -1e+21, -1e+22, -1e+23, -1e+24, -1e+25, -1e+26,
            -1e+27, -1e+28, -1e+29, -1e+30};
        final int ulps = 1;
        for (int n = 1; n < 30; n += 1) {
            final double x = Math.pow(10.0, -n);
            TestUtils.assertEquals(expected[n - 1], Digamma.value(x), ulps,
                null, () -> Double.toString(x));
        }
    }

    @Test
    void testDigammaZero() {
        // There is a branch-cut (discontinuity between +/- infinity) at poles 0, -1, -2, ..., -n.
        // Currently this function returns NaN for the negative integer poles.
        Assertions.assertEquals(Double.NaN, Digamma.value(-1));
        Assertions.assertEquals(Double.NaN, Digamma.value(-101));
        // At zero it would be possible to return +/-infinity using the sign.
        // Presently the function returns NaN at zero.
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, Digamma.value(Double.MIN_VALUE));
        Assertions.assertEquals(Double.NaN, Digamma.value(0.0));
        Assertions.assertEquals(Double.NaN, Digamma.value(-0.0));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, Digamma.value(-Double.MIN_VALUE));
    }

    @Test
    void testDigammaNonRealArgs() {
        // See NUMBERS-219
        Assertions.assertEquals(Double.NaN, Digamma.value(Double.NaN));
        Assertions.assertEquals(Double.POSITIVE_INFINITY, Digamma.value(Double.POSITIVE_INFINITY));
        Assertions.assertEquals(Double.NaN, Digamma.value(Double.NEGATIVE_INFINITY));
    }
}

