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

import java.util.SplittableRandom;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Tests for edge cases using {@link SpecialMath}.
 */
class SpecialMathTest {
    @ParameterizedTest
    @ValueSource(doubles = {1, 3, 13, Integer.MAX_VALUE, (1L << 53) - 1})
    void testIsOddWithOddNumbers(double x) {
        Assertions.assertTrue(x >= 0 && Math.floor(x) == x, "isOdd requires positive integer argument");
        Assertions.assertTrue(SpecialMath.isOdd(x));
    }

    @ParameterizedTest
    @ValueSource(doubles = {0, 2, 42, Integer.MAX_VALUE + 1L, (1L << 53) - 2,
        // Numbers above -2^53 are even
        (1L << 53) + 1, (1L << 53) + 2,
        0x1p54, 1e100, Double.MAX_VALUE, Double.POSITIVE_INFINITY})
    void testIsOddWithEvenNumbers(double x) {
        Assertions.assertTrue(x >= 0 && Math.floor(x) == x, "isOdd requires positive integer argument");
        Assertions.assertFalse(SpecialMath.isOdd(x));
    }

    @Test
    void testAddAsInt() {
        final int[] values = {Integer.MIN_VALUE, -13, -1, 0, 1, 42, Integer.MAX_VALUE};
        for (final int a : values) {
            for (final int b : values) {
                Assertions.assertEquals(toInteger((long) a + b), SpecialMath.addAsInt(a, b));
            }
        }
    }

    @Test
    void testToInt() {
        final long[] longSpecial = {
            Long.MIN_VALUE, Integer.MIN_VALUE, 0, Integer.MAX_VALUE, Long.MAX_VALUE
        };
        for (final long r : longSpecial) {
            Assertions.assertEquals(toInteger(r), SpecialMath.toInt(r));
        }
        for (int b = 0; b < 64; b++) {
            long r = 1L << b;
            Assertions.assertEquals(toInteger(r), SpecialMath.toInt(r));
            // Set the sign bit
            r |= Long.MIN_VALUE;
            Assertions.assertEquals(toInteger(r), SpecialMath.toInt(r));
        }
        final SplittableRandom rng = new SplittableRandom();
        // ints are passed through unchanged
        rng.ints(100).forEach(r -> Assertions.assertEquals(r, SpecialMath.toInt(r)));
        rng.longs(100).forEach(r -> Assertions.assertEquals(toInteger(r), SpecialMath.toInt(r)));
    }

    private static int toInteger(long r) {
        // Clip to the range of an integer
        r = Math.max(Integer.MIN_VALUE, r);
        r = Math.min(Integer.MAX_VALUE, r);
        return (int) r;
    }
}
