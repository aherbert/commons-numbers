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
package org.apache.commons.numbers.examples.jmh.arrays;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Tests for {@link DoublePrecision}.
 */
public class DoublePrecisionTest {
    @Test
    public void testSplitAssumptions() {
        // The multiplier used to split the double value into high and low parts.
        final double scale = (1 << 27) + 1;
        // The upper limit above which a number may overflow during the split into a high part.
        final double limit = 0x1.0p996;
        Assertions.assertTrue(Double.isFinite(limit * scale));
        Assertions.assertTrue(Double.isFinite(-limit * scale));
        // Cannot make the limit the next power up
        Assertions.assertEquals(Double.POSITIVE_INFINITY, limit * 2 * scale);
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, -limit * 2 * scale);
    }
}
