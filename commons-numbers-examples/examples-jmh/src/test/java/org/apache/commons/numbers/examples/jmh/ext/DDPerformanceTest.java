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
package org.apache.commons.numbers.examples.jmh.ext;

import java.util.List;
import java.util.stream.Collectors;
import org.apache.commons.numbers.examples.jmh.ext.DDPerformance.DoubleInt;
import org.apache.commons.numbers.examples.jmh.ext.DDPerformance.DoubleIntFunction;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/**
 * Tests the implementations in {@link DDPerformance}.
 */
class DDPerformanceTest {
    @ParameterizedTest
    @CsvSource({
        "300, 300, 0, 1",
        "1000, 500, 0.001, 0.5",
        //"10000, 50, 0.001, 0.15",
        //"100000, 5, 0.001, 0.05",
    })
    void testKSFunction(int samples, int values, double lx, double ux) {
        final List<DoubleIntFunction> funs = DDPerformance.getImplementations()
            .map(DDPerformance.KSMethod::createFunction)
            .collect(Collectors.toList());
        final DoubleInt[] data = DDPerformance.KSData.createData(samples, values, lx, ux);
        final DoubleIntFunction ref = funs.get(0);
//        DD3.reset();
        for (final DoubleInt d : data) {
            final double x = d.getX();
            final int n = d.getN();
            final double expected = ref.apply(x, n);
            for (int i = 1; i < funs.size(); i++) {
                final double actual = funs.get(i).apply(x, n);
                Assertions.assertEquals(expected, actual, () -> String.format("(%s,%s)", x, n));
//                printf("%d, %s : dd2=%s%n", samples, x,
//                    //java.util.Arrays.toString(DD2.getCounters()),
//                    java.util.Arrays.toString(DD3.getCounters())
//                    );
            }
        }
//        printf("%d, %d, %s, %s : dd2=%s%n", samples, values, lx, ux,
//            //java.util.Arrays.toString(DD2.getCounters()),
//            java.util.Arrays.toString(DD3.getCounters())
//            );
    }

    @Test
    void testNaN() {
        // Demonstrate that operations on inf creates NaN.
        // For this reason special handling of single operations to return the IEEE correct result
        // for overflow may not be required. It is possible to include a multiply that is safe
        // against intermediate overflow. But this may never be used. Instead it may be better
        // to document the class as unsuitable for computations that approach +/- inf, and document
        // that the multiply is safe when the exponent is < 996.
        for (DD3 x : new DD3[] {DD3.create(Double.POSITIVE_INFINITY)}) {
            for (DD3 a : new DD3[] {DD3.create(0), DD3.create(1)}) {
                printf("fastAdd(%s, %s) = %s%n", x, a, x.fastAdd(a));
                printf("add(%s, %s) = %s%n", x, a, x.add(a));
                printf("uncheckedMultiply(%s, %s) = %s%n", x, a, x.uncheckedMultiply(a));
                printf("multiply(%s, %s) = %s%n", x, a, x.multiply(a));
            }
            for (double a : new double[] {0, 1}) {
                printf("fastAdd(%s, %s) = %s%n", x, a, x.fastAdd(a));
                printf("add(%s, %s) = %s%n", x, a, x.add(a));
                printf("multiply(%s, %s) = %s%n", x, a, x.multiply(a));
            }
        }
    }

    private static void printf(String format, Object... args) {
        // @CHECKSTYLE: stop regex
        System.out.printf(format, args);
        // @CHECKSTYLE: resume regex
    }
}
