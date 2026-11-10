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

import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.numbers.fraction.BigFraction;
import org.junit.jupiter.api.Assertions;

/**
 * Provide Bernoulli numbers. This is a helper class used in testing the zeta function.
 */
final class BernoulliNumbers {
    /** Resource with the Bernoulli numerators. */
    private static final String NUM = "bernoulli_n.txt";
    /** Resource with the Bernoulli denominators. */
    private static final String DENOM = "bernoulli_d.txt";
    /** Bernoulli number B_2n. */
    private static final BigFraction[] B2N;
    /** Bernoulli number B_1 = -0.5. */
    private static final BigFraction B1 = BigFraction.of(-1, 2);

    static {
        final List<BigInteger> n = read(NUM);
        final List<BigInteger> d = read(DENOM);
        final int size = Math.min(n.size(), d.size());
        B2N = new BigFraction[size];
        for (int i = 0; i < size; i++) {
            B2N[i] = BigFraction.of(n.get(i), d.get(i));
        }
    }

    /** No instances. */
    private BernoulliNumbers() {}

    /**
     * Read the Bernoulli resource.
     *
     * @param filename resource name
     * @return Bernoulli numbers.
     */
    private static List<BigInteger> read(String filename) {
        final ArrayList<BigInteger> data = new ArrayList<BigInteger>();
        try (DataReader in = new DataReader(filename)) {
            while (in.next()) {
                try {
                    final BigDecimal n = in.getBigDecimal(1);
                    data.add(n.toBigIntegerExact());
                } catch (final NumberFormatException | ArithmeticException ex) {
                    Assertions.fail("Failed to load data: " + Arrays.toString(in.getFields()), ex);
                }
            }
        } catch (final IOException ex) {
            Assertions.fail("Failed to load data: " + filename, ex);
        }
        return data;
    }

    /**
     * Provide the Bernoulli number B_n.
     *
     * @param n N.
     * @return B_n
     * @throws IllegalArgumentException if {@code n < 0} or {@code n} is even and larger than the
     * supported range using precomputed numbers.
     */
    static BigFraction value(int n) {
        if (isSupported(n)) {
            // Odd numbers above n=1 are zero
            if ((n & 1) == 1) {
                return n == 1 ? B1 : BigFraction.ZERO;
            }
            // Tabulated even numbers
            return B2N[n >>> 1];
        }
        throw new IllegalArgumentException("Unsupported Bernoulli number: " + n);
    }

    /**
     * Check if the Bernoulli number B_n is supported
     *
     * @param n N.
     * @return B_n
     * @returns false if {@code n < 0} or {@code n} is even and larger than the
     * supported range using precomputed numbers.
     */
    static boolean isSupported(int n) {
        return n >= 0 && ((n & 1) == 1 || (n >>> 1) < B2N.length);
    }
}
