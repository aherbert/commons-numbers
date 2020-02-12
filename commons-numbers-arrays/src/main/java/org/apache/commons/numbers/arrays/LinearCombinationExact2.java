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
package org.apache.commons.numbers.arrays;

/**
 * Computes linear combinations accurately.
 * This class computes the sum of the products of two sequences of numbers
 * <code>a<sub>i</sub> b<sub>i</sub></code> to high accuracy.
 * It does so by using extended precision multiplication and addition algorithms to
 * preserve accuracy and reduce cancellation effects.
 *
 * <p>It is based on the paper by
 * <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
 * Shewchuk (1997): Arbitrary Precision Floating-Point Arithmetic</a>.
 *
 * <p>Due to the methods used to increase precision it is possible that the computation
 * overflows when the standard dot product does not. This may occur when an individual
 * product exceeds the magnitude of {@link Double#MAX_VALUE}{@code  / 2}.
 * In this case the standard precision result will be returned.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Dot_product">Dot product</a>
 */
public final class LinearCombinationExact2 {
    /*
     * Caveat:
     *
     * The code below uses many additions/subtractions that may
     * appear redundant. However, they should NOT be simplified, as they
     * do use IEEE754 floating point arithmetic rounding properties.
     *
     * Algorithms are based on computing the product or sum of two values x and y in
     * extended precision. The standard result is stored using a double (high part z) and
     * the round-off error (or low part zz) is stored in a second double, e.g:
     * x * y = (z, zz); z + zz = x * y
     * x + y = (z, zz); z + zz = x + y
     *
     * To sum multiple (z, zz) results each is added to the current expansion using
     * an algorithm that expands the split extended precision representation of the number
     * (see Shewchuk (1997) Grow-Expansion and Expansion-Sum [1]).
     */

    /**
     * The multiplier used to split the double value into high and low parts. From
     * Dekker (1971): "The constant should be chosen equal to 2^(p - p/2) + 1,
     * where p is the number of binary digits in the mantissa". Here p is 53
     * and the multiplier is {@code 2^27 + 1}.
     */
    private static final double MULTIPLIER = 1.34217729E8;

    /** The upper limit above which a number may overflow during the split into a high part.
     * Assuming the multiplier is above 2^27 and the maximum exponent is 1023 then a safe
     * limit is a value with an exponent of (1023 - 28) = 2^995. */
    private static final double SAFE_UPPER = 0x1.0p995;

    /** The scale to use when down-scaling during a split into a high part.
     * This must be larger than the multiplier and a power of 2 for exact scaling. */
    private static final double DOWN_SCALE = 0x1.0p-30;

    /** The scale to use when re-scaling during a split into a high part.
     * This is the inverse of {@link #DOWN_SCALE}. */
    private static final double UP_SCALE = 0x1.0p30;

    //static int count1, count2;

    /** Private constructor. */
    private LinearCombinationExact2() {
        // intentionally empty.
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a Factors.
     * @param b Factors.
     * @return \( \sum_i a_i b_i \).
     * @throws IllegalArgumentException if the sizes of the arrays are different.
     */
    public static double value(double[] a,
                               double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("Dimension mismatch: " + a.length + " != " + b.length);
        }

        final int len = a.length;

        if (len == 0) {
            return 0;
        }
        if (len == 1) {
            // Revert to scalar multiplication.
            return a[0] * b[0];
        }

        // Expansion e of length m
        // This may grow to contain 'length' split numbers.
        final double[] e = new double[len * 2];
        // working space for the merge of two expansions
        final double[] g = new double[len * 2];

        // Compute all the split products
        // Add pairs of split products using an optimised expansion sum.
        // Store the count and size of each expanded number:
        int n = 0;
        int[] size = new int[(len + 1) / 2];
        for (int i = 1; i < len; i += 2) {
            size[n] = sumProduct(a[i - 1], b[i - 1], a[i], b[i], e, n * 4);
            n++;
        }
        // Final trailing split product for an odd length
        if ((len & 0x1) == 0x1) {
            size[n] = product(a[len - 1], b[len - 1], e, n * 4);
            n++;
        }

        // e contains n expansions of size 4 and maybe a trailing size 2 expansion:
        // ---- ---- ---- ---- --
        // Perform a distillation sum of the expansions by
        // adding pairs until there are no more pairs.
        // A trailing odd expansion is added to the previous one:
        // -------- -------- --  <-- [1] add pairs
        // -------- ----------   <-- [1] add final trailing number
        // ------------------    <-- [2] add pairs

        // Each addition uses a fast-expansion sum which merges e and f in to a
        // single sequence g of increasing magnitude. This is summed and
        // the output is placed back into e:
        // e_m merge f_n = g_(m+n) => e_(m+n)
        // Due to zero elimination the sum of a pair may be shorter.

        // Since we add pairs the start of each number is spaced by powers of 2.
        // The numbers may be shorter so we store the length of each number.
        // Each iteration the count reduces to floor(n / 2) as the trailing unpaired
        // number is added to the previous one.

        int spacing = 4;
        while (n > 1) {
            int j = 0;
            // i is the index of the second part of the pair
            for (int i = 1; i < n; i += 2) {
                // Numbers start at regular spaced intervals
                final int f0 = i * spacing;
                size[j++] = fastExpansionSum(e, f0 - spacing, size[i - 1], f0, size[i], g);
            }
            // Final trailing number
            if ((n & 0x1) == 0x1) {
                // Index of trailing number
                final int i = n - 1;
                // Add to preceding number:
                // - - - - -   <-- start of iteration with a final trailing number
                // -- -- -     <-- currently here adding the final two numbers
                // -- ---
                // If the preceding number is the result of a pair addition the number
                // now occurs at the next power of 2 spacing.
                j--;
                size[j] = fastExpansionSum(e, j * spacing * 2, size[j], i * spacing, size[i], g);
            }
            n /= 2;
            spacing <<= 1;
        }

        // Sum the expansion
        return getSum(a, b, sum(e, size[0]));
    }

    /**
     * Perform the fast expansion sum of e + f, where e and f are in the same array
     * at positions given. The result is written to e and the size of the new expansion
     * returned.
     *
     * @param e Expansion data.
     * @param e0 Start of expansion e.
     * @param m Length of expansion e.
     * @param f0 Start of expansion f.
     * @param n Length of expansion f.
     * @param g Temporary working space
     * @return the size of the result
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 13</a>
     */
//    private static int fastExpansionSum(double[] e, int e0, int m, int f0, int n, double[] g) {
//        // TODO - change so that merge does the zero elimination and returns the length.
//        // remove zero elimination from the fast expansion sum.
//        // Thus zero elimination is only done if the expansion is to be used
//        // in another sum.
//
//        // Combine e and f to a new expansion sorted in increasing order of magnitude.
//        merge(e, e0, e0 + m, f0, f0 + n, g);
//        final int len = m + n;
//        // Fast-two-sum the initial pair
//        int ei = e0;
//        double q = g[1] + g[0];
//        double r = fastSumLow(g[1], g[0], q);
//        if (r != 0) {
//            e[ei++] = r;
//        }
//        // Two-sum the remaining pairs
//        for (int i = 2; i < len; i++) {
//            final double x = q + g[i];
//            r = sumLow(q, g[i], x);
//            q = x;
//            if (r != 0) {
//                e[ei++] = r;
//            }
//        }
//        // if the sum is zero only add this if no others have been added
//        if (q != 0 || ei == e0) {
//            e[ei++] = q;
//        }
//        return ei - e0;
//    }
    private static int fastExpansionSum(double[] e, int e0, int m, int f0, int n, double[] g) {
        // Combine e and f to a new expansion sorted in increasing order of magnitude.
        // The merge performs zero elimination.
        final int len = merge(e, e0, e0 + m, f0, f0 + n, g);
        if (len < 2) {
            System.arraycopy(g, 0, e, e0, len);
            return len;
        }
        // Fast-two-sum the initial pair
        int ei = e0;
        double q = g[1] + g[0];
        e[ei++] = fastSumLow(g[1], g[0], q);
        // Two-sum the remaining pairs
        for (int i = 2; i < len; i++) {
            final double x = q + g[i];
            e[ei++] = sumLow(q, g[i], x);
            q = x;
        }
        e[ei++] = q;
        return ei - e0;
    }
    /**
     * Merge the expansions e and f into g, sorted in order of increasing magnitude.
     * This performs zero elimination to remove zeros from the sequence g.
     *
     * @param e Expansion data.
     * @param e0 Start of expansion e.
     * @param em End of expansion e.
     * @param f0 Start of expansion f.
     * @param fn End of expansion f.
     * @param g Output sequence.
     * @return the length of the output sequence
     */
//    private static void merge(double[] e, int e0, int em, int f0, int fn, double[] g) {
//        int ei = e0;
//        int fi = f0;
//        int gi = 0;
//        // TODO
//        // Make more efficient. This can be done by looking ahead in the lower
//        // sequence until it is not lower anymore. If the expansions do not
//        // contain zeros then this can be done using an exponential search:
//        // https://en.wikipedia.org/wiki/Exponential_search
//        //
//        // The merge flip flops between the two sequences.
//        // Set the low sequence as e or f.
//        // Advance along the low sequence until is it not lower.
//        // Bulk copy the sequence.
//        // Then swap to the other as the low sequence and repeat.
//        // The current high point in the high sequence can be cached as an absolute.
////        while (ei < em && fi < fn) {
////            if (Math.abs(e[ei]) < Math.abs(e[fi])) {
////                g[gi++] = e[ei++];
////            } else {
////                g[gi++] = e[fi++];
////            }
////        }
//        // Assume that the start is within the array bounds.
//        // Cache the absolute value to reduce the number of calls to Math.abs.
//        double ev = Math.abs(e[ei]);
//        double fv = Math.abs(e[fi]);
//        if (ei < em && fi < fn) {
//            // Stop the loop when either counter reaches the length of the expansion.
//            for (;;) {
//                if (ev < fv) {
//                    g[gi++] = e[ei++];
//                    if (ei == em) {
//                        break;
//                    }
//                    ev = Math.abs(e[ei]);
//                } else {
//                    g[gi++] = e[fi++];
//                    if (fi == fn) {
//                        break;
//                    }
//                    fv = Math.abs(e[fi]);
//                }
//            }
//        }
//        if (ei < em) {
//            System.arraycopy(e, ei, g, gi, em - ei);
//        } else {
//            System.arraycopy(e, fi, g, gi, fn - fi);
//        }
//    }
    private static int merge(double[] e, int e0, int em, int f0, int fn, double[] g) {
        int ei = e0;
        int fi = f0;
        int gi = 0;
        // Assume that the start is within the array bounds.
        // Cache the absolute values of each sequence.
        double ev = Math.abs(e[ei]);
        double fv = Math.abs(e[fi]);
        if (ei < em && fi < fn) {
            // Stop the loop when either counter reaches the length of the expansion.
            for (;;) {
                if (ev < fv) {
                    if (ev != 0) {
                        g[gi++] = e[ei];
                    }
                    ei++;
                    if (ei == em) {
                        break;
                    }
                    ev = Math.abs(e[ei]);
                } else {
                    if (fv != 0) {
                        g[gi++] = e[fi];
                    }
                    fi++;
                    if (fi == fn) {
                        break;
                    }
                    fv = Math.abs(e[fi]);
                }
            }
        }
        if (ei < em) {
            while (ei < em) {
                if (e[ei] != 0) {
                    g[gi++] = e[ei];
                }
                ei++;
            }
        } else {
            while (fi < fn) {
                if (e[fi] != 0) {
                    g[gi++] = e[fi];
                }
                fi++;
            }
        }
        return gi;
    }

    /**
     * Compute the sum of the two products and store the result in the expansion at
     * the given index. Interspersed zeros are removed and the length returned.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param e Expansion
     * @param ei Expansion index
     * @return the length of the new expansion
     */
    private static int sumProduct(double a1, double b1, double a2, double b2, double[] e, int ei) {
        // Expansion e
        double e1 = a1 * b1;
        double e0 = productLow(a1, b1, e1);

        // Expansion f
        final double f1 = a2 * b2;
        final double f0 = productLow(a2, b2, f1);

        // Inline an expansion sum to avoid sorting e and f into a sequence g.
        double x;
        double a;
        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        double e2 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);

//        // Store but remove interspersed zeros
//        int en = ei;
//        if (e0 != 0) {
//            e[en++] = e0;
//        }
//        if (e1 != 0) {
//            e[en++] = e1;
//        }
//        if (e2 != 0) {
//            e[en++] = e2;
//        }
//        e[en++] = x;
//        return en - ei;

        // Zero elimination is performed in the merge functon
        e[ei++] = e0;
        e[ei++] = e1;
        e[ei++] = e2;
        e[ei] = x;
        return 4;
    }

    /**
     * Compute the product and store the result in the expansion at
     * the given index. A trailing low zero is removed.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param e Expansion
     * @param ei Expansion index
     * @return the length of the new expansion
     */
    private static int product(double a1, double b1, double[] e, int ei) {
        double e1 = a1 * b1;
        double e0 = productLow(a1, b1, e1);
//        // remove zeros
//        if (e0 != 0) {
//            e[ei] = e0;
//            e[ei + 1] = e1;
//            return 2;
//        }
//        e[ei] = e1;
//        return 1;

        // Zero elimination is performed in the merge functon
        e[ei] = e0;
        e[ei + 1] = e1;
        return 2;
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @return \( a_1 b_1 + a_2 b_2 \)
     *
     * @see #value(double, double, double, double, double, double)
     * @see #value(double, double, double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2) {
        // Expansion e
        double e1 = a1 * b1;
        double e0 = productLow(a1, b1, e1);

        // Expansion f
        double f1 = a2 * b2;
        double f0 = productLow(a2, b2, f1);

        // Inline an expansion sum to avoid sorting e and f into a sequence g.
        // q is the running scalar product, x & a are working variables
        final double q = e1 + f1;
        double x;
        double a;

        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        double e2 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);

        // Final summation
        return getSum(q, e0 + e1 + e2 + x);
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param a3 First factor of the third term.
     * @param b3 Second factor of the third term.
     * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 \)
     *
     * @see #value(double, double, double, double)
     * @see #value(double, double, double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2,
                               double a3, double b3) {
        // Expansion e
        double e1 = a1 * b1;
        double e0 = productLow(a1, b1, e1);

        // Expansion f
        double f1 = a2 * b2;
        double f0 = productLow(a2, b2, f1);

        // Inline an expansion sum to avoid sorting e and f into a sequence g.
        // q is the running scalar product, x & a are working variables
        double q = e1 + f1;
        double x;
        double a;

        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        double e2 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        double e3 = x;

        // Second sum of the expansion f creates e[0-5]
        f1 = a3 * b3;
        f0 = productLow(a3, b3, f1);
        q += f1;

        // TODO - sort f into e
        double e4;
        double e5;

        // Binary search for insertion point f0 into e

        // Binary search for insertion point f1 into e must start at position after f0 insertion

        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        a = x;
        x = e3 + a;
        e3 = sumLow(e3, a, x);
        e4 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        a = x;
        x = e3 + a;
        e3 = sumLow(e3, a, x);
        a = x;
        x = e4 + a;
        e4 = sumLow(e4, a, x);
        e5 = x;

        // Final summation
        return getSum(q, e0 + e1 + e2 + e3 + e4 + e5);
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @param a2 First factor of the second term.
     * @param b2 Second factor of the second term.
     * @param a3 First factor of the third term.
     * @param b3 Second factor of the third term.
     * @param a4 First factor of the fourth term.
     * @param b4 Second factor of the fourth term.
     * @return \( a_1 b_1 + a_2 b_2 + a_3 b_3 + a_4 b_4 \)
     *
     * @see #value(double, double, double, double)
     * @see #value(double, double, double, double, double, double)
     * @see #value(double[], double[])
     */
    public static double value(double a1, double b1,
                               double a2, double b2,
                               double a3, double b3,
                               double a4, double b4) {
        // TODO
        // Create 2 two-products
        // Sort them into one and then fast expansion sum.
        // The sort may be impractical in complexity.

        // Initial product creates the expansion e.
        // This is merged with each product expansion f using a linear expansion.
        double e0;
        double e1;
        double e2;
        double e3;
        double e4;
        double e5;
        double e6;
        double e7;
        double f0;
        double f1;

        // q is the running scalar product, x & a are working variables
        double q;
        double x;
        double a;
        e1 = a1 * b1;
        e0 = productLow(a1, b1, e1);
        q = e1;

        // First sum of the expansion f creates e[0-3]
        f1 = a2 * b2;
        f0 = productLow(a2, b2, f1);
        q += f1;
        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        e2 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        e3 = x;

        // Second sum of the expansion f creates e[0-5]
        f1 = a3 * b3;
        f0 = productLow(a3, b3, f1);
        q += f1;
        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        a = x;
        x = e3 + a;
        e3 = sumLow(e3, a, x);
        e4 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        a = x;
        x = e3 + a;
        e3 = sumLow(e3, a, x);
        a = x;
        x = e4 + a;
        e4 = sumLow(e4, a, x);
        e5 = x;

        // Second sum of the expansion f creates e[0-7]
        f1 = a4 * b4;
        f0 = productLow(a4, b4, f1);
        q += f1;
        // f0 into e
        x = e0 + f0;
        e0 = sumLow(e0, f0, x);
        a = x;
        x = e1 + a;
        e1 = sumLow(e1, a, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        a = x;
        x = e3 + a;
        e3 = sumLow(e3, a, x);
        a = x;
        x = e4 + a;
        e4 = sumLow(e4, a, x);
        a = x;
        x = e5 + a;
        e5 = sumLow(e5, a, x);
        e6 = x;
        // f1 into e
        x = e1 + f1;
        e1 = sumLow(e1, f1, x);
        a = x;
        x = e2 + a;
        e2 = sumLow(e2, a, x);
        a = x;
        x = e3 + a;
        e3 = sumLow(e3, a, x);
        a = x;
        x = e4 + a;
        e4 = sumLow(e4, a, x);
        a = x;
        x = e5 + a;
        e5 = sumLow(e5, a, x);
        a = x;
        x = e6 + a;
        e6 = sumLow(e6, a, x);
        e7 = x;

        // Final summation
        return getSum(q, e0 + e1 + e2 + e3 + e4 + e5 + e6 + e7);
    }

    /**
     * Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) creates
     * a big value from which to derive the two split parts.
     * <pre>
     * c = (2^s + 1) * a
     * a_big = c - a
     * a_hi = c - a_big
     * a_lo = a - a_hi
     * a = a_hi + a_lo
     * </pre>
     *
     * <p>The multiplicand allows a p-bit value to be split into
     * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
     * Combined they have (p􏰔-1) bits of significand but the sign bit of {@code a_lo}
     * contains a bit of information. The constant is chosen so that s is ceil(p/2) where
     * the precision p for a double is 53-bits (1-bit of the mantissa is assumed to be
     * 1 for a non sub-normal number) and s is 27.
     *
     * @param value Value.
     * @return the high part of the value.
     * @see <a href="https://doi.org/10.1007/BF01397083">
     * Dekker (1971) A floating-point technique for extending the available precision</a>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 17</a>
     */
    private static double highPart(double value) {
        // Avoid overflow
        if (value >= SAFE_UPPER || value <= -SAFE_UPPER) {
            // Do scaling.
            final double x = value * DOWN_SCALE;
            final double c = MULTIPLIER * x;
            final double hi = (c - (c - x)) * UP_SCALE;
            if (Double.isInfinite(hi)) {
                // Number is too large.
                // This occurs if value is infinite or close to Double.MAX_VALUE.
                // Note that multiplication by (2^s+1) can cause hi to have an exponent
                // 1 greater than input value which will overflow if the exponent is already +1023.
                // Revert to the raw upper 26 bits of the 53-bit mantissa (including the assumed
                // leading 1 bit). This conversion will result in the low part being a
                // 27-bit significand and the potential loss of bits during addition and
                // multiplication. (Contrast to the Dekker split which creates two 26-bit
                // numbers with a bit moved to the sign of low.)
                // The conversion will maintain Infinite in the high part where the resulting
                // low part (value - high) is NaN.
                return Double.longBitsToDouble(Double.doubleToRawLongBits(value) & ((-1L) << 27));
            }
            return hi;
        }
        // normal conversion
        final double c = MULTIPLIER * value;
        return c - (c - value);
    }

    /**
     * Compute the low part of the double length number {@code (z,zz)} for the exact
     * product of {@code x} and {@code y}. Each factor is split into high and low
     * parts to compute the lower part of the product result. The standard precision
     * product must be provided.
     *
     * <p>Note: This uses the high part as {@code x * y} and not
     * {@code hx * hy + hx * ty + tx * hy} as specified in Dekker's original paper.
     * See Shewchuk (1997) for working examples.
     *
     * <p>Warning: Dekker's split can produce high parts that are larger in magnitude than
     * the input number as the high part is a 26-bit approximation of the number. Thus it is
     * possible that the standard product {@code x * y} does not overflow but the extended
     * precision sub-product {@code hx * hy} does overflow.
     *
     * @param x First factor.
     * @param y Second factor.
     * @param xy Product of the factors (x * y).
     * @return the low part of the product double length number
     * @see <a href="https://doi.org/10.1007/BF01397083">
     * Dekker (1971) A floating-point technique for extending the available precision</a>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    private static double productLow(double x, double y, double xy) {
        // Split the numbers using Dekker's algorithm
        final double hx = highPart(x);
        final double lx = x - hx;

        final double hy = highPart(y);
        final double ly = y - hy;

        // Compute the multiply low part:
        // err1 = xy - hx * hy
        // err2 = err1 - lx * hy
        // err3 = err2 - hx * ly
        // low = lx * ly - err3
        return lx * ly - (((xy - hx * hy) - lx * hy) - hx * ly);
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Dekker's two-sum algorithm. The values are required to be ordered by magnitude
     * {@code |a| >= |b|}. The standard precision sum must be provided.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param sum Sum of the parts (a + b).
     * @return <code>b - (sum - a)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 6</a>
     */
    private static double fastSumLow(double a, double b, double sum) {
        // bVitual = sum - a
        // b - bVirtual == b round-off
        return b - (sum - a);
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Knuth's two-sum algorithm. The values are not required to be ordered by magnitude.
     * The standard precision sum must be provided.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param sum Sum of the parts (a + b).
     * @return <code>(b - (sum - (sum - b))) + (a - (sum - b))</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 7</a>
     */
    private static double sumLow(double a, double b, double sum) {
        final double bVirtual = sum - a;
        // sum - bVirtual == aVirtual.
        // a - aVirtual == a round-off
        // b - bVirtual == b round-off
        return (a - (sum - bVirtual)) + (b - bVirtual);
    }

    /**
     * Grow the expansion. This maintains the increasing non-overlapping expansion
     * by two-summing the new value through the entire expansion from the given start.
     * Ignores newly created zero parts to limit the expansion growth.
     *
     * @param expansion Expansion
     * @param length Expansion size.
     * @param start Start point to begin the merge. To be used for optimised
     * expansion sum.
     * @param value Value to add.
     * @return the new size.
     */
    private static int growExpansion(double[] expansion, int length, int start, double value) {
        double a = value;
        int size = start;
        for (int i = start; i < length; i++) {
            final double b = expansion[i];
            final double x = a + b;
            final double error = sumLow(a, b, x);
            // Store the smaller magnitude if non-zero.
            // This limits the size of the growing expansion.
            if (error != 0) {
                expansion[size++] = error;
            }
            // Carry the larger magnitude up to the next iteration.
            a = x;
        }
        expansion[size++] = a;
        return size;
    }

    /**
     * Sum to the data.
     *
     * @param p Data to sum.
     * @param length Length of the data.
     * @return the sum
     */
    private static double sum(double[] p, int length) {
        double sum = 0;
        for (int i = 0; i < length; i++) {
            sum += p[i];
        }
        return sum;
    }

    /**
     * Gets the final sum. This checks the high precision sum is finite, otherwise
     * returns the standard precision sum for the IEEE754 result.
     *
     * <p>The standard of high precision sum may be non-finite due to input infinite
     * or NaN numbers or overflow in the summation. However the high precision sum
     * can also be non-finite when the standard dot sum is finite. This occurs when
     * the split product had a component high-part that overflowed during
     * computation of the hx * hy partial result. In all cases returning the
     * standard dot sum ensures the IEEE754 result.
     *
     * @param a Factors.
     * @param b Factors.
     * @param hpSum High precision sum.
     * @return the sum
     */
    private static double getSum(double[] a, double[] b, double hpSum) {
        if (!Double.isFinite(hpSum)) {
            // Either we have split infinite numbers, some coefficients were NaNs,
            // or the split product overflowed.
            // Return the naive implementation for the IEEE754 result
            double sum = 0;
            for (int i = 0; i < a.length; i++) {
                sum += a[i] * b[i];
            }
            return sum;
        }
        return hpSum;
    }

    /**
     * Gets the final sum. This checks the high precision sum is finite, otherwise
     * returns the standard precision sum for the IEEE754 result.
     *
     * <p>The standard of high precision sum may be non-finite due to input infinite
     * or NaN numbers or overflow in the summation. However the high precision sum
     * can also be non-finite when the standard dot sum is finite. This occurs when
     * the split product had a component high-part that overflowed during
     * computation of the hx * hy partial result. In all cases returning the
     * standard dot sum ensures the IEEE754 result.
     *
     * @param sum Standard dot sum.
     * @param dotKSum High precision dotK sum.
     * @return the sum
     */
    private static double getSum(double sum, double dotKSum) {
        if (!Double.isFinite(dotKSum)) {
            // Either we have split infinite numbers, some coefficients were NaNs,
            // or the split product overflowed.
            // Return the naive implementation for the IEEE754 result
            return sum;
        }
        return dotKSum;
    }
}
