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

import java.util.Arrays;

/**
 * Support class for sorting arrays.
 *
 * <p>Optimal sorting networks are used for small fixed size array sorting.
 *
 * <p>Note: Requires that the floating-point data contains no NaN values; sorting
 * does not respect the order of signed zeros imposed by {@link Double#compare(double, double)}.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Sorting_network">Sorting network (Wikipedia)</a>
 * @see <a href="https://bertdobbelaere.github.io/sorting_networks.html">Sorting Networks (Bert Dobbelaere)</a>
 *
 * @since 1.2
 */
final class Sorting {

    /** No instances. */
    private Sorting() {}

    /**
     * Sorts an array using an insertion sort.
     *
     * @param x Data array.
     * @param left Lower bound (inclusive).
     * @param right Upper bound (inclusive).
     */
    static void sort(double[] x, int left, int right) {
        for (int i = left; ++i <= right;) {
            final double v = x[i];
            // Move preceding higher elements above (if required)
            if (v < x[i - 1]) {
                int j = i;
                while (--j >= left && v < x[j]) {
                    x[j + 1] = x[j];
                }
                x[j + 1] = v;
            }
        }
    }

    /**
     * Sorts the given indices in an array.
     *
     * <p>Note: Requires that the range contains no NaN values. It does not respect the
     * order of signed zeros.
     *
     * <p>Assumes all indices are valid and distinct.
     *
     * <p>Data are arranged such that:
     * <pre>{@code
     * a != b != c
     * data[a] < data[b] < data[c]
     * }</pre>
     *
     * <p>If indices are duplicated elements will <em>not</em> be correctly ordered.
     * However in this case data will contain the same values and may be partially ordered.
     *
     * @param x Data array.
     * @param a Index.
     * @param b Index.
     * @param c Index.
     */
    static void sort3(double[] x, int a, int b, int c) {
        // Decision tree avoiding swaps:
        // Order [(0,2)]
        // Move point 1 above point 2 or below point 0
        final double u = x[a];
        final double v = x[b];
        final double w = x[c];
        if (w < u) {
            if (v < w) {
                x[a] = v;
                x[b] = w;
                x[c] = u;
                return;
            }
            if (u < v) {
                x[a] = w;
                x[b] = u;
                x[c] = v;
                return;
            }
            // z < y < z
            x[a] = w;
            x[c] = u;
            return;
        }
        if (v < u) {
            // y < x < z
            x[a] = v;
            x[b] = u;
            return;
        }
        if (w < v) {
            // x < z < y
            x[b] = w;
            x[c] = v;
        }
        // x < y < z
    }

    /**
     * Sorts the given indices in an array.
     *
     * <p>Note: Requires that the range contains no NaN values. It does not respect the
     * order of signed zeros.
     *
     * <p>Assumes all indices are valid and distinct.
     *
     * <p>Data are arranged such that:
     * <pre>{@code
     * a != b != c != d != e
     * data[a] < data[b] < data[c] < data[d] < data[e]
     * }</pre>
     *
     * <p>If indices are duplicated elements will <em>not</em> be correctly ordered.
     * However in this case data will contain the same values and may be partially ordered.
     *
     * @param x Data array.
     * @param a Index.
     * @param b Index.
     * @param c Index.
     * @param d Index.
     * @param e Index.
     */
    static void sort5(double[] x, int a, int b, int c, int d, int e) {
        // Uses an optimal sorting network from Knuth's Art of Computer Programming.
        // 9 comparisons.
        // Order pairs:
        // [(0,3),(1,4)]
        // [(0,2),(1,3)]
        // [(0,1),(2,4)]
        // [(1,2),(3,4)]
        // [(2,3)]
        if (x[e] < x[b]) {
            final double u = x[e];
            x[e] = x[b];
            x[b] = u;
        }
        if (x[d] < x[a]) {
            final double v = x[d];
            x[d] = x[a];
            x[a] = v;
        }

        if (x[d] < x[b]) {
            final double u = x[d];
            x[d] = x[b];
            x[b] = u;
        }
        if (x[c] < x[a]) {
            final double v = x[c];
            x[c] = x[a];
            x[a] = v;
        }

        if (x[e] < x[c]) {
            final double u = x[e];
            x[e] = x[c];
            x[c] = u;
        }
        if (x[b] < x[a]) {
            final double v = x[b];
            x[b] = x[a];
            x[a] = v;
        }

        if (x[e] < x[d]) {
            final double u = x[e];
            x[e] = x[d];
            x[d] = u;
        }
        if (x[c] < x[b]) {
            final double v = x[c];
            x[c] = x[b];
            x[b] = v;
        }

        if (x[d] < x[c]) {
            final double u = x[d];
            x[d] = x[c];
            x[c] = u;
        }
    }

    /**
     * Place the lower median of 4 elements in {@code b}; the smaller element in
     * {@code a}; and the larger two elements in {@code c, d}.
     *
     * @param x Values
     * @param a Index.
     * @param b Index.
     * @param c Index.
     * @param d Index.
     */
    static void lowerMedian4(double[] x, int a, int b, int c, int d) {
        // 4 comparisons
        if (x[d] < x[b]) {
            final double u = x[d];
            x[d] = x[b];
            x[b] = u;
        }
        if (x[c] < x[a]) {
            final double v = x[c];
            x[c] = x[a];
            x[a] = v;
        }
        // a--c
        // b--d
        if (x[b] < x[a]) {
            final double v = x[a];
            final double u = x[c];
            x[a] = x[b];
            x[c] = x[d];
            x[b] = v;
            x[d] = u;
        }
        // a--c
        //    b--d
        if (x[c] < x[b]) {
            final double u = x[c];
            x[c] = x[b];
            x[b] = u;
        }
    }

    /**
     * Place the lower median of 4 elements in {@code b}; the smaller element in
     * {@code a}; and the larger two elements in {@code c, d}.
     *
     * @param x Values
     * @param a Index.
     * @param b Index.
     * @param c Index.
     * @param d Index.
     */
    static void lowerMedian4d(double[] x, int a, int b, int c, int d) {
        // 4 comparisons
        if (x[d] < x[a]) {
            final double u = x[d];
            x[d] = x[a];
            x[a] = u;
        }
        if (x[c] < x[b]) {
            final double v = x[c];
            x[c] = x[b];
            x[b] = v;
        }
        // a--d
        // b--c
        if (x[b] < x[a]) {
            final double xb = x[a];
            x[a] = x[b];
            x[b] = xb;
            //    b--d
            // a--c
            if (x[c] < xb) {
                x[b] = x[c];
                x[c] = xb;
                // fully sorted here
            }
            // else full sort requires c:d ordering
            // Not fully sorted for 6 of 24 permutations
        } else if (x[d] < x[b]) {
            // a--d
            //       b--c
            final double v = x[d];
            // Do a full sort for 1 additional swap
            x[d] = x[c];
            x[c] = x[b];
            x[b] = v;
            // minimum swaps to put the lower median at b
            //x[d] = x[b];
            //x[b] = v;
        }
    }

    /**
     * Place the upper median of 4 elements in {@code c}; the smaller two elements in
     * {@code a,b}; and the larger element in {@code d}.
     *
     * @param x Values
     * @param a Index.
     * @param b Index.
     * @param c Index.
     * @param d Index.
     */
    static void upperMedian4(double[] x, int a, int b, int c, int d) {
        // 4 comparisons
        if (x[d] < x[b]) {
            final double u = x[d];
            x[d] = x[b];
            x[b] = u;
        }
        if (x[c] < x[a]) {
            final double v = x[c];
            x[c] = x[a];
            x[a] = v;
        }
        // a--c
        // b--d
        if (x[d] < x[c]) {
            final double v = x[a];
            final double u = x[c];
            x[a] = x[b];
            x[c] = x[d];
            x[b] = v;
            x[d] = u;
        }
        // a--c
        //    b--d
        if (x[c] < x[b]) {
            final double u = x[c];
            x[c] = x[b];
            x[b] = u;
        }
    }

    /**
     * Place the upper median of 4 elements in {@code c}; the smaller two elements in
     * {@code a,b}; and the larger element in {@code d}.
     *
     * @param x Values
     * @param a Index.
     * @param b Index.
     * @param c Index.
     * @param d Index.
     */
    static void upperMedian4d(double[] x, int a, int b, int c, int d) {
        // 4 comparisons
        if (x[d] < x[a]) {
            final double u = x[d];
            x[d] = x[a];
            x[a] = u;
        }
        if (x[c] < x[b]) {
            final double v = x[c];
            x[c] = x[b];
            x[b] = v;
        }
        // a--d
        // b--c
        if (x[d] < x[c]) {
            final double xc = x[d];
            x[d] = x[c];
            x[c] = xc;
            // a--c
            //    b--d
            if (xc < x[b]) {
                x[c] = x[b];
                x[b] = xc;
                // fully sorted here
            }
            // else full sort requires a:b ordering
            // Not fully sorted for 6 of 24 permutations
        } else if (x[c] < x[a]) {
            //       a--d
            // b--c
            final double v = x[a];
            // Do a full sort for 1 additional swap
            x[a] = x[b];
            x[b] = x[c];
            x[c] = v;
            // minimum swaps to put the lower median at b
            //x[d] = x[b];
            //x[b] = v;
        }
    }

    /**
     * Sort the unique indices in-place to the start of the array. The number of
     * unique indices is returned.
     *
     * <p>Uses an insertion sort modified to ignore duplicates. Use on small {@code n}.
     *
     * <p>Warning: Requires {@code n > 0}. The array contents after the count of unique
     * indices {@code c} is unchanged (i.e. {@code [c, n)}. This may change the count of
     * each unique index in the entire array.
     *
     * @param data Indices.
     * @param n Number of indices.
     * @return the number of indices
     */
    static int insertionSortIndices(int[] data, int n) {
        int unique = 1;
        // Do an insertion sort but only compare the current set of unique values.
        for (int i = 0; ++i < n;) {
            final int v = data[i];
            int j = unique - 1;
            if (v > data[j]) {
                // Insert at end
                data[unique] = v;
                unique++;
            } else if (v < data[j]) {
                // Find insertion point in the unique indices
                do {
                    --j;
                } while (j >= 0 && v < data[j]);
                // Insertion point = j + 1
                // Insert if at start or non-duplicate
                if (j < 0 || v != data[j]) {
                    // Move (j, unique) to (j+1, unique+1)
                    for (int k = unique; --k > j;) {
                        data[k + 1] = data[k];
                    }
                    data[j + 1] = v;
                    unique++;
                }
            }
        }
        return unique;
    }

    /**
     * Sort the unique indices in-place to the start of the array. The number of
     * unique indices is returned.
     *
     * <p>Uses an Order(1) data structure to ignore duplicates.
     *
     * <p>Warning: Requires {@code n > 0}. The array contents after the count of unique
     * indices is unchanged. This may change the count of each unique index in the
     * entire array.
     *
     * @param x Indices.
     * @param n Number of indices.
     * @return the number of indices
     */
    static int sortIndices(int[] x, int n) {
        // Duplicates are checked using a HashIndexSet.
        // Storage (bytes) = 4 * next-power-of-2(n*2) => 2-4 times n
        final HashIndexSet set = new HashIndexSet(n);
        int i = 0;
        int last = 0;
        set.add(x[0]);
        while (++i < n) {
            final int v = x[i];
            if (set.add(v)) {
                x[++last] = v;
            }
        }
        Arrays.sort(x, 0, ++last);
        return last;
    }
}
