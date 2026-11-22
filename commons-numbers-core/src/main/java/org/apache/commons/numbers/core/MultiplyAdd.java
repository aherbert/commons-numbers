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
package org.apache.commons.numbers.core;

import java.math.BigDecimal;

/**
 * Computes a multiply-add operation {@code a * b + c}.
 *
 * <p>This class provides an API to re-implement algorithms with a configurable
 * multiply-add operation. The class can allow switching between a standard precision
 * chained operation with 2 roundings; and a fused-multiply-add (FMA) operation with a
 * single rounding.
 *
 * <p>A software implementation of FMA is provided. This has a significant
 * performance overhead. Hardware supported FMA within Java HotSpot virtual machines (VMs)
 * can be detected at runtime and accessed via {@code Math.fma}.
 * 
 * <p>Use cases for this class include:
 * <ul>
 * <li>Assessing if an existing algorithm will increase accuracy using FMA. For
 *     studies targeting accuracy evaluation the performance penalty of the software
 *     implementation is not critical. Use {@link #FMA}. 
 * <li>Providing Java 8 code access to FMA when running within a Java 9+ HotSpot VM,
 *     otherwise use standard chained operators. Performance will not by impacted,
 *     accuracy will be improved when hardware support is detected. Use {@link #HYBRID}.
 * </ul>
 *
 * @since 1.4
 */
public enum MultiplyAdd {
    /**
     * Return the standard precision result {@code a * b + c}.
     * The result will have 2 roundings; intermediate overflow is not handled.
     */
    STD {
        @Override
        public double apply(double a, double b, double c) {
            return a * b + c;
        }
    },
    /**
     * Return the fused-multiply-add (FMA) result {@code a * b + c}.
     * The result will have {@code a * b} computed to unlimited precision without overflow
     * or underflow and the addend {@code c} added before a single rounding to the result.
     *
     * <p>This uses a hardware supported FMA instruction when available,
     * otherwise a software fallback ({@link #SOFTWARE_FMA}). The results are identical;
     * the performance changes with the software implementation.
     */
    FMA {
        @Override
        public double apply(double a, double b, double c) {
            // TODO - this should have its own method and property so
            // org.apache.commons.numbers.core.DD.fma can be ignored.
            if (FMASupport.isHardwareFMA()) {
                return FMASupport.fma(a, b, c);
            }
            return SOFTWARE_FMA.apply(a, b, c);
        }
    },
    /**
     * Return the fused-multiply-add (FMA) result {@code a * b + c}.
     * The result will have {@code a * b} computed to unlimited precision without overflow
     * or underflow and the addend {@code c} added before a single rounding to the result.
     *
     * <p>This uses a software implementation.
     */
    SOFTWARE_FMA {
        @Override
        public double apply(double a, double b, double c) {
            final int ia = getScale(a);
            final int ib = getScale(b);
            final int ic = getScale(c);
            if (ia > Double.MAX_EXPONENT || ib > Double.MAX_EXPONENT || ic > Double.MAX_EXPONENT) {
                // If any value is not finite, or is zero return the IEEE result
                return a * b + c;
            }
            // TODO: A faster software implementation:
            // - Scale a and b using the integer scaling, e.g. Math.scalb(a, -ia).
            // - Compute the the 106-bit result (x, xx) using DD two-prod.
            // - Align the mantissas for the 106-bit product and the 53-bit added;
            //   this orders terms x and c by magnitude.
            // - Use sticky summation of two smallest terms to create y;
            //   Combine the terms using add. If the high part y has the lowest
            //   bit set to 0 it can be adjusted +/- 1 ULP to carry extra magnitude from
            //   a non-zero fast-two-sum low part.
            // - Add y to largest term; the sticky bit correctly rounds the result.
            // - Rescale the result
            return new BigDecimal(a).multiply(new BigDecimal(b)).add(new BigDecimal(c)).doubleValue();
        }
    },
    /**
     * Return the multiply-add result {@code a * b + c}.
     *
     * <p>This uses a hardware supported FMA instruction when available ({@link #FMA}),
     * otherwise a software fallback ({@link #SOFTWARE_FMA}). The results may change when
     * the FMA instruction is used; performance is identical.
     */
    HYBRID {
        @Override
        public double apply(double a, double b, double c) {
            // TODO - this should have its own method and property so
            // org.apache.commons.numbers.core.DD.fma can be ignored.
            if (FMASupport.isHardwareFMA()) {
                return FMASupport.fma(a, b, c);
            }
            return STD.apply(a, b, c);
        }
    };

    /** Exponent offset in IEEE754 representation. */
    private static final int EXPONENT_OFFSET = 1023;
    /** Mask to extract the 52-bit mantissa from a long representation of a double. */
    private static final long MANTISSA_MASK = 0x000f_ffff_ffff_ffffL;

    /**
     * Compute the multipy-add {@code a * b + c}.
     *
     * @param a Factor
     * @param b Factor
     * @param c Addend
     * @return {@code a * b + c}
     */
    public abstract double apply(double a, double b, double c);

    /**
     * Returns a scale suitable for use with {@link Math#scalb(double, int)} to normalise
     * the number to the interval {@code [1, 2)}.
     *
     * <p>The scale is typically the largest unbiased exponent used in the representation of the
     * two numbers. In contrast to {@link Math#getExponent(double)} this handles
     * sub-normal numbers by computing the number of leading zeros in the mantissa
     * and shifting the unbiased exponent. The result is that for all finite, non-zero,
     * numbers {@code x}, the magnitude of {@code scalb(x, -getScale(x))} is
     * always in the range {@code [1, 2)}.
     *
     * <p>This method is a functional equivalent of the c function ilogb(double) adapted for
     * two input arguments.
     *
     * <p>The result is to be used to scale a number using {@link Math#scalb(double, int)}.
     * Hence the special case of a zero argument is handled using the return value for NaN
     * as zero cannot be scaled. This is different from {@link Math#getExponent(double)}.
     *
     * <p>Special cases:
     *
     * <ul>
     * <li>If either argument is NaN or infinite, then the result is
     * {@link Double#MAX_EXPONENT} + 1.</li>
     * <li>If both arguments are zero, then the result is
     * {@link Double#MAX_EXPONENT} + 1.</li>
     * </ul>
     *
     * @param x the value
     * @return The maximum unbiased exponent of the value to be used for scaling
     * @see Math#getExponent(double)
     * @see Math#scalb(double, int)
     * @see <a href="https://www.cplusplus.com/reference/cmath/ilogb/">ilogb</a>
     */
    private static int getScale(double x) {
        // Only interested in the exponent and mantissa so remove the sign bit
        final long bits = Double.doubleToRawLongBits(x) & Long.MAX_VALUE;
        // Get the unbiased exponent
        int exp = (int) (bits >>> 52) - EXPONENT_OFFSET;

        // No case to distinguish nan/inf
        // Handle sub-normal numbers
        if (exp == Double.MIN_EXPONENT - 1) {
            // Special case for zero, return as nan/inf to indicate scaling is not possible
            if (bits == 0) {
                return Double.MAX_EXPONENT + 1;
            }
            // A sub-normal number has an exponent below -1022. The amount below
            // is defined by the number of shifts of the most significant bit in
            // the mantissa that is required to get a 1 at position 53 (i.e. as
            // if it were a normal number with assumed leading bit)
            final long mantissa = bits & MANTISSA_MASK;
            exp -= Long.numberOfLeadingZeros(mantissa << 12);
        }
        return exp;
    }
}
