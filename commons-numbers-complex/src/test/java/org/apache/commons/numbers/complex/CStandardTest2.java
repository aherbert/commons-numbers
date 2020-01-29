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

package org.apache.commons.numbers.complex;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.function.ToDoubleFunction;

/**
 * Tests the standards defined by the C.99 standard for complex numbers
 * defined in ISO/IEC 9899, Annex G.
 *
 * @see <a href="http://www.open-std.org/JTC1/SC22/WG14/www/standards">
 *    ISO/IEC 9899 - Programming languages - C</a>
 */
public class CStandardTest2 {
    private static final double inf = Double.POSITIVE_INFINITY;
    private static final double negInf = Double.NEGATIVE_INFINITY;
    private static final double nan = Double.NaN;
    private static final BigDecimal TWO = BigDecimal.valueOf(2);

    /** 60 shifted 20-bits to align with the exponent of the upper 32-bits of a double. */
    private static final int EXP_60 = 0x3c_00000;
    /** Represents an exponent of 500 in unbiased form shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_500 = 0x5f3_00000;
    /** Represents an exponent of 1024 in unbiased form (infinite or nan)
     * shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_1024 = 0x7ff_00000;
    /** Represents an exponent of -500 in unbiased form shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_NEG_500 = 0x20b_00000;
    /** Represents an exponent of -1022 in unbiased form (sub-normal number)
     * shifted 20-bits to align with the upper 32-bits of a double. */
    private static final int EXP_NEG_1022 = 0x001_00000;

    private static final char[] TMP = new char[64];
    private static final char[] DIGITS = {'0', '1'};

    private static long hp;
    private static long hp1;
    private static long lp;
    private static long better;
    private static long worse;
    private static long total;
    private static long same;
    private static long correct;

    /**
     * Assert {@link Complex#abs()} is functionally equivalent to using
     * {@link Math#hypot(double, double)}. If the results differ the true result
     * is computed using BigDecimal. The test fails if the result is further than
     * the provided ULPs from the true result.
     *
     * <p>This can be used to assert that the custom implementation of abs() is no worse than
     * {@link Math#hypot(double, double)} which aims to be within 1 ULP of the exact result.
     *
     * <p>Note: This method will not handle an input complex that is infinite or nan.
     *
     * @param z the complex
     * @param ulps the maximum allowed ULPs from the exact result
     */
    private static void assertAbs(Complex z, int ulps) {
        double x = z.getReal();
        double y = z.getImaginary();
        final double e = Math.hypot(x, y);
        final double o = z.abs();
        //if (e == o) {
        //    return;
        //}

        // TODO: For benchmarking use LinearCombination not BigDecimal ...
        // Add a counter for deviations from the correct answer.
        // Add scaling onl when necessary

        // Do scaling for sub-normals. Assume x > y
        final double exact = hypot(x, y);
//        final int scale = Math.getExponent(y);
//        final double xs = Math.scalb(x, -scale);
//        final double ys = Math.scalb(y, -scale);
//        final double exact = Math.scalb(
//            Math.sqrt(new BigDecimal(xs).pow(2).add(new BigDecimal(ys).pow(2)).doubleValue()), scale);

        final double standard = hypot2(x, y);
        total++;

        if (exact == o) {
            correct++;
        }
        if (exact == e && exact == standard) {
            return;
        }

        final long ebits = Double.doubleToLongBits(e);
        final long obits = Double.doubleToLongBits(o);
        final long bits = Double.doubleToLongBits(exact);
        final long bits2 = Double.doubleToLongBits(standard);
        final long ulpe = Math.abs(bits - ebits);
        final long ulpo = Math.abs(bits - obits);
        final long ulpst = Math.abs(bits - bits2);

        //if (ulpo < Math.abs(bits - bits2) ) {
        if (ulpst != 0) {
            // standard computation is not enough
            hp++;
        }

        // CHECKSTYLE: stop Regexp

        // Is the computation better than standard precision
        if (ulpo < ulpst) {
            hp1++;
        } else if (ulpo > ulpst) {
            lp++;
//          System.out.printf("%s.abs(). %.2f  Exact %s : abs() %s (%d ulps). hypot() %s (%d ulps)%n%s%n%s%n%s%n%s%n",
//          z, Math.max(x, y) / Math.min(x, y),
//          exact, o, ulpo, e, ulpe,
//          format(x), format(y),
//          format(x * 0x1.0p1022), format(y * 0x1.0p1022)
//  );
        }

        if (ulpo < ulpe) {
            better++;
        } else if (ulpo > ulpe) {
//            System.out.printf("%s.abs(). %.2f  Exact %s : abs() %s (%d ulps). hypot() %s (%d ulps)%n%s%n%s%n%s%n%s%n",
//                    z, Math.max(x, y) / Math.min(x, y),
//                    exact, o, ulpo, e, ulpe,
//                    format(x), format(y),
//                    format(x * 0x1.0p1022), format(y * 0x1.0p1022)
//            );
            worse++;
        } else {
            same++;
        }

        if (exact == o) {
            return;
        }
        // Distance from the exact result should be comparable
//        final long ebits = Double.doubleToLongBits(e);
//        final long obits = Double.doubleToLongBits(o);
//        final long bits = Double.doubleToLongBits(exact);
//        final long ulpe = Math.abs(bits - ebits);
//        final long ulpo = Math.abs(bits - obits);
        // Used to debug differences in the functions.
//        System.out.printf("%s.abs(). Exact %s : abs() %s (%d ulps). hypot() %s (%d ulps)%n",
//            z, exact, o, ulpo, e, ulpe);
        Assertions.assertTrue(ulpo <= ulps, () ->
            String.format("%s.abs(). Expected %s, was %s (%d ulps). hypot %s (%d ulps)",
                z, exact, o, ulpo, e, ulpe));
    }
    // CHECKSTYLE: resume Regexp

    private static String format(double x) {
        // Full length representation with no sign bit
        long bits = Double.doubleToRawLongBits(x);
        // 52 bit mantissa
        for (int i = 63; i >= 11; i--) {
            TMP[i] = DIGITS[(int)(bits & 0x1)];
            bits >>>= 1;
        }
        TMP[11] = ' ';
        // 11 bit exponent
        for (int i = 10; i >= 0; i--) {
            TMP[i] = DIGITS[(int)(bits & 0x1)];
            bits >>>= 1;
        }
        return Math.getExponent(x) + " " + String.valueOf(TMP);
    }

    // https://stackoverflow.com/questions/13649703/square-root-of-bigdecimal-in-java
    private static BigDecimal sqrt(BigDecimal a, final int scale) {
        BigDecimal x0 = BigDecimal.ZERO;
        BigDecimal x1 = new BigDecimal(Math.sqrt(a.doubleValue()));
        while (!x0.equals(x1)) {
            x0 = x1;
            x1 = a.divide(x0, MathContext.DECIMAL128);
            x1 = x1.add(x0);
            x1 = x1.divide(TWO, MathContext.DECIMAL128);

        }
        return x1;
    }

    private static double hypot(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This can incorrectly order small sub-normal numbers which differ
        // only in the lower 32-bits for the x^2+y^2 sum which requires x to be
        // greater than y. This version performs a second reorder
        // if the upper parts are equal. This can be done after the test for large magnitude
        // differences. Note: An alternative using the entire 63-bit unsigned raw bits for the
        // magnitude comparison are marginally slower in performance tests. For alignment
        // with the fdlibm reference this implementation uses the upper 32-bits.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by direct writing to the exponent bits.
        // The high part was maintained through scaling for using in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation is not tied to
        // manipulating bits to represent the magnitude of the high part of the split number.
        // The split is computed dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases. In the event of worse performance
        // the result is 1 ULP from the exact result.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fffffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fffffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_60) {
            /* x/y > 2**60 */
            // The addition will return nan for finite/nan combinations.
            // Note: if both inputs are inf or nan this will not occur as the exponent is the
            // same and inf+nan are handled later.
            return a + b;
        }

        // Second re-order in the rare event the upper 32-bits are the same.
        // This reorder will ignore inf-nan combinations. These are rejected below.
        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // No bit manipulation on high and low 32-bits to check for a zero mantissa (Inf).
                // Revert to Java functions for the IEEE754 result.
                return a == Double.POSITIVE_INFINITY || b == Double.POSITIVE_INFINITY ?
                    Double.POSITIVE_INFINITY :
                    /* for sNaN */
                    a + b;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p1022;
                b *= 0x1.0p1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p600;
                b *= 0x1.0p600;
                rescale = 0x1.0p-600;
            }
        }
        //double x2y2 = new BigDecimal(a).pow(2).add(new BigDecimal(b).pow(2)).doubleValue();
        double x2y2 = value(a, b);
        return Math.sqrt(x2y2) * rescale;
    }

    private static double hypot2(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This can incorrectly order small sub-normal numbers which differ
        // only in the lower 32-bits for the x^2+y^2 sum which requires x to be
        // greater than y. This version performs a second reorder
        // if the upper parts are equal. This can be done after the test for large magnitude
        // differences. Note: An alternative using the entire 63-bit unsigned raw bits for the
        // magnitude comparison are marginally slower in performance tests. For alignment
        // with the fdlibm reference this implementation uses the upper 32-bits.
        //
        // 2. This stores the re-scaling factor for use in a multiplication.
        // The original computed scaling by direct writing to the exponent bits.
        // The high part was maintained through scaling for using in the high
        // precision sum x^2 + y^2. This version has no requirement for that
        // as the high part is extracted using Dekker's method.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Since the exponent difference between a and b
        // is below 60, if a's exponent is above 500 then b's cannot be below
        // -500 even after scaling by -600 in the first conditional:
        // ((>500 - 60) - 600) > -500
        //
        // 4. The original handling of sub-normal numbers does not scale the representation
        // of the upper 32-bits (ha, hb). This is because sub-normals can have effective
        // exponents in the range -1023 to -1074 and the exponent increase must be dynamically
        // computed for an addition to the exponent bits. This is neglected in the fdlibm
        // version. The effect is the extended precision computation reverts to a
        // computation as if t1 and y1 were zero (see x2y2() below) and the result is the
        // standard precision result x^2 + y^2. In contrast this implementation is not tied to
        // manipulating bits to represent the magnitude of the high part of the split number.
        // The split is computed dynamically on the scaled number. The effect is increased
        // precision for the majority of sub-normal cases. In the event of worse performance
        // the result is 1 ULP from the exact result.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // For clarity the values have been refactored to constants.
        //
        // Scaling factors are written using Java's binary floating point notation 0xN.NpE
        // where the number following the p is the exponent. These are only used with
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & 0x7fffffff;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & 0x7fffffff;

        // Order by approximate magnitude (lower bits excluded)
        double a;
        double b;
        if (hb > ha) {
            a = y;
            b = x;
            final int j = ha;
            ha = hb;
            hb = j;
        } else {
            a = x;
            b = y;
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_60) {
            /* x/y > 2**60 */
            // The addition will return nan for finite/nan combinations.
            // Note: if both inputs are inf or nan this will not occur as the exponent is the
            // same and inf+nan are handled later.
            return a + b;
        }

        // Second re-order in the rare event the upper 32-bits are the same.
        // This reorder will ignore inf-nan combinations. These are rejected below.
        if (ha == hb && a < b) {
            final double tmp = a;
            a = b;
            b = tmp;
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // No bit manipulation on high and low 32-bits to check for a zero mantissa (Inf).
                // Revert to Java functions for the IEEE754 result.
                return a == Double.POSITIVE_INFINITY || b == Double.POSITIVE_INFINITY ?
                    Double.POSITIVE_INFINITY :
                    /* for sNaN */
                    a + b;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p1022;
                b *= 0x1.0p1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p600;
                b *= 0x1.0p600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(a * a + b * b) * rescale;
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @return \( a1^2 + b1^2 \)
     */
    public static double value(double a1, double b1) {
        // Initial product creates the expansion e.
        // This is merged with each product expansion f using a linear expansion sum.
        double e0;
        double e1;
        double e2;
        double e3;
        double f0;
        double f1;

        // q is the running scalar product, x & a are working variables
        double x;
        double a;
        e1 = a1 * a1;
        e0 = productLow(a1, e1);

        // First sum of the expansion f creates e[0-3]
        f1 = b1 * b1;
        f0 = productLow(b1, f1);
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

        // Final summation
        return e0 + e1 + e2 + e3;
    }

    private static double productLow(double x, double xy) {
        // Split the numbers using Dekker's algorithm
        final double hx = highPart(x);
        final double lx = x - hx;

        // Compute the multiply low part:
        // err1 = xy - hx * hy
        // err2 = err1 - lx * hy
        // err3 = err2 - hx * ly
        // low = lx * ly - err3
        return lx * lx - ((xy - hx * hx) - 2 * lx * hx);
    }

    private static double sumLow(double a, double b, double sum) {
        final double bVirtual = sum - a;
        // sum - bVirtual == aVirtual.
        // a - aVirtual == a round-off
        // b - bVirtual == b round-off
        return (a - (sum - bVirtual)) + (b - bVirtual);
    }

    private static double highPart(double value) {
        // normal conversion
        final double c = 1.34217729E8 * value;
        return c - (c - value);
    }

    /**
     * Assert {@link Complex#abs()} functions as per {@link Math#hypot(double, double)}.
     * The two numbers for {@code z = x + iy} are generated from the two function.
     *
     * <p>The functions should not generate numbers that are infinite or nan.
     *
     * @param rng Source of randomness
     * @param fx Function to generate x
     * @param fy Function to generate y
     * @param samples Number of samples
     */
    private static void assertAbs(UniformRandomProvider rng,
                                  ToDoubleFunction<UniformRandomProvider> fx,
                                  ToDoubleFunction<UniformRandomProvider> fy,
                                  int samples) {
        for (int i = 0; i < samples; i++) {
            double x = fx.applyAsDouble(rng);
            double y = fy.applyAsDouble(rng);
            assertAbs(Complex.ofCartesian(x, y), 1);
        }
    }

    /**
     * Creates a sub-normal number with up to 52-bits in the mantissa.
     *
     * @param rng Source of randomness
     * @return the number
     */
    private static double createSubNormalNumber52(UniformRandomProvider rng) {
        return Double.longBitsToDouble(rng.nextLong() >>> 12);
    }

    /**
     * Creates a sub-normal number with up to 32-bits in the mantissa.
     *
     * @param rng Source of randomness
     * @return the number
     */
    private static double createSubNormalNumber32(UniformRandomProvider rng) {
        return Double.longBitsToDouble(rng.nextLong() >>> 32);
    }

    /**
     * Creates a number in the range {@code [1, 2)} with up to 52-bits in the mantissa.
     *
     * @param rng Source of randomness
     * @return the number
     */
    private static double createFixedExponentNumber(UniformRandomProvider rng) {
        return Double.longBitsToDouble((rng.nextLong() >>> 12) | (1023L << 52));
    }

    /**
     * Utility to create a Complex.
     *
     * @param real the real
     * @param imaginary the imaginary
     * @return the complex
     */
    private static Complex complex(double real, double imaginary) {
        return Complex.ofCartesian(real, imaginary);
    }

    /**
     * Create a number in the range {@code [-5,5)}.
     *
     * @param rng the random generator
     * @return the number
     */
    private static double next(UniformRandomProvider rng) {
        // Note: [0, 1) minus 1 is [1, 0). This occurs half the time to create [-1, 1).
        return (rng.nextDouble() - rng.nextInt(1)) * 5;
    }

    /**
     * ISO C Standard G.6 (6) for abs().
     * Defined by ISO C Standard F.9.4.3 hypot function.
     */
    @Test
    public void testAbs() {
        //assertAbs(complex(1.58917270352193E-308, 1.885544663211271E-308), 1);
        //assertAbs(complex(1.58917270352193E-308 * 8, 1.885544663211271E-308 * 8), 1);

        Assertions.assertEquals(inf, complex(inf, nan).abs());
        Assertions.assertEquals(inf, complex(negInf, nan).abs());
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP, 5267354726736L);
        for (int i = 0; i < 10; i++) {
            final double x = next(rng);
            final double y = next(rng);
            Assertions.assertEquals(complex(x, y).abs(), complex(y, x).abs());
            Assertions.assertEquals(complex(x, y).abs(), complex(x, -y).abs());
            Assertions.assertEquals(Math.abs(x), complex(x, 0.0).abs());
            Assertions.assertEquals(Math.abs(x), complex(x, -0.0).abs());
            Assertions.assertEquals(inf, complex(inf, y).abs());
            Assertions.assertEquals(inf, complex(negInf, y).abs());
        }

        // Test verses Math.hypot due to the use of a custom implementation.
        // First test edge cases. Ignore negatives as the sign is simply removed.
        final double[] parts = {0, Double.MIN_VALUE, Double.MIN_NORMAL, Double.MAX_VALUE,
            Double.POSITIVE_INFINITY, Double.NaN};
        for (final double x : parts) {
            for (final double y : parts) {
                Assertions.assertEquals(Math.hypot(x, y), complex(x, y).abs());
            }
        }

        // Test with a range of numbers.
        // Sub-normals require special handling so we use different variations of these.
        int samples = 1 << 30; //Integer.MAX_VALUE;
        //assertAbs(rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber52, samples);
        //report("sub-normal 52");

        // Numbers on different scales
        assertAbs(rng, CStandardTest2::createFixedExponentNumber, CStandardTest2::createFixedExponentNumber, samples);
        report("range0");
        assertAbs(rng, CStandardTest2::createFixedExponentNumber,
            r -> createFixedExponentNumber(r) * 0x1.0p+1, samples);
        report("range1");
        assertAbs(rng, CStandardTest2::createFixedExponentNumber,
            r -> createFixedExponentNumber(r) * 0x1.0p+2, samples);
        report("range2");
        assertAbs(rng, CStandardTest2::createFixedExponentNumber,
            r -> createFixedExponentNumber(r) * 0x1.0p+3, samples);
        report("range3");
        assertAbs(rng, CStandardTest2::createFixedExponentNumber,
            r -> createFixedExponentNumber(r) * 0x1.0p+4, samples);
        report("range4");
        assertAbs(rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber52, samples);
        report("sub-normal 52");
        assertAbs(rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber32, samples);
        report("sub-normal 52/32");
        assertAbs(rng, CStandardTest2::createSubNormalNumber32, CStandardTest2::createSubNormalNumber32, samples);
        report("sub-normal 32");
        assertAbs(rng, UniformRandomProvider::nextDouble, UniformRandomProvider::nextDouble, samples);
        report("uniform");
        // The distribution of floating-point numbers has the log-uniform distribution as
        // its limiting distribution. We do not use that to force a full 52-bit random mantissa
        // which increase the likelihood for high precision.
        final ToDoubleFunction<UniformRandomProvider> luGenerator = new ToDoubleFunction<UniformRandomProvider>() {
            @Override
            public double applyAsDouble(UniformRandomProvider rng) {
                return Math.scalb(createFixedExponentNumber(rng), rng.nextInt(30));
            }
        };
        assertAbs(rng, luGenerator, luGenerator, samples);
        report("log-uniform");
        // Complex cis numbers
        final ToDoubleFunction<UniformRandomProvider> cisGenerator = new ToDoubleFunction<UniformRandomProvider>() {
            private double tmp = Double.NaN;
            @Override
            public double applyAsDouble(UniformRandomProvider rng) {
                if (Double.isNaN(tmp)) {
                    double u = rng.nextDouble() * Math.PI;
                    tmp = Math.cos(u);
                    return Math.sin(u);
                }
                final double r = tmp;
                tmp = Double.NaN;
                return r;
            }
        };
        assertAbs(rng, cisGenerator, cisGenerator, samples);
        report("cis");
    }

    private static void report(String type) {
        // CHECKSTYLE: stop all
        System.out.printf("%-17s  %d / %d (%.3f) : high-precision %d (%.3f) : hp %d (%.3f) : lp %d (%.3f) : better %d (%.3f) : worse %d (%.3f) : same %d (%.3f)%n",
                type,
                correct, total, (double) correct / total,
                hp, (double) hp / total,
                hp1, (double) hp1 / total,
                lp, (double) lp / total,
                better, (double) better / total,
                worse, (double) worse / total,
                same, (double) same / total
                );
        // CHECKSTYLE: resume all
        better = worse = total = same = correct = hp = hp1 = lp = 0;
    }

    @Test
    public void testLower32() {
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
        final long exp = Double.doubleToLongBits(1.0);
        final long total2 = 1L << 34;
        long diff = 0;
        for (long l = total2; l-- > 0;) {
            // 52-bits
            final long bits = rng.nextLong() >>> 12;
            // A different lower 32-bits
            final int lower = rng.nextInt();
            final long a = exp | bits;
            final long b = exp | (bits ^ (lower & 0xffff_ffffL));
            final double x = Double.longBitsToDouble(a);
            final double y = Double.longBitsToDouble(b);
            //final double e = x2y2(x, y);
            //final double o = x2y2(y, x);
            //final double e = fma(x, y);
            //final double o = fma(y, x);
            final double e = x2y2Dekker(x, y);
            final double o = x2y2Dekker(y, x);
            if (e != o) {
                diff++;
                // CHECKSTYLE: stop all
                System.out.printf("{%s, %s}%n", x, y);
            }
        }
        System.out.printf("Lower diff %d / %d (%s)%n", diff, total2, (double) diff / total2);

        // If this is OK then do a full 2^52 combinations test: 2^20 uppers * 2^32 lowers.
        // Do this multi threaded.
    }

    @Test
    public void testHypotUnderScaling() {
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
        final long total2 = 1L << 34;
        long diff = 0;
        for (long l = total2; l-- > 0;) {
            // Create a known sub-normal with 52-bits
            final double x = Double.longBitsToDouble(rng.nextLong() >>> 52);
            final double y = Double.longBitsToDouble(rng.nextLong() >>> 52);
            final double e = Math.hypot(x, y) * 0x1.0p+300;
            final double o = Math.hypot(x * 0x1.0p+300, y * 0x1.0p+300);
            if (e != o) {
                // CHECKSTYLE: stop all
                System.out.printf("{%s, %s} = %s %s%n", x, y, e, o);
                final double e1 = complex(x, y).abs() * 0x1.0p+300;
                final double o1 = complex(x * 0x1.0p+300, y * 0x1.0p+300).abs();
                if (e != o && e1 == o1) {
                    System.out.printf("{%s, %s}%n", x, y);
                    if (diff++ == 10) {
                        break;
                    }
                }
            }
        }
    }

    /**
     * High accuracy {@code x^2 + y^2}.
     *
     * @param x the x
     * @param y the y
     * @return {@code x^2 + y^2}
     */
    private static double x2y2(double x, double y) {
        // Copied from Complex for the case where x and y are close
        // 2y > x > y
        final double w = x - y;
        final double t = x + x;
        final double y1 = splitHigh(y);
        final double y2 = y - y1;
        final double t1 = splitHigh(t);
        final double t2 = t - t1;
        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));
    }

    /**
     * High accuracy {@code x^2 + y^2}.
     *
     * @param x the x
     * @param y the y
     * @return {@code x^2 + y^2}
     */
    private static double x2y2Dekker(double x, double y) {
        final double xx = x * x;
        final double yy = y * y;
        final double xHigh = splitHigh(x);
        final double xLow = x - xHigh;
        final double yHigh = splitHigh(y);
        final double yLow = y - yHigh;
        final double x2Low = squareLow(xLow, xHigh, xx);
        final double y2Low = squareLow(yLow, yHigh, yy);

        final double r = xx + yy;
        return xx - r + yy + y2Low + x2Low + r;
    }

    /**
     * High accuracy {@code x^2 + y^2}.
     *
     * @param x the x
     * @param y the y
     * @return {@code x^2 + y^2}
     */
    private static double fma(double x, double y) {
        // Copied from Complex for the case where x and y are close
        // 2y > x > y
        final double yy = y * y;
        final double xx = x * x;
        final double xHigh = splitHigh(x);
        final double xLow = x - xHigh;
        final double x2Low = squareLow(xLow, xHigh, xx);
        // Two sum y^2 into the expansion of x^2
        final double r = xx + yy;
        return xx - r + yy + x2Low + r;
    }

    /**
     * Compute the round-off from the square of a split number with {@code low} and {@code high}
     * components. Uses Dekker's algorithm for split multiplication modified for a square product.
     *
     * @param low Low part of number.
     * @param high High part of number.
     * @param square Square of the number.
     * @return <code>low * low - ((product - high * high) - 2 * low * high)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    private static double squareLow(double low, double high, double square) {
        return low * low - ((square - high * high) - 2 * low * high);
    }

    /**
     * Compute the upper part of the number using the Dekker split.
     *
     * @param a the number
     * @return the upper part
     */
    private static double splitHigh(double a) {
        final double c = 1.34217729E8 * a;
        return c - (c - a);
    }

    @Test
    public void testUniform() {
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
        final long total2 = 1L << 20;
        long diff = 0;
        for (long l = total2; l-- > 0;) {
            double x = rng.nextDouble();
            double y = rng.nextDouble();
            x = Math.max(x, y);
            y = Math.min(x, y);
            if (x < 2 * y) {
                diff++;
                // CHECKSTYLE: stop all
                //System.out.printf("{%s, %s}%n", x, y);
            }
        }
        System.out.printf("Lower diff %d / %d (%s)%n", diff, total2, (double) diff / total2);
    }
}
