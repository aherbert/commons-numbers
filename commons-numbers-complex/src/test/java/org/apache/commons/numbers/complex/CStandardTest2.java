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

import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Formatter;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.DoubleSupplier;
import java.util.function.Function;
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

    private interface DoubleDoubleBiFunction {
        double apply(double a, double b);
    }

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
        final double exact = hypotExact(x, y);
//        final int scale = Math.getExponent(y);
//        final double xs = Math.scalb(x, -scale);
//        final double ys = Math.scalb(y, -scale);
//        final double exact = Math.scalb(
//            Math.sqrt(new BigDecimal(xs).pow(2).add(new BigDecimal(ys).pow(2)).doubleValue()), scale);

        final double standard = hypotStandard(x, y);
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

    private static class UlpChecker implements Callable<UlpChecker> {
        final DoubleSupplier fx;
        final DoubleSupplier fy;
        final int samples;

        long same;
        long d1, d2, d3, d4;
        long s1, s2, s3, s4;
        // Special counter to compare hypot implementations.
        long better, worse;
        BigDecimal m1 = BigDecimal.ZERO, m2 = BigDecimal.ZERO, m3 = BigDecimal.ZERO, m4 = BigDecimal.ZERO;
        Complex z1 = Complex.ZERO, z2 = Complex.ZERO, z3 = Complex.ZERO, z4 = Complex.ZERO;

        UlpChecker(DoubleSupplier fx,
            DoubleSupplier fy, int samples) {
            this.fx = fx;
            this.fy = fy;
            this.samples = samples;
        }

        UlpChecker add(UlpChecker c) {
            if (c != null) {
                same += c.same;
                d1 += c.d1;
                d2 += c.d2;
                d3 += c.d3;
                d4 += c.d4;
                s1 += c.s1;
                s2 += c.s2;
                s3 += c.s3;
                s4 += c.s4;
                better += c.better;
                worse += c.worse;
                if (m1.compareTo(c.m1) < 0) {
                    m1 = c.m1;
                    z1 = c.z1;
                }
                if (m2.compareTo(c.m2) < 0) {
                    m2 = c.m2;
                    z2 = c.z2;
                }
                if (m3.compareTo(c.m3) < 0) {
                    m3 = c.m3;
                    z3 = c.z3;
                }
                if (m4.compareTo(c.m4) < 0) {
                    m4 = c.m4;
                    z4 = c.z4;
                }
            }
            return this;
        }

        @Override
        public UlpChecker call() {
            for (int i = 0; i < samples; i++) {
                check();
            }
            return this;
        }

        /**
         * Check the difference between the computation of {@code x^2 + y^2} using
         * {@link Math#fma(double, double, double)} and
         * {@link Math#hypot(double, double)}. If the results differ the true result
         * is computed using high precision.BigDecimal and the ULPs computed from the
         * extended precision result.
         *
         * <p>This can be used to assert that the custom implementation of abs() is no worse than
         * {@link Math#hypot(double, double)} which aims to be within 1 ULP of the exact result.
         *
         * <p>Note: This method will not handle an input complex that is infinite or nan.
         */
        void check() {
            double x = Math.abs(fx.getAsDouble());
            double y = Math.abs(fy.getAsDouble());

            if (x < y) {
                double tmp = x;
                x = y;
                y = tmp;
            }

            // High precision dekker method as reference
            double e = x2y2DekkerSqrt(x, y);

            // high precision sum then standard sqrt
            double o1 = Math.sqrt(x2y2Dekker(x, y));
            //double o1 = Math.sqrt(Math.fma(x, x, y * y));
            // Alternative using FMA
            double o2 = Math.sqrt(Math.fma(x, x, y * y));
            //double o2 = hypot3(x, y);
            // Reference
            double o3 = Math.hypot(x, y);
            // Use the variant of hypot to contrast with FMA.
            double o4 = Math.sqrt(x2y2Hypot(x, y));
            //double o4 = Complex.ofCartesian(x, y).abs();

            // Only interested in differences from realistic best result without BigDecimal or
            // extended precision sqrt (although we could add this using Dekker's method)
            if (e == o1 && e == o2 && e == o3 && e == o4) {
                same++;
                return;
            }

            // Verify the exact result
            BigDecimal bdexact = x2y2BigDecimal(x, y).sqrt(MathContext.DECIMAL128);

            score(x, y, o1, o2, o3, o4, bdexact);
        }
        // CHECKSTYLE: resume Regexp

        void score(double x, double y, double o1, double o2, double o3, double o4, BigDecimal bdexact) {
            final double exact = bdexact.doubleValue();
            if (exact != o1) {
                d1++;
            }
            if (exact != o2) {
                d2++;
            }
            // Assume
            // o3 = Math.hypot
            // o4 == hypot reimplemented
            if (exact != o3) {
                d3++;
                if (o4 == exact) {
                    better++;
                }
            }
            if (exact != o4) {
                d4++;
                if (o3 == exact) {
                    worse++;
                }
            }

            // Compute ULPs relative to the answer as a double.
            // Kahan (eecs.berkeley.edu/~wkahan/LOG10HAF.TXT) defines ulp(x) as
            // "the gap between the two floating-point numbers nearest x".
            // This specifies ulp(1) as epsilon/2 not epsilon, i.e. the ulp distance
            // of a power of 2 is the dsitance to the next lower number:
            // ---|------|------------|------------|
            //  1-eps/2  1          1+eps       1+2eps
            // If we use this definition a computation of the next float above 1 would be
            // 2 ulps. Do this and see the effect on the results.

            // Determine if the BigDecimal was rounded up and use the next down.
            // This will have two-fold effect on ulps when testing cis numbers where the radius
            // is computed above one.
            BigDecimal exactUlp;
            if (new BigDecimal(exact).compareTo(bdexact) >= 0) {
                exactUlp = new BigDecimal(Math.ulp(Math.nextDown(exact)));
            } else {
                exactUlp = new BigDecimal(Math.ulp(exact));
            }

            final BigDecimal ulpo1 = new BigDecimal(o1).subtract(bdexact).abs().divide(exactUlp, 5, RoundingMode.HALF_UP);
            BigDecimal ulpo2;
            BigDecimal ulpo3;
            BigDecimal ulpo4;
            if (o2 == o1) {
                ulpo2 = ulpo1;
            } else {
                ulpo2 = new BigDecimal(o2).subtract(bdexact).abs().divide(exactUlp, 5, RoundingMode.HALF_UP);
            }
            if (o3 == o1) {
                ulpo3 = ulpo1;
            } else if (o3 == o2) {
                ulpo3 = ulpo2;
            } else {
                ulpo3 = new BigDecimal(o3).subtract(bdexact).abs().divide(exactUlp, 5, RoundingMode.HALF_UP);
            }
            if (o4 == o1) {
                ulpo4 = ulpo1;
            } else if (o4 == o2) {
                ulpo4 = ulpo2;
            } else if (o4 == o3) {
                ulpo4 = ulpo3;
            } else {
                ulpo4 = new BigDecimal(o4).subtract(bdexact).abs().divide(exactUlp, 5, RoundingMode.HALF_UP);
            }
            // CHECKSTYLE: stop Regexp
//            System.out.printf("%s,%s = %s : dekker %s (%s) %s %s %s : fma %s (%s) : hypot %s (%s) : hypot2 %s (%s)%n", 
//                    x, y, bdexact, o1, ulpo1, 
//                    new BigDecimal(o1).subtract(bdexact).abs(),
//                    new BigDecimal(o1).subtract(bdexact).abs().divide(exactUlp, 5, RoundingMode.HALF_UP),
//                    new BigDecimal(o1).subtract(bdexact).abs().divide(exactUlp, 10, RoundingMode.HALF_UP),
//                    o2, ulpo2, o3, ulpo3, o4, ulpo4);

            // Accrue differences
            s1 += ulpo1.scaleByPowerOfTen(5).longValue();
            if (m1.compareTo(ulpo1) < 0) {
                m1 = ulpo1;
                z1 = Complex.ofCartesian(x, y);
            }
            s2 += ulpo2.scaleByPowerOfTen(5).longValue();
            if (m2.compareTo(ulpo2) < 0) {
                m2 = ulpo2;
                z2 = Complex.ofCartesian(x, y);
            }
            s3 += ulpo3.scaleByPowerOfTen(5).longValue();
            if (m3.compareTo(ulpo3) < 0) {
                m3 = ulpo3;
                z3 = Complex.ofCartesian(x, y);
            }
            s4 += ulpo4.scaleByPowerOfTen(5).longValue();
            if (m4.compareTo(ulpo4) < 0) {
                m4 = ulpo4;
                z4 = Complex.ofCartesian(x, y);
            }
        }
    }

    private static class ScaledUlpChecker extends UlpChecker {
        ScaledUlpChecker(DoubleSupplier fx, DoubleSupplier fy, int samples) {
            super(fx, fy, samples);
        }

        /**
         * {@inheritDoc}
         *
         * <p>Note: This method will handle an input complex that is extreme valued, e.g. sub-normal.
         */
        @Override
        void check() {
            double x = Math.abs(fx.getAsDouble());
            double y = Math.abs(fy.getAsDouble());

            if (x < y) {
                double tmp = x;
                x = y;
                y = tmp;
            }

            // Do scaling once
            double rescale = 1.0;
            if (Math.min(x, y) < 0x1.0p-500) {
                if (Math.min(x, y) < Double.MIN_NORMAL) {
                    x *= 0x1.0p+1022;
                    y *= 0x1.0p+1022;
                    rescale = 0x1.0p-1022;
                } else {
                    x *= 0x1.0p+600;
                    y *= 0x1.0p+600;
                    rescale = 0x1.0p-600;
                }
            } else if (Math.max(x, y) > 0x1.0p500) {
                x *= 0x1.0p-600;
                y *= 0x1.0p-600;
                rescale = 0x1.0p+600;
            }

            // High precision dekker method as reference
            double e = x2y2DekkerSqrt(x, y) * rescale;

            // high precision sum then standard sqrt
            double o1 = Math.sqrt(x2y2Dekker(x, y)) * rescale;
            //double o1 = Math.sqrt(Math.fma(x, x, y * y)) * rescale;
            // Alternative using FMA
            double o2 = Math.sqrt(Math.fma(x, x, y * y)) * rescale;
            //double o2 = hypot3(x, y) * rescale;
            // Reference
            double o3 = Math.hypot(x, y) * rescale;
            // Use the variant of hypot to contrast with FMA.
            double o4 = Math.sqrt(x2y2Hypot(x, y)) * rescale;
            //double o4 = Complex.ofCartesian(x, y).abs() * rescale;

            // Only interested in differences from realistic best result without BigDecimal or
            // extended precision sqrt (although we could add this using Dekker's method)
            if (e == o1 && e == o2 && e == o3 && e == o4) {
                same++;
                return;
            }

            // Reset scale
            x *= rescale;
            y *= rescale;

            // Exact
            BigDecimal bdexact = x2y2BigDecimal(x, y).sqrt(MathContext.DECIMAL128);

            //score(x, y, o1, o2, o3, o4, bdexact);
            score(x, y, o3, o3, o3, o4, bdexact);
        }
    }

    private static class HypotChecker extends UlpChecker {
        HypotChecker(DoubleSupplier fx,
            DoubleSupplier fy, int samples) {
            super(fx, fy, samples);
        }

        /**
         * Check the difference between the computation of {@code x^2 + y^2} using
         * {@link Math#fma(double, double, double)} and
         * {@link Math#hypot(double, double)}. If the results differ the true result
         * is computed using high precision.BigDecimal and the ULPs computed from the
         * extended precision result.
         *
         * <p>This can be used to assert that the custom implementation of abs() is no worse than
         * {@link Math#hypot(double, double)} which aims to be within 1 ULP of the exact result.
         *
         * <p>Note: This method will not handle an input complex that is infinite or nan.
         */
        void check() {
            double x = Math.abs(fx.getAsDouble());
            double y = Math.abs(fy.getAsDouble());

            // Reference
            double o3 = Math.hypot(x, y);
            // Use the variant of hypot to contrast with FMA.
            double o4 = Complex.ofCartesian(x, y).abs();

            // Only interested in differences from realistic best result without BigDecimal or
            // extended precision sqrt (although we could add this using Dekker's method)
            //if (e == o1 && e == o2 && e == o3 && e == o4) {
            if (o3 == o4) {
                same++;
                return;
            }

            // Verify the exact result
            BigDecimal bdexact = x2y2BigDecimal(x, y).sqrt(MathContext.DECIMAL128);

            //score(x, y, o1, o2, o3, o4, bdexact);
            score(x, y, o3, o3, o3, o4, bdexact);
        }
    }

    private static class CisNumberGenerator implements DoubleSupplier {
        final UniformRandomProvider rng;
        double tmp = Double.NaN;

        CisNumberGenerator(UniformRandomProvider rng) {
            this.rng = rng;
        }

        @Override
        public double getAsDouble() {
            if (Double.isNaN(tmp)) {
                // Random angle
                double u = rng.nextDouble() * Math.PI;
                tmp = Math.cos(u);
                return Math.sin(u);
            }
            final double r = tmp;
            tmp = Double.NaN;
            return r;
        }
    }

    private static class PolarNumberGenerator extends CisNumberGenerator {
        PolarNumberGenerator(UniformRandomProvider rng) {
            super(rng);
        }

        @Override
        public double getAsDouble() {
            if (Double.isNaN(tmp)) {
                // Random angle and magnitude
                double u = rng.nextDouble() * Math.PI;
                double v = rng.nextDouble() * 2;
                tmp = Math.cos(u) * v;
                return Math.sin(u) * v;
            }
            final double r = tmp;
            tmp = Double.NaN;
            return r;
        }
    }

    private static class CisNumberGenerator2 extends CisNumberGenerator {
        final ZigguratNormalizedGaussianSampler sampler;

        CisNumberGenerator2(UniformRandomProvider rng) {
            super(rng);
            sampler = ZigguratNormalizedGaussianSampler.of(rng);
        }

        @Override
        public double getAsDouble() {
            if (Double.isNaN(tmp)) {
                double x;
                double y;
                double norm;
                do {
                    x = sampler.sample();
                    y = sampler.sample();
                    //norm = x * x + y * y;
                    // Do this as well as we can
                    norm = x2y2DekkerSqrt(x, y);
                } while (norm == 0);
                tmp = x / norm;
                return y / norm;
            }
            final double r = tmp;
            tmp = Double.NaN;
            return r;
        }
    }

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

    private static double hypotExact(double x, double y) {
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
        double x2y2 = x2y2Exact(a, b);
        return Math.sqrt(x2y2) * rescale;
    }

    private static double hypotStandard(double x, double y) {
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

        return Math.sqrt(x * x + y * y) * rescale;
        //return Math.sqrt(x2y2Hypot(a, b)) * rescale;
    }

    private static double hypot3(double x, double y) {
        // The mask is used to remove the sign.
        long bitsx = Double.doubleToRawLongBits(x) & 0x7fff_ffff_ffff_ffffL;
        long bitsy = Double.doubleToRawLongBits(y) & 0x7fff_ffff_ffff_ffffL;

        // Order by magnitude
        double a;
        double b;
        if (bitsx < bitsy) {
            a = y;
            b = x;
            final long tmp = bitsx;
            bitsx = bitsy;
            bitsy = tmp;
        } else {
            a = x;
            b = y;
        }

        // Compute absolutes
        a = Math.abs(a);
        b = Math.abs(b);

        // inf/nan handling
        if (bitsx >= 0x7ff0_0000_0000_0000L) {
            // a is inf/nan. Return inf is b is infinite
            return bitsy == 0x7ff0_0000_0000_0000L ? b : a;
        }

        // Finite numbers
        // Do scaling towards 1 but avoid x^2+y^2 being too close to 1.
        // This forces the result of sqrt to be a different order of magnitude.
        // Dropping the last 3 bits limits the exponent to [-1023, -1016], [+8, +1017] (unbiased)
        long exponent = bitsx & 0x7f80_0000_0000_0000L;
        double scale = Double.longBitsToDouble(0x7fc0_0000_0000_0000L - exponent);
        double rescale = Double.longBitsToDouble(0x0020_0000_0000_0000L + exponent);

        a *= scale;
        b *= scale;

        // a and b in [2^-52, 2^6)
        return rescale * Math.sqrt(fma(a, b));
    }

    private static double x2y2Hypot(double x, double y) {
        // Below we use a Dekker split to create two 26-bit numbers avoiding use of
        // conversion to raw bits. The original fdlibm version used a split of
        // the upper 32-bits and lower 32-bits. This results in the upper part
        // being a 21-bit precision number (including assumed leading 1) and the
        // The different split can create 1 ULP differences
        // in the result of Math.sqrt(x2y2(x, y)) verses the JDK Math.hypot(x, y)
        // implementation of the fdlibm function. Replacing the Dekker split with
        // Double.longBitsToDouble(Double.doubleToRawLongBits(a) & 0xffffffff00000000L)
        // computes the same result as Math.hypot(x, y).
        // Testing with billions of numbers show that either split achieves the target of
        // 1 ULP from the result of a high precision computation
        // (e.g. see o.a.c.numbers.arrays.LinearCombination) and standard Math.sqrt().
        // Performance comparisons show the Dekker split is faster.
        //
        // This could be replaced by Math.fma(x, x, y * y).
        // Tests with FMA on JDK 9 show precision is not significantly changed from this
        // implementation. Using simulated fused multiply addition (FMA)
        // with split precision summation is slower. A full high precision computation
        // is slower again and is a poor trade-off to gain 1 ULP accuracy.
        final double w = x - y;
        if (w > y) {
            final double t1 = Double.longBitsToDouble(Double.doubleToRawLongBits(x) & 0xffffffff00000000L);
            //final double t1 = splitHigh(x);
            final double t2 = x - t1;
            return t1 * t1 - (y * (-y) - t2 * (x + t1));
        }
        // 2y > x > y
        final double t = x + x;
        final double y1 = Double.longBitsToDouble(Double.doubleToRawLongBits(y) & 0xffffffff00000000L);
        //final double y1 = splitHigh(y);
        final double y2 = y - y1;
        final double t1 = Double.longBitsToDouble(Double.doubleToRawLongBits(t) & 0xffffffff00000000L);
        //final double t1 = splitHigh(t);
        final double t2 = t - t1;
        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));
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
     * Then modifies the exponent by the given amount.
     *
     * @param rng Source of randomness
     * @param exponent Amount to change the exponent (in range [-1023, 1023])
     * @return the number
     */
    private static double createFixedExponentNumber(UniformRandomProvider rng, int exponent) {
        return Double.longBitsToDouble((rng.nextLong() >>> 12) | ((1023L + exponent) << 52));
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
        assertAbs(rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 0), samples);
        report("range0");
        assertAbs(rng, r -> createFixedExponentNumber(r, 0),
            r -> createFixedExponentNumber(r, 1), samples);
        report("range1");
        assertAbs(rng, r -> createFixedExponentNumber(r, 0),
            r -> createFixedExponentNumber(r, 2), samples);
        report("range2");
        assertAbs(rng, r -> createFixedExponentNumber(r, 0),
            r -> createFixedExponentNumber(r, 3), samples);
        report("range3");
        assertAbs(rng, r -> createFixedExponentNumber(r, 0),
            r -> createFixedExponentNumber(r, 4), samples);
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
        final ToDoubleFunction<UniformRandomProvider> luGenerator = new ToDoubleFunction<>() {
            @Override
            public double applyAsDouble(UniformRandomProvider rng) {
                return createFixedExponentNumber(rng, rng.nextInt(30));
            }
        };
        assertAbs(rng, luGenerator, luGenerator, samples);
        report("log-uniform");
        // Complex cis numbers
        final ToDoubleFunction<UniformRandomProvider> cisGenerator = new ToDoubleFunction<>() {
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
    public void testFma() throws InterruptedException, ExecutionException {
        final JumpableUniformRandomProvider rng =
            (JumpableUniformRandomProvider) RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP, 526735472636L);

        // Can run max of approximately 1 billion per run.
        int samples = 1 << 16;
        int runs = 8;

        // Sub-normals are slow to compute the correct answer
        int subNormalSamples = samples >>> 0;

        ExecutorService es = Executors.newFixedThreadPool(
            Math.min(Runtime.getRuntime().availableProcessors(), runs));
        try {
            // Cis and uniform are2 slower. Is this just because we have to use BigDecimal
            // to check the answer more often. Or use sin/cos to create a cis number?
            // These are cases where the branch prediction
            // in hypot cannot learn which to choose.
            checkFma("cis", rng, CisNumberGenerator::new, samples, runs, es);
            checkFma("cis2", rng, CisNumberGenerator2::new, samples, runs, es);
            checkFma("polar", rng, PolarNumberGenerator::new, samples, runs, es);
            checkFma("uniform", rng, UniformRandomProvider::nextDouble, UniformRandomProvider::nextDouble, samples, runs, es);
            checkFma("random", rng, r -> {
                ZigguratNormalizedGaussianSampler s = ZigguratNormalizedGaussianSampler.of(r);
                return s::sample;
            }, samples, runs, es);

            checkFma("range_p0_p0", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 0), samples, runs, es);
            checkFma("range_p0_p1", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 1), samples, runs, es);
            checkFma("range_p0_p2", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 2), samples, runs, es);
            checkFma("range_p0_p3", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 3), samples, runs, es);
            checkFma("range_p0_p4", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 4), samples, runs, es);

            checkFma("range_m1_m1", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -1), samples, runs, es);
            checkFma("range_m1_m2", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -2), samples, runs, es);
            checkFma("range_m1_m3", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -3), samples, runs, es);
            checkFma("range_m1_m4", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -4), samples, runs, es);
            checkFma("range_m1_m5", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -5), samples, runs, es);

            checkFma("range_p128_p128", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 128), samples, runs, es);
            checkFma("range_p128_p129", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 129), samples, runs, es);
            checkFma("range_p128_p130", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 130), samples, runs, es);
            checkFma("range_p128_p131", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 131), samples, runs, es);
            checkFma("range_p128_p132", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 132), samples, runs, es);

            checkFma("range_m128_m128", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -128), samples, runs, es);
            checkFma("range_m128_m129", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -129), samples, runs, es);
            checkFma("range_m128_m130", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -130), samples, runs, es);
            checkFma("range_m128_m131", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -131), samples, runs, es);
            checkFma("range_m128_m132", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -132), samples, runs, es);

            // Cannot check these due to scaling. Maybe add a scaled version for BigDecimal.
            // Or just state that ULPs have less meaning for sub-normals.

            checkFmaScaled("m1010 sub-52", rng, r -> createFixedExponentNumber(r, -1010), CStandardTest2::createSubNormalNumber52, samples, runs, es);
            checkFmaScaled("m1021 sub-52", rng, r -> createFixedExponentNumber(r, -1022), CStandardTest2::createSubNormalNumber52, subNormalSamples, runs, es);
            checkFmaScaled("m1022 sub-52", rng, r -> createFixedExponentNumber(r, -1022), CStandardTest2::createSubNormalNumber52, subNormalSamples, runs, es);
            checkFmaScaled("sub-52 sub-52", rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber52, subNormalSamples, runs, es);
            checkFmaScaled("sub-52 sub-32", rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber32, subNormalSamples, runs, es);
            // This is not hard so omit
            checkFmaScaled("sub-normal 32/32", rng, CStandardTest2::createSubNormalNumber32, CStandardTest2::createSubNormalNumber32, subNormalSamples, runs, es);

            // The distribution of floating-point numbers has the log-uniform distribution as
            // its limiting distribution. We do not use that to force a full 52-bit random mantissa
            // which increases the likelihood for high precision.
            // This pointless as a big difference in magnitude is easy.
            final ToDoubleFunction<UniformRandomProvider> luGenerator = new ToDoubleFunction<>() {
                @Override
                public double applyAsDouble(UniformRandomProvider rng) {
                    // Do not use scalb. Directly write the exponent.
                    return createFixedExponentNumber(rng, rng.nextInt(30));
                }
            };
            checkFma("log-uniform", rng, luGenerator, luGenerator, samples, runs, es);
        } finally {
            es.shutdown();
        }
    }

    private static void checkFma(String type, JumpableUniformRandomProvider rng,
        ToDoubleFunction<UniformRandomProvider> fx,
        ToDoubleFunction<UniformRandomProvider> fy, int samples, int runs, ExecutorService es) throws InterruptedException, ExecutionException {
        ArrayList<Future<UlpChecker>> list = new ArrayList<>();
        for (int i = 0; i < runs; i++) {
            UniformRandomProvider r = rng.jump();
            DoubleSupplier f1 = () -> fx.applyAsDouble(r);
            DoubleSupplier f2 = () -> fy.applyAsDouble(r);
            list.add(es.submit(new UlpChecker(f1, f2, samples)));
        }
        report(type, list);
    }

    private static void checkFmaScaled(String type, JumpableUniformRandomProvider rng,
        ToDoubleFunction<UniformRandomProvider> fx,
        ToDoubleFunction<UniformRandomProvider> fy, int samples, int runs, ExecutorService es) throws InterruptedException, ExecutionException {
        ArrayList<Future<UlpChecker>> list = new ArrayList<>();
        for (int i = 0; i < runs; i++) {
            UniformRandomProvider r = rng.jump();
            DoubleSupplier f1 = () -> fx.applyAsDouble(r);
            DoubleSupplier f2 = () -> fy.applyAsDouble(r);
            list.add(es.submit(new ScaledUlpChecker(f1, f2, samples)));
        }
        report(type, list);
    }

    private static void checkFma(String type, JumpableUniformRandomProvider rng,
        Function<UniformRandomProvider, DoubleSupplier> f, int samples, int runs, ExecutorService es) throws InterruptedException, ExecutionException {
        ArrayList<Future<UlpChecker>> list = new ArrayList<>();
        for (int i = 0; i < runs; i++) {
            DoubleSupplier fx = f.apply(rng.jump());
            list.add(es.submit(new UlpChecker(fx, fx, samples)));
        }
        report(type, list);
    }

    private static void report(String type, ArrayList<Future<UlpChecker>> list) throws InterruptedException, ExecutionException {
        UlpChecker r = null;
        for (Future<UlpChecker> f : list) {
            r = f.get().add(r);
        }
        if (r == null) {
            return;
        }

        // CHECKSTYLE: stop all
        long tot = (long) r.samples * list.size();
        double t = tot * 100000.0;
        System.out.printf("%-17s  same %d / %d (%.5f) : dekker %.5f (%.5f) %10d (%.5f) : fma %.5f (%.5f) %10d (%.5f) : hypot %.5f (%.5f) %10d (%.5f) : hypot2 %.5f (%.5f) %10d (%.5f) : better %d, worse %d : %s : %s : %s : %s%n",
                type,
                r.same, tot, (double) r.same / tot,
                r.s1 / t, r.m1, r.d1, (double) (tot - r.d1) / tot,
                r.s2 / t, r.m2, r.d2, (double) (tot - r.d2) / tot,
                r.s3 / t, r.m3, r.d3, (double) (tot - r.d3) / tot,
                r.s4 / t, r.m4, r.d4, (double) (tot - r.d4) / tot,
                r.better, r.worse,
                r.z1, r.z2, r.z3, r.z4);
        // CHECKSTYLE: resume all
    }

    @Test
    public void testHypot() throws InterruptedException, ExecutionException {
        final JumpableUniformRandomProvider rng =
            (JumpableUniformRandomProvider) RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP, 526735472636L);

        // Can run max of approximately 1 billion per run.
        int samples = 1 << 16;
        int runs = 8;

        ExecutorService es = Executors.newFixedThreadPool(
            Math.min(Runtime.getRuntime().availableProcessors(), runs));
        try {
            // Cis and uniform are2 slower. Is this just because we have to use BigDecimal
            // to check the answer more often. Or use sin/cos to create a cis number?
            // These are cases where the branch prediction
            // in hypot cannot learn which to choose.
            checkHypot("cis", rng, CisNumberGenerator::new, samples, runs, es);
            checkHypot("cis2", rng, CisNumberGenerator2::new, samples, runs, es);
            checkHypot("polar", rng, PolarNumberGenerator::new, samples, runs, es);
            checkHypot("uniform", rng, UniformRandomProvider::nextDouble, UniformRandomProvider::nextDouble, samples, runs, es);
            checkHypot("random", rng, r -> {
                ZigguratNormalizedGaussianSampler s = ZigguratNormalizedGaussianSampler.of(r);
                return s::sample;
            }, samples, runs, es);

            checkHypot("range_p0_p0", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 0), samples, runs, es);
            checkHypot("range_p0_p1", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 1), samples, runs, es);
            checkHypot("range_p0_p2", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 2), samples, runs, es);
            checkHypot("range_p0_p3", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 3), samples, runs, es);
            checkHypot("range_p0_p4", rng, r -> createFixedExponentNumber(r, 0), r -> createFixedExponentNumber(r, 4), samples, runs, es);

            checkHypot("range_m1_m1", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -1), samples, runs, es);
            checkHypot("range_m1_m2", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -2), samples, runs, es);
            checkHypot("range_m1_m3", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -3), samples, runs, es);
            checkHypot("range_m1_m4", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -4), samples, runs, es);
            checkHypot("range_m1_m5", rng, r -> createFixedExponentNumber(r, -1), r -> createFixedExponentNumber(r, -5), samples, runs, es);

            checkHypot("range_p128_p128", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 128), samples, runs, es);
            checkHypot("range_p128_p129", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 129), samples, runs, es);
            checkHypot("range_p128_p130", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 130), samples, runs, es);
            checkHypot("range_p128_p131", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 131), samples, runs, es);
            checkHypot("range_p128_p132", rng, r -> createFixedExponentNumber(r, 128), r -> createFixedExponentNumber(r, 132), samples, runs, es);

            checkHypot("range_m128_m128", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -128), samples, runs, es);
            checkHypot("range_m128_m129", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -129), samples, runs, es);
            checkHypot("range_m128_m130", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -130), samples, runs, es);
            checkHypot("range_m128_m131", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -131), samples, runs, es);
            checkHypot("range_m128_m132", rng, r -> createFixedExponentNumber(r, -128), r -> createFixedExponentNumber(r, -132), samples, runs, es);

            // Cannot check these due to scaling. Maybe add a scaled version for BigDecimal.
            // Or just state that ULPs have less meaning for sub-normals.

            checkHypot("m1010 sub-52", rng, r -> createFixedExponentNumber(r, -1010), CStandardTest2::createSubNormalNumber52, samples, runs, es);
            checkHypot("m1021 sub-52", rng, r -> createFixedExponentNumber(r, -1022), CStandardTest2::createSubNormalNumber52, samples, runs, es);
            checkHypot("m1022 sub-52", rng, r -> createFixedExponentNumber(r, -1022), CStandardTest2::createSubNormalNumber52, samples, runs, es);
            checkHypot("sub-52 sub-52", rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber52, samples, runs, es);
            checkHypot("sub-52 sub-32", rng, CStandardTest2::createSubNormalNumber52, CStandardTest2::createSubNormalNumber32, samples, runs, es);
            // This is not hard so omit
            checkHypot("sub-normal 32/32", rng, CStandardTest2::createSubNormalNumber32, CStandardTest2::createSubNormalNumber32, samples, runs, es);

            // The distribution of floating-point numbers has the log-uniform distribution as
            // its limiting distribution. We do not use that to force a full 52-bit random mantissa
            // which increases the likelihood for high precision.
            // This pointless as a big difference in magnitude is easy.
            final ToDoubleFunction<UniformRandomProvider> luGenerator = new ToDoubleFunction<>() {
                @Override
                public double applyAsDouble(UniformRandomProvider rng) {
                    // Do not use scalb. Directly write the exponent.
                    return createFixedExponentNumber(rng, rng.nextInt(30));
                }
            };
            checkHypot("log-uniform", rng, luGenerator, luGenerator, samples, runs, es);
        } finally {
            es.shutdown();
        }
    }

    private static void checkHypot(String type, JumpableUniformRandomProvider rng,
        ToDoubleFunction<UniformRandomProvider> fx,
        ToDoubleFunction<UniformRandomProvider> fy, int samples, int runs, ExecutorService es) throws InterruptedException, ExecutionException {
        ArrayList<Future<UlpChecker>> list = new ArrayList<>();
        for (int i = 0; i < runs; i++) {
            UniformRandomProvider r = rng.jump();
            DoubleSupplier f1 = () -> fx.applyAsDouble(r);
            DoubleSupplier f2 = () -> fy.applyAsDouble(r);
            list.add(es.submit(new HypotChecker(f1, f2, samples)));
        }
        report(type, list);
    }

    private static void checkHypot(String type, JumpableUniformRandomProvider rng,
        Function<UniformRandomProvider, DoubleSupplier> f, int samples, int runs, ExecutorService es) throws InterruptedException, ExecutionException {
        ArrayList<Future<UlpChecker>> list = new ArrayList<>();
        for (int i = 0; i < runs; i++) {
            DoubleSupplier fx = f.apply(rng.jump());
            list.add(es.submit(new HypotChecker(fx, fx, samples)));
        }
        report(type, list);
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
     * High accuracy {@code x^2 + y^2} when {@code 2y > x > y}.
     *
     * @param x the x
     * @param y the y
     * @return {@code x^2 + y^2}
     */
    private static double x2y2TwoFactor(double x, double y) {
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
     * High accuracy {@code sqrt(x^2 + y^2)}.
     *
     * @param x the x
     * @param y the y
     * @return {@code x^2 + y^2}
     */
    private static double x2y2DekkerSqrt(double x, double y) {
        final double xx = x * x;
        final double yy = y * y;
        final double xHigh = splitHigh(x);
        final double xLow = x - xHigh;
        final double yHigh = splitHigh(y);
        final double yLow = y - yHigh;
        final double x2Low = squareLow(xLow, xHigh, xx);
        final double y2Low = squareLow(yLow, yHigh, yy);

        double r = xx + yy;
        double s = xx - r + yy + y2Low + x2Low;
        double z = r + s;
        double zz = r - z + s;

        // High precision sqrt
        final double c = Math.sqrt(z);
        final double cHigh = splitHigh(c);
        final double cLow = c - cHigh;
        final double u = c * c;
        final double uu = squareLow(cLow, cHigh, u);
        final double cc = (z - u - uu + zz) * 0.5 / c;
        r = c + cc;
        return c - r + cc + r;
    }

    /**
     * Compute the sum of the products of two sequences of factors.
     *
     * @param a1 First factor of the first term.
     * @param b1 Second factor of the first term.
     * @return \( a1^2 + b1^2 \)
     */
    private static double x2y2Exact(double a1, double b1) {
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

    private static BigDecimal x2y2BigDecimal(double x, double y) {
        return new BigDecimal(x).pow(2).add(new BigDecimal(y).pow(2));
    }

    /**
     * High accuracy {@code x^2 + y^2}.
     *
     * @param x the x
     * @param y the y
     * @return {@code x^2 + y^2}
     */
    private static double fma(double x, double y) {
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
     * @return <code>low * low - (((product - high * high) - low * high) - high * low)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    private static double squareLow(double low, double high, double square) {
        // 4 add, 3 multiply (7)
        final double lh = low * high;
        return low * low - (((square - high * high) - lh) - lh);
    }

    /**
     * Compute the round-off from the square of a split number with {@code low} and
     * {@code high} components. Uses Dekker's algorithm for split multiplication modified
     * for a square product.
     *
     * @param low Low part of number.
     * @param high High part of number.
     * @param value The number.
     * @param square Square of the number.
     * @return <code>low * low - ((product - high * high) - 2 * low * high)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    private static double squareLow(double low, double high, double value, double square) {
        // This works for sub-normal numbers as the products are smaller so rounding discards
        // bits in most computations.
        // 3 add, 2 multiply (5)
        return ((high * high - square) + low * (high + value));
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

    /**
     * Compute the upper part of the number using the Dekker split.
     *
     * @param a the number
     * @return the upper part
     */
    private static double splitHigh2(double a) {
        // Use 2^s for a (s-1) bit lower number. Using s=33 should equal the masking split.
        final double c = 0x1.0p33 * a;
        return c - (c - a);
    }

    @Test
    public void testUniform() {
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
        final long total2 = 1L << 20;
        long diff = 0;
        for (long l = total2; l-- > 0;) {
            double a = rng.nextDouble();
            double b = rng.nextDouble();
            double x = Math.max(a, b);
            double y = Math.min(a, b);
            if (x < 2 * y) {
                diff++;
                // CHECKSTYLE: stop all
                //System.out.printf("{%s, %s}%n", x, y);
            }
        }
        System.out.printf("Lower diff %d / %d (%s)%n", diff, total2, (double) diff / total2);
    }

    @Test
    public void testBadHighPrecision() {
        testBadHighPrecision(0.7600300504761309, 0.6498879306259271);
        //testBadHighPrecision(-0.8272978928538699, 0.561763470225278);
        // sub-normal numbers
        //testBadHighPrecision(Double.longBitsToDouble(0x0000_0856_3784_6824_3842L),
        //                     Double.longBitsToDouble(0x0000_0863_7846_8243_6842L));
    }

    private static void testBadHighPrecision(double x, double y) {
        final double xx = x;
        final double yy = y;
        new UlpChecker(() -> xx, () -> yy, 1).call();

        x = Math.abs(x);
        y = Math.abs(y);
        if (x < y) {
            final double tmp = x;
            x = y;
            y = tmp;
        }
        double rescale = 1.0;
        if (Math.min(x, y) < 0x1.0p-500) {
            if (Math.min(x, y) < Double.MIN_NORMAL) {
                x *= 0x1.0p+1022;
                y *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                x *= 0x1.0p+600;
                y *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        } else if (Math.max(x, y) > 0x1.0p500) {
            x *= 0x1.0p-600;
            y *= 0x1.0p-600;
            rescale = 0x1.0p+600;
        }
        double v1 = x2y2Dekker(x, y);
        double v2 = x2y2BigDecimal(x, y).doubleValue();
        double v3 = x2y2Hypot(x, y);
        System.out.printf("%s, %s, %s", v1, v2, v3);
        System.out.printf(" => %s, %s, %s", v1 * rescale, v2 * rescale, v3 * rescale);
        v1 = Math.sqrt(v1) * rescale;
        v2 = x2y2BigDecimal(x, y).sqrt(MathContext.DECIMAL64).doubleValue() * rescale;
        v3 = Math.sqrt(v3) * rescale;
        double v4 = x2y2DekkerSqrt(x, y) * rescale;
        System.out.printf(" => %s, %s, %s, %s%n", v1, v2, v3, v4);
    }

    @Test
    public void testSquareLow() {
        // Test the faster fdlibm version for the round-off from the square.
        // Does this always match the Dekker version. It saves 2 flops (7 => 5)
        //final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP);
        final long total2 = 1L << 20;
        long diff = 0;
        long start = Double.doubleToRawLongBits(0.0);
        //long start = Double.doubleToRawLongBits(Math.nextDown(2.0));
        for (long l = total2; l-- > 0;) {
            //final double x = Double.longBitsToDouble(rng.nextLong() >>> 12);
            final double x = Double.longBitsToDouble(start++);
            final double xx = x * x;
            final double xHigh = splitHigh(x);
            final double xLow = x - xHigh;
            final double x2Low1 = squareLow(xLow, xHigh, xx);
            final double x2Low2 = squareLow(xLow, xHigh, x, xx);
            if (x2Low1 != x2Low2) {
                diff++;
            }
        }
        System.out.printf("%d / %d  = %s%n", diff, total2, (double) diff / total2);
    }

    @Test
    public void testSqrtLog() {
        //testSqrtLog("Dekker",  CStandardTest2::x2y2Dekker);
        testSqrtLog("hypot2",  CStandardTest2::x2y2Hypot);
        //testSqrtLog("fmaSim",  CStandardTest2::fma);
        //testSqrtLog("fma",  (x, y) -> Math.fma(x, x, y * y));
        //testSqrtLog("standard", (x, y) -> x * x + y * y);
    }

    private static void testSqrtLog(String name, DoubleDoubleBiFunction f) {
        // Find cases where hypot is different and see if the sqrt/log computation is different
        // at different scales.
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP, 567464367L);
        final long total2 = 1L << 30;
        long diff = 0;
        long diffSqrt = 0;
        long maxSqrt = 0;
        long diffLog1 = 0;
        long diffLog2 = 0;
        long diffLog3 = 0;
        long maxLog1 = 0;
        long maxLog2 = 0;
        long maxLog3 = 0;
        long d;
        for (long l = total2; l-- > 0;) {
            double x;
            double y;
            x = createFixedExponentNumber(rng, 2);
            y = createFixedExponentNumber(rng, 2);
            //double u = rng.nextDouble() * Math.PI;
            //x = Math.abs(Math.sin(u));
            //y = Math.cos(u);
            if (x < y) {
                final double tmp = x;
                x = y;
                y = tmp;
            }
            double s = f.apply(x, y);
            double s1 = Math.sqrt(s);
            double s2 = Math.hypot(x, y);
            if (s1 != s2) {
                diff++;
                // CHECKSTYLE: stop all
                //System.out.printf("{%s, %s}%n", x, y);
                double r1 = Math.sqrt(2 * (s1 + x));
                double r2 = Math.sqrt(2 * (s2 + x));
                if (r1 != r2) {
                    diffSqrt++;
                    d = Math.abs(Double.doubleToRawLongBits(r1) - Double.doubleToRawLongBits(r2));
                    if (maxSqrt < d) {
                        maxSqrt = d;
                    }
                }
                // Cannot use log when input is between 0.5 and 2.0
                double l1 = 0.5 * Math.log(s);
                double l2 = Math.log(s2);
                // alternative
                double l3 = Math.log(s1);
                if (l1 != l2) {
                    diffLog1++;
                    d = Math.abs(Double.doubleToRawLongBits(l1) - Double.doubleToRawLongBits(l2));
                    if (maxLog1 < d) {
                        maxLog1 = d;
                    }
                }
                if (l1 != l3) {
                    diffLog2++;
                    d = Math.abs(Double.doubleToRawLongBits(l1) - Double.doubleToRawLongBits(l3));
                    if (maxLog2 < d) {
                        maxLog2 = d;
                    }
                }
                if (l2 != l3) {
                    diffLog3++;
                    d = Math.abs(Double.doubleToRawLongBits(l1) - Double.doubleToRawLongBits(l3));
                    if (maxLog3 < d) {
                        maxLog3 = d;
                    }
                }
            }
        }
        System.out.printf("%-10s  diff %10d / %d (%.5f) : sqrt %10d (%.5f) (%3d) : log1 %10d (%.5f) (%3d) : log2 %10d (%.5f) (%3d) : log3 %10d (%.5f) (%3d)%n",
                name, diff, total2, (double) diff / total2,
                diffSqrt, (double) diffSqrt / total2, maxSqrt,
                diffLog1, (double) diffLog1 / total2, maxLog1,
                diffLog2, (double) diffLog2 / total2, maxLog2,
                diffLog3, (double) diffLog3 / total2, maxLog3
                );
    }

    @Test
    public void testMultiplyDistribution() throws IOException {
        // at different scales.
        final UniformRandomProvider rng = RandomSource.create(RandomSource.XO_RO_SHI_RO_128_PP, 567464367L);
        final long total2 = 1L << 16;
        try (BufferedWriter out = Files.newBufferedWriter(Paths.get("/tmp/1.dat"));
            Formatter f = new Formatter(out)) {
            ZigguratNormalizedGaussianSampler s = ZigguratNormalizedGaussianSampler.of(rng);
            for (long l = total2; l-- > 0;) {
                //Complex c1 = Complex.ofPolar(rng.nextDouble() * 10, rng.nextDouble() * Math.PI);
                //Complex c2 = Complex.ofPolar(rng.nextDouble() * 10, rng.nextDouble() * Math.PI);
                Complex c1 = Complex.ofCartesian(s.sample(), s.sample());
                Complex c2 = Complex.ofCartesian(s.sample(), s.sample());
                //Complex z = Complex.ofCartesian(s.sample(), s.sample());
                Complex z = c1.multiply(c2);
                double x = Math.abs(z.real());
                double y = Math.abs(z.imag());
                // Determine if 2y > x > y > 0
                if (x < y) {
                    final double tmp = x;
                    x = y;
                    y = tmp;
                }
                double w = x - y;
                f.format("%s %s %s %s %s %d%n", x, y, z.abs(), Math.log(x), Math.log(y), (w > y) ? 1 : 0);
            }
        }
    }
}
