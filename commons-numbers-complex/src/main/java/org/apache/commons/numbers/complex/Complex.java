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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Cartesian representation of a complex number, i.e. a number which has both a
 * real and imaginary part.
 *
 * <p>This class is immutable. All arithmetic will create a new instance for the
 * result.</p>
 *
 * <p>Arithmetic in this class conforms to the C99 standard for complex numbers
 * defined in ISO/IEC 9899, Annex G. Methods have been named using the equivalent
 * method in ISO C99. The behaviour for special cases is listed as defined in C99.</p>
 *
 * <p>For functions \( f \) which obey the conjugate equality \( conj(f(z)) = f(conj(z)) \),
 * the specifications for the upper half-plane imply the specifications for the lower
 * half-plane.</p>
 *
 * <p>For functions that are either odd, \( f(z) = -f(-z) \), or even, \( f(z) =  f(-z) \)
 * the specifications for the first quadrant imply the specifications for the other three
 * quadrants.</p>
 *
 * <p>Special cases of <a href="http://mathworld.wolfram.com/BranchCut.html">branch cuts</a>
 * for multivalued functions adopt the principle value convention from C99. Specials cases
 * from C99 that raise the "invalid" or "divide-by-zero"
 * <a href="https://en.cppreference.com/w/c/numeric/fenv/FE_exceptions">floating-point
 * exceptions</a> return the documented value without an explicit mechanism to notify
 * of the exception case, that is no exceptions are thrown during computations in-line with
 * the convention of the corresponding single-valued functions in {@link java.lang.Math}.
 * These cases are documented in the method special cases as "invalid" or "divide-by-zero"
 * floating-point operation.
 * Note: Invalid floating-point exception cases will result in a complex number where the
 * cardinality of NaN component parts has increased as a real or imaginary part could
 * not be computed and is set to NaN.
 *
 * @see <a href="http://www.open-std.org/JTC1/SC22/WG14/www/standards">
 *    ISO/IEC 9899 - Programming languages - C</a>
 */
public final class Complex implements Serializable  {
    /**
     * A complex number representing \( i \), the square root of \( -1 \).
     *
     * <p>\( (0 + i 1) \).
     */
    public static final Complex I = new Complex(0, 1);
    /**
     * A complex number representing one.
     *
     * <p>\( (1 + i 0) \).
     */
    public static final Complex ONE = new Complex(1, 0);
    /**
     * A complex number representing zero.
     *
     * <p>\( (0 + i 0) \).
     */
    public static final Complex ZERO = new Complex(0, 0);

    /** A complex number representing {@code NaN + i NaN}. */
    private static final Complex NAN = new Complex(Double.NaN, Double.NaN);
    /** &pi;/2. */
    private static final double PI_OVER_2 = 0.5 * Math.PI;
    /** &pi;/4. */
    private static final double PI_OVER_4 = 0.25 * Math.PI;
    /** Mask an integer number to even by discarding the lowest bit. */
    private static final int MASK_INT_TO_EVEN = ~0x1;
    /** Natural logarithm of 2 (ln(2)). */
    private static final double LN_2 = Math.log(2);
    /** Base 10 logarithm of 10 divided by 2 (log10(e)/2). */
    private static final double LOG_10E_O_2 = Math.log10(Math.E) / 2;
    /** Base 10 logarithm of 2 (log10(2)). */
    private static final double LOG10_2 = Math.log10(2);
    /** {@code 1/2}. */
    private static final double HALF = 0.5;
    /** {@code sqrt(2)}. */
    private static final double ROOT2 = Math.sqrt(2);
    /** Half the number of bits of precision of the mantissa of a {@code double} rounded up: {@code 27}. */
    private static final double HALF_PRECISION = 27;
    /** The bit representation of {@code -0.0}. */
    private static final long NEGATIVE_ZERO_LONG_BITS = Double.doubleToLongBits(-0.0);
    /** Exponent offset in IEEE754 representation. */
    private static final long EXPONENT_OFFSET = 1023L;
    /**
     * Largest double-precision floating-point number such that
     * {@code 1 + EPSILON} is numerically equal to 1. This value is an upper
     * bound on the relative error due to rounding real numbers to double
     * precision floating-point numbers.
     *
     * <p>In IEEE 754 arithmetic, this is 2<sup>-53</sup>.
     * Copied from o.a.c.numbers.Precision.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Machine_epsilon">Machine epsilon</a>
     */
    private static final double EPSILON = Double.longBitsToDouble((EXPONENT_OFFSET - 53L) << 52);
    /** Mask to remove the sign bit from a long. */
    private static final long UNSIGN_MASK = 0x7fffffffffffffffL;
    /** Mask to remove the sign bit from a long. */
    private static final int UNSIGN_INT_MASK = 0x7fffffff;
    /** Mask to extract the 52-bit mantissa from a long representation of a double. */
    private static final long MANTISSA_MASK = 0x000fffffffffffffL;
    /** The multiplier used to split the double value into hi and low parts. This must be odd
     * and a value of 2^s + 1 in the range {@code p/2 <= s <= p-1} where p is the number of
     * bits of precision of the floating point number. Here {@code s = 27}.*/
    private static final double MULTIPLIER = 1.34217729E8;

    /**
     * Crossover point to switch computation for asin/acos factor A.
     * This has been updated from the 1.5 value used by Hull et al to 10
     * as used in boost::math::complex.
     * @see <a href="https://svn.boost.org/trac/boost/ticket/7290">Boost ticket 7290</a>
     */
    private static final double A_CROSSOVER = 10;
    /** Crossover point to switch computation for asin/acos factor B. */
    private static final double B_CROSSOVER = 0.6471;
    /**
     * The safe maximum double value {@code x} to avoid loss of precision in asin/acos.
     * Equal to sqrt(M) / 8 in Hull, et al (1997) with M the largest normalised floating-point value.
     */
    private static final double SAFE_MAX = Math.sqrt(Double.MAX_VALUE) / 8;
    /**
     * The safe minimum double value {@code x} to avoid loss of precision/underflow in asin/acos.
     * Equal to sqrt(u) * 4 in Hull, et al (1997) with u the smallest normalised floating-point value.
     */
    private static final double SAFE_MIN = Math.sqrt(Double.MIN_NORMAL) * 4;
    /**
     * The safe maximum double value {@code x} to avoid loss of precision in atanh.
     * Equal to sqrt(M) / 2 with M the largest normalised floating-point value.
     */
    private static final double SAFE_UPPER = Math.sqrt(Double.MAX_VALUE) / 2;
    /**
     * The safe minimum double value {@code x} to avoid loss of precision/underflow in atanh.
     * Equal to sqrt(u) * 2 with u the smallest normalised floating-point value.
     */
    private static final double SAFE_LOWER = Math.sqrt(Double.MIN_NORMAL) * 2;
    /**
     * A safe maximum double value {@code m} where {@code e^m} is not infinite.
     * This can be used when functions require approximations of sinh(x) or cosh(x)
     * when x is large using exp(x):
     * <pre>
     * sinh(x) = (e^x - e^-x) / 2 = sign(x) * e^|x| / 2
     * cosh(x) = (e^x + e^-x) / 2 = e^|x| / 2 </pre>
     *
     * <p>This value can be used to approximate e^x using a product:
     *
     * <pre>
     * e^x = product_n (e^m) * e^(x-nm)
     * n = (int) x/m
     * e.g. e^2000 = e^m * e^m * e^(2000 - 2m) </pre>
     *
     * <p>The value should be below ln(max_value) ~ 709.783.
     * The value m is set to an integer for less error when subtracting m and chosen as
     * even (m=708) as it is used as a threshold in tanh with m/2.
     *
     * <p>The value is used to compute e^x multiplied by a small number avoiding
     * overflow (sinh/cosh) or a small number divided by e^x without underflow due to
     * infinite e^x (tanh). The following conditions are used:
     * <pre>
     * 0.5 * e^m * Double.MIN_VALUE * e^m * e^m = Infinity
     * 2.0 / e^m / e^m = 0.0 </pre>
     */
    private static final double SAFE_EXP = 708;
    /**
     * The value of Math.exp(SAFE_EXP): e^708.
     * To be used in overflow/underflow safe products of e^m to approximate e^x where x > m.
     */
    private static final double EXP_M = Math.exp(SAFE_EXP);

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

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20180201L;

    /**
     * The size of the buffer for {@link #toString()}.
     *
     * <p>The longest double will require a sign, a maximum of 17 digits, the decimal place
     * and the exponent, e.g. for max value this is 24 chars: -1.7976931348623157e+308.
     * Set the buffer size to twice this and round up to a power of 2 thus
     * allowing for formatting characters. The size is 64.
     */
    private static final int TO_STRING_SIZE = 64;
    /** The minimum number of characters in the format. This is 5, e.g. {@code "(0,0)"}. */
    private static final int FORMAT_MIN_LEN = 5;
    /** {@link #toString() String representation}. */
    private static final char FORMAT_START = '(';
    /** {@link #toString() String representation}. */
    private static final char FORMAT_END = ')';
    /** {@link #toString() String representation}. */
    private static final char FORMAT_SEP = ',';
    /** The minimum number of characters before the separator. This is 2, e.g. {@code "(0"}. */
    private static final int BEFORE_SEP = 2;

    /** The imaginary part. */
    private final double imaginary;
    /** The real part. */
    private final double real;

    /**
     * Define a constructor for a Complex.
     * This is used in functions that implement trigonomic identities.
     */
    @FunctionalInterface
    private interface ComplexConstructor {
        /**
         * Create a complex number given the real and imaginary parts.
         *
         * @param real Real part.
         * @param imaginary Imaginary part.
         * @return {@code Complex} object.
         */
        Complex create(double real, double imaginary);
    }

    /**
     * Define a unary operation on a double.
     * This is used in the log() and log10() functions.
     */
    @FunctionalInterface
    private interface UnaryOperation {
        /**
         * Apply an operation to a value.
         *
         * @param value The value.
         * @return The result.
         */
        double apply(double value);
    }

    /**
     * Private default constructor.
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     */
    private Complex(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @return {@code Complex} number.
     */
    public static Complex ofCartesian(double real, double imaginary) {
        return new Complex(real, imaginary);
    }

    /**
     * Creates a complex number from its polar representation using modulus {@code rho} (\( \rho \))
     * and phase angle {@code theta} (\( \theta \)).
     *
     * \[ x = \rho \cos(\theta) \\
     *    y = \rho \sin(\theta) \]
     *
     * <p>Requires that {@code rho} is non-negative and non-NaN and {@code theta} is finite;
     * otherwise returns a complex with NaN real and imaginary parts. A value of {@code -0.0} is
     * considered negative and an invalid modulus.
     *
     * <p>A non-NaN complex number constructed using this method will satisfy the following
     * to within floating-point error when {@code theta} is in the range
     * \( -\pi\ \lt \theta \leq \pi \):
     *
     * <pre>
     *  Complex.ofPolar(rho, theta).abs() == rho
     *  Complex.ofPolar(rho, theta).arg() == theta</pre>
     *
     * <p>If {@code rho} is infinite then the resulting parts may be infinite or NaN
     * following the rules for double arithmetic, for example:</p>
     *
     * <ul>
     * <li>{@code ofPolar(}\( -0.0 \){@code , }\( 0 \){@code ) = }\( \text{NaN} + i \text{NaN} \)
     * <li>{@code ofPolar(}\( 0.0 \){@code , }\( 0 \){@code ) = }\( 0 + i 0 \)
     * <li>{@code ofPolar(}\( 1 \){@code , }\( 0 \){@code ) = }\( 1 + i 0 \)
     * <li>{@code ofPolar(}\( 1 \){@code , }\( \pi \){@code ) = }\( -1 + i \sin(\pi) \)
     * <li>{@code ofPolar(}\( \infty \){@code , }\( \pi \){@code ) = }\( -\infty + i \infty \)
     * <li>{@code ofPolar(}\( \infty \){@code , }\( 0 \){@code ) = }\( -\infty + i \text{NaN} \)
     * <li>{@code ofPolar(}\( \infty \){@code , }\( -\frac{\pi}{4} \){@code ) = }\( \infty - i \infty \)
     * <li>{@code ofPolar(}\( \infty \){@code , }\( 5\frac{\pi}{4} \){@code ) = }\( -\infty - i \infty \)
     * </ul>
     *
     * <p>This method is the functional equivalent of the C++ method {@code std::polar}.
     *
     * @param rho The modulus of the complex number.
     * @param theta The argument of the complex number.
     * @return {@code Complex} number.
     * @see <a href="http://mathworld.wolfram.com/PolarCoordinates.html">Polar Coordinates</a>
     */
    public static Complex ofPolar(double rho, double theta) {
        // Require finite theta and non-negative, non-nan rho
        if (!Double.isFinite(theta) || negative(rho) || Double.isNaN(rho)) {
            return NAN;
        }
        final double x = rho * Math.cos(theta);
        final double y = rho * Math.sin(theta);
        return new Complex(x, y);
    }

    /**
     * Create a complex cis number. This is also known as the complex exponential:
     *
     * \[ \text{cis}(x) = e^{ix} = \cos(x) + i \sin(x) \]
     *
     * @param x {@code double} to build the cis number.
     * @return {@code Complex} cis number.
     * @see <a href="http://mathworld.wolfram.com/Cis.html">Cis</a>
     */
    public static Complex ofCis(double x) {
        return new Complex(Math.cos(x), Math.sin(x));
    }

    /**
     * Returns a {@code Complex} instance representing the specified string {@code s}.
     *
     * <p>If {@code s} is {@code null}, then a {@code NullPointerException} is thrown.
     *
     * <p>The string must be in a format compatible with that produced by
     * {@link #toString() Complex.toString()}.
     * The format expects a start and end parentheses surrounding two numeric parts split
     * by a separator. Leading and trailing spaces are allowed around each numeric part.
     * Each numeric part is parsed using {@link Double#parseDouble(String)}. The parts
     * are interpreted as the real and imaginary parts of the complex number.
     *
     * <p>Examples of valid strings and the equivalent {@code Complex} are shown below:
     *
     * <pre>
     * "(0,0)"             = Complex.ofCartesian(0, 0)
     * "(0.0,0.0)"         = Complex.ofCartesian(0, 0)
     * "(-0.0, 0.0)"       = Complex.ofCartesian(-0.0, 0)
     * "(-1.23, 4.56)"     = Complex.ofCartesian(-123, 4.56)
     * "(1e300,-1.1e-2)"   = Complex.ofCartesian(1e300, -1.1e-2)</pre>
     *
     * @param s String representation.
     * @return {@code Complex} number.
     * @throws NullPointerException if the string is null.
     * @throws NumberFormatException if the string does not contain a parsable complex number.
     * @see Double#parseDouble(String)
     * @see #toString()
     */
    public static Complex parse(String s) {
        final int len = s.length();
        if (len < FORMAT_MIN_LEN) {
            throw parsingException("Expected format",
                FORMAT_START + "real" + FORMAT_SEP + "imaginary" + FORMAT_END, null);
        }

        // Confirm start: '('
        if (s.charAt(0) != FORMAT_START) {
            throw parsingException("Expected start", FORMAT_START, null);
        }

        // Confirm end: ')'
        if (s.charAt(len - 1) != FORMAT_END) {
            throw parsingException("Expected end", FORMAT_END, null);
        }

        // Confirm separator ',' is between at least 2 characters from
        // either end: "(x,x)"
        // Count back from the end ignoring the last 2 characters.
        final int sep = s.lastIndexOf(FORMAT_SEP, len - 3);
        if (sep < BEFORE_SEP) {
            throw parsingException("Expected separator between two numbers", FORMAT_SEP, null);
        }

        // Should be no more separators
        if (s.indexOf(FORMAT_SEP, sep + 1) != -1) {
            throw parsingException("Incorrect number of parts, expected only 2 using separator",
                FORMAT_SEP, null);
        }

        // Try to parse the parts

        final String rePart = s.substring(1, sep);
        final double re;
        try {
            re = Double.parseDouble(rePart);
        } catch (final NumberFormatException ex) {
            throw parsingException("Could not parse real part", rePart, ex);
        }

        final String imPart = s.substring(sep + 1, len - 1);
        final double im;
        try {
            im = Double.parseDouble(imPart);
        } catch (final NumberFormatException ex) {
            throw parsingException("Could not parse imaginary part", imPart, ex);
        }

        return ofCartesian(re, im);
    }

    /**
     * Returns {@code true} if either the real <em>or</em> imaginary component of the complex number is NaN
     * <em>and</em> the complex number is not infinite.
     *
     * <p>Note that in contrast to {@link Double#isNaN()}:
     * <ul>
     *   <li>There is more than one complex number that can return {@code true}.
     *   <li>Different representations of NaN can be distinguished by the
     *       {@link #equals(Object) Complex.equals(Object)} method.
     * </ul>
     *
     * @return {@code true} if this instance contains NaN and no infinite parts.
     * @see Double#isNaN(double)
     * @see #isInfinite()
     * @see #equals(Object) Complex.equals(Object)
     */
    public boolean isNaN() {
        if (Double.isNaN(real) || Double.isNaN(imaginary)) {
            return !isInfinite();
        }
        return false;
    }

    /**
     * Returns {@code true} if either real or imaginary component of the complex number is infinite.
     *
     * <p>Note: A complex or imaginary value with at least one infinite part is regarded
     * as an infinity (even if its other part is a NaN).</p>
     *
     * @return {@code true} if this instance contains an infinite value.
     * @see Double#isInfinite(double)
     */
    public boolean isInfinite() {
        return Double.isInfinite(real) || Double.isInfinite(imaginary);
    }

    /**
     * Returns {@code true} if both real and imaginary component of the complex number are finite.
     *
     * @return {@code true} if this instance contains finite values.
     * @see Double#isFinite(double)
     */
    public boolean isFinite() {
        return Double.isFinite(real) && Double.isFinite(imaginary);
    }

    /**
     * Returns projection of this complex number onto the Riemann sphere.
     *
     * <p>\( z \) projects to \( z \), except that all complex infinities (even those
     * with one infinite part and one NaN part) project to positive infinity on the real axis.
     *
     * If \( z \) has an infinite part, then {@code z.proj()} shall be equivalent to:
     *
     * <pre>return Complex.ofCartesian(Double.POSITIVE_INFINITY, Math.copySign(0.0, z.imag());</pre>
     *
     * @return \( z \) projected onto the Riemann sphere.
     * @see #isInfinite()
     * @see <a href="http://pubs.opengroup.org/onlinepubs/9699919799/functions/cproj.html">
     * IEEE and ISO C standards: cproj</a>
     */
    public Complex proj() {
        if (isInfinite()) {
            return new Complex(Double.POSITIVE_INFINITY, Math.copySign(0.0, imaginary));
        }
        return this;
    }

    /**
     * Returns the absolute value of this complex number. This is also called complex norm, modulus,
     * or magnitude.
     *
     * <p>\[ \text{abs}(x + i y) = \sqrt{(x^2 + y^2)} \]
     *
     * <p>If either component is infinite then the result is positive infinity. If either
     * component is NaN and this is not {@link #isInfinite() infinite} then the result is NaN.
     *
     * <p>This code follows the
     * <a href="http://www.iso-9899.info/wiki/The_Standard">ISO C Standard</a>, Annex G,
     * in calculating the returned value using an equivalent function to {@code hypot(x, y)}
     * for complex \( x + i y \). Results are expected to be within 1 ULP of the exact
     * result but due to the implementation may differ from {@link Math#hypot(double, double)}.
     *
     * @return The absolute value.
     * @see #isInfinite()
     * @see #isNaN()
     * @see Math#hypot(double, double)
     * @see <a href="http://mathworld.wolfram.com/ComplexModulus.html">Complex modulus</a>
     */
    public double abs() {
        // Specialised implementation of hypot
        return hypot(real, imaginary);
    }

    /**
     * Returns the squared norm value of this complex number. This is also called the absolute
     * square.
     *
     * <p>\[ \text{norm}(x + i y) = x^2 + y^2 \]
     *
     * <p>If either component is infinite then the result is positive infinity. If either
     * component is NaN and this is not {@link #isInfinite() infinite} then the result is NaN.
     *
     * <p>This method will return the square of {@link #abs()}. It can be used as a faster
     * alternative for ranking by magnitude although overflow to infinity will create equal
     * ranking for values that may be still distinguished by {@code abs()}.
     *
     * @return The square norm value.
     * @see #isInfinite()
     * @see #isNaN()
     * @see #abs()
     * @see <a href="http://mathworld.wolfram.com/AbsoluteSquare.html">Absolute square</a>
     */
    public double norm() {
        if (isInfinite()) {
            return Double.POSITIVE_INFINITY;
        }
        return real * real + imaginary * imaginary;
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this + addend)}.
     * Implements the formula:
     *
     * <p>\[ (a + i b) + (c + i d) = (a + c) + i (b + d) \]
     *
     * @param  addend Value to be added to this complex number.
     * @return {@code this + addend}.
     * @see <a href="http://mathworld.wolfram.com/ComplexAddition.html">Complex Addition</a>
     */
    public Complex add(Complex addend) {
        return new Complex(real + addend.real,
                           imaginary + addend.imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this + addend)},
     * with {@code addend} interpreted as a real number.
     * Implements the formula:
     *
     * <p>\[ (a + i b) + c = (a + c) + i b \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * real-only and complex numbers.</p>
     *
     * <p>Note: This method preserves the sign of the imaginary component \( b \) if it is {@code -0.0}.
     * The sign would be lost if adding \( (c + i 0) \) using
     * {@link #add(Complex) add(Complex.ofCartesian(addend, 0))} since
     * {@code -0.0 + 0.0 = 0.0}.
     *
     * @param addend Value to be added to this complex number.
     * @return {@code this + addend}.
     * @see #add(Complex)
     * @see #ofCartesian(double, double)
     */
    public Complex add(double addend) {
        return new Complex(real + addend, imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this + addend)},
     * with {@code addend} interpreted as an imaginary number.
     * Implements the formula:
     *
     * <p>\[ (a + i b) + i d = a + i (b + d) \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * imaginary-only and complex numbers.</p>
     *
     * <p>Note: This method preserves the sign of the real component \( a \) if it is {@code -0.0}.
     * The sign would be lost if adding \( (0 + i d) \) using
     * {@link #add(Complex) add(Complex.ofCartesian(0, addend))} since
     * {@code -0.0 + 0.0 = 0.0}.
     *
     * @param addend Value to be added to this complex number.
     * @return {@code this + addend}.
     * @see #add(Complex)
     * @see #ofCartesian(double, double)
     */
    public Complex addImaginary(double addend) {
        return new Complex(real, imaginary + addend);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/ComplexConjugate.html">conjugate</a>
     * \( \overline{z} \) of this complex number \( z \).
     *
     * <p>\[ z           = x + i y \\
     *      \overline{z} = x - i y \]
     *
     * @return The conjugate (\( \overline{z} \)) of this complex number.
     */
    public Complex conj() {
        return new Complex(real, -imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this / divisor)}.
     * Implements the formula:
     *
     * <p>\[ \frac{a + i b}{c + i d} = \frac{(ac + bd) + i (bc - ad)}{c^2+d^2} \]
     *
     * <p>Re-calculates NaN result values to recover infinities as specified in C99 standard G.5.1.
     *
     * @param divisor Value by which this complex number is to be divided.
     * @return {@code this / divisor}.
     * @see <a href="http://mathworld.wolfram.com/ComplexDivision.html">Complex Division</a>
     */
    public Complex divide(Complex divisor) {
        return divide(real, imaginary, divisor.real, divisor.imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is:
     * <pre>
     * <code>
     *   a + i b     (ac + bd) + i (bc - ad)
     *   -------  =  -----------------------
     *   c + i d            c<sup>2</sup> + d<sup>2</sup>
     * </code>
     * </pre>
     *
     * <p>Recalculates to recover infinities as specified in C99
     * standard G.5.1. Method is fully in accordance with
     * C++11 standards for complex numbers.</p>
     *
     * <p>Note: In the event of divide by zero this method produces the same result
     * as dividing by a real-only zero using {@link #divide(double)}.
     *
     * @param re1 Real component of first number.
     * @param im1 Imaginary component of first number.
     * @param re2 Real component of second number.
     * @param im2 Imaginary component of second number.
     * @return (a + i b) / (c + i d).
     * @see <a href="http://mathworld.wolfram.com/ComplexDivision.html">Complex Division</a>
     * @see #divide(double)
     */
    private static Complex divide(double re1, double im1, double re2, double im2) {
        double a = re1;
        double b = im1;
        double c = re2;
        double d = im2;
        int ilogbw = 0;
        // Get the exponent to scale the divisor parts to the range [1, 2).
        final int exponent = getScale(c, d);
        if (exponent <= Double.MAX_EXPONENT) {
            ilogbw = exponent;
            c = Math.scalb(c, -ilogbw);
            d = Math.scalb(d, -ilogbw);
        }
        final double denom = c * c + d * d;
        // Divisor is in the range [1, 2).
        // Avoid overflow if a or b are very big. Here we do not need to
        // handle sub-normal exponents so just get the maximum exponent.
        if (getMaxExponent(a, b) > Double.MAX_EXPONENT - 2) {
            ilogbw -= 2;
            a /= 4;
            b /= 4;
        }
        double x = Math.scalb((a * c + b * d) / denom, -ilogbw);
        double y = Math.scalb((b * c - a * d) / denom, -ilogbw);
        // Recover infinities and zeros that computed as NaN+iNaN
        // the only cases are nonzero/zero, infinite/finite, and finite/infinite, ...
        // --------------
        // Modification from the listing in ISO C99 G.5.1 (8):
        // Prevent overflow in (a * c + b * d) and (b * c - a * d).
        // It is only the sign that is important. not the magnitude.
        // --------------
        if (Double.isNaN(x) && Double.isNaN(y)) {
            if ((denom == 0.0) &&
                    (!Double.isNaN(a) || !Double.isNaN(b))) {
                // nonzero/zero
                // This case produces the same result as divide by a real-only zero
                // using divide(+/-0.0).
                x = Math.copySign(Double.POSITIVE_INFINITY, c) * a;
                y = Math.copySign(Double.POSITIVE_INFINITY, c) * b;
            } else if ((Double.isInfinite(a) || Double.isInfinite(b)) &&
                    Double.isFinite(c) && Double.isFinite(d)) {
                // infinite/finite
                a = boxInfinity(a);
                b = boxInfinity(b);
                x = Double.POSITIVE_INFINITY * computeACplusBD(a, b, c, d);
                y = Double.POSITIVE_INFINITY * computeBCminusAD(a, b, c, d);
            } else if ((Double.isInfinite(c) || Double.isInfinite(d)) &&
                    Double.isFinite(a) && Double.isFinite(b)) {
                // finite/infinite
                c = boxInfinity(c);
                d = boxInfinity(d);
                x = 0.0 * computeACplusBD(a, b, c, d);
                y = 0.0 * computeBCminusAD(a, b, c, d);
            }
        }
        return new Complex(x, y);
    }

    /**
     * Compute {@code a*c + b*d} without overflow.
     * It is assumed: either {@code a} an\( b \)b} or {@code c} and {@code d} are
     * either zero or one (i.e. a boxed infinity); and the sign of the result is important,
     * not the value.
     *
     * @param a the a
     * @param b the b
     * @param c the c
     * @param d the d
     * @return The result
     */
    private static double computeACplusBD(double a, double b, double c, double d) {
        final double ac = a * c;
        final double bd = b * d;
        final double result = ac + bd;
        return Double.isFinite(result) ?
            result :
            // Overflow. Just divide by 2 as it is the sign of the result that matters.
            ac * 0.5 + bd * 0.5;
    }

    /**
     * Compute {@code b*c - a*d} without overflow.
     * It is assumed: either {@code a} and {@code b} or {@code c} and {@code d} are
     * either zero or one (i.e. a boxed infinity); and the sign of the result is important,
     * not the value.
     *
     * @param a the a
     * @param b the b
     * @param c the c
     * @param d the d
     * @return The result
     */
    private static double computeBCminusAD(double a, double b, double c, double d) {
        final double bc = b * c;
        final double ad = a * d;
        final double result = bc - ad;
        return Double.isFinite(result) ?
            result :
            // Overflow. Just divide by 2 as it is the sign of the result that matters.
            bc * 0.5 - ad * 0.5;
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this / divisor)},
     * with {@code divisor} interpreted as a real number.
     * Implements the formula:
     *
     * <p>\[ \frac{a + i b}{c} = \frac{a}{c} + i \frac{b}{c} \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * real-only and complex numbers.</p>
     *
     * <p>Note: This method should be preferred over using
     * {@link #divide(Complex) divide(Complex.ofCartesian(divisor, 0))}. Division
     * can generate signed zeros if {@code this} complex has zeros for the real
     * and/or imaginary component, or the divisor is infinity. The summation of signed zeros
     * in {@link #divide(Complex)} may create zeros in the result that differ in sign
     * from the equivalent call to divide by a real-only number.
     *
     * @param  divisor Value by which this complex number is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(Complex)
     */
    public Complex divide(double divisor) {
        return new Complex(real / divisor, imaginary / divisor);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this / divisor)},
     * with {@code divisor} interpreted as an imaginary number.
     * Implements the formula:
     *
     * <p>\[ \frac{a + i b}{id} = \frac{b}{d} - i \frac{a}{d} \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * imaginary-only and complex numbers.</p>
     *
     * <p>Note: This method should be preferred over using
     * {@link #divide(Complex) divide(Complex.ofCartesian(0, divisor))}. Division
     * can generate signed zeros if {@code this} complex has zeros for the real
     * and/or imaginary component, or the divisor is infinity. The summation of signed zeros
     * in {@link #divide(Complex)} may create zeros in the result that differ in sign
     * from the equivalent call to divide by an imaginary-only number.
     *
     * <p>Warning: This method will generate a different result from
     * {@link #divide(Complex) divide(Complex.ofCartesian(0, divisor))} if the divisor is zero.
     * In this case the divide method using a zero-valued Complex will produce the same result
     * as dividing by a real-only zero. The output from dividing by imaginary zero will create
     * infinite and NaN values in the same component parts as the output from
     * {@code this.divide(Complex.ZERO).multiplyImaginary(1)}, however the sign
     * of some infinity values may be negated.
     *
     * @param  divisor Value by which this complex number is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(Complex)
     * @see #divide(double)
     */
    public Complex divideImaginary(double divisor) {
        return new Complex(imaginary / divisor, -real / divisor);
    }

    /**
     * Test for equality with another object. If the other object is a {@code Complex} then a
     * comparison is made of the real and imaginary parts; otherwise {@code false} is returned.
     *
     * <p>If both the real and imaginary parts of two complex numbers
     * are exactly the same the two {@code Complex} objects are considered to be equal.
     * For this purpose, two {@code double} values are considered to be
     * the same if and only if the method {@link Double #doubleToLongBits(double)}
     * returns the identical {@code long} value when applied to each.
     *
     * <p>Note that in most cases, for two instances of class
     * {@code Complex}, {@code c1} and {@code c2}, the
     * value of {@code c1.equals(c2)} is {@code true} if and only if
     *
     * <pre>
     *  {@code c1.getReal() == c2.getReal() && c1.getImaginary() == c2.getImaginary()}</pre>
     *
     * <p>also has the value {@code true}. However, there are exceptions:
     *
     * <ul>
     *  <li>
     *   Instances that contain {@code NaN} values in the same part
     *   are considered to be equal for that part, even though {@code Double.NaN==Double.NaN}
     *   has the value {@code false}.
     *  </li>
     *  <li>
     *   Instances that share a {@code NaN} value in one part
     *   but have different values in the other part are <em>not</em> considered equal.
     *  </li>
     *  <li>
     *   Instances that contain different representations of zero in the same part
     *   are <em>not</em> considered to be equal for that part, even though {@code -0.0==0.0}
     *   has the value {@code true}.
     *  </li>
     * </ul>
     *
     * <p>The behavior is the same as if the components of the two complex numbers were passed
     * to {@link java.util.Arrays#equals(double[], double[]) Arrays.equals(double[], double[])}:
     *
     * <pre>
     *  Arrays.equals(new double[]{c1.getReal(), c1.getImaginary()},
     *                new double[]{c2.getReal(), c2.getImaginary()}); </pre>
     *
     * @param other Object to test for equality with this instance.
     * @return {@code true} if the objects are equal, {@code false} if object
     * is {@code null}, not an instance of {@code Complex}, or not equal to
     * this instance.
     * @see java.lang.Double#doubleToLongBits(double)
     * @see java.util.Arrays#equals(double[], double[])
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof Complex) {
            final Complex c = (Complex) other;
            return equals(real, c.real) &&
                equals(imaginary, c.imaginary);
        }
        return false;
    }

    /**
     * Gets a hash code for the complex number.
     *
     * <p>The behavior is the same as if the components of the complex number were passed
     * to {@link java.util.Arrays#hashCode(double[]) Arrays.hashCode(double[])}:
     *
     * <pre>
     *  {@code Arrays.hashCode(new double[] {getReal(), getImaginary()})}</pre>
     *
     * @return A hash code value for this object.
     * @see java.util.Arrays#hashCode(double[]) Arrays.hashCode(double[])
     */
    @Override
    public int hashCode() {
        return 31 * (31 + Double.hashCode(real)) + Double.hashCode(imaginary);
    }

    /**
     * Gets the imaginary part.
     *
     * @return The imaginary part.
     */
    public double getImaginary() {
        return imaginary;
    }

    /**
     * Gets the imaginary part (C++ grammar).
     *
     * @return The imaginary part.
     * @see #getImaginary()
     */
    public double imag() {
        return getImaginary();
    }

    /**
     * Gets the real part.
     *
     * @return The real part.
     */
    public double getReal() {
        return real;
    }

     /**
     * Gets the real part (C++ grammar).
     *
     * @return The real part.
     * @see #getReal()
     */
    public double real() {
        return getReal();
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}.
     * Implements the formula:
     *
     * <p>\[ (a + i b)(c + i d) = (ac - bd) + i (ad + bc) \]
     *
     * <p>Recalculates to recover infinities as specified in C99 standard G.5.1.
     *
     * @param  factor Value to be multiplied by this complex number.
     * @return {@code this * factor}.
     * @see <a href="http://mathworld.wolfram.com/ComplexMultiplication.html">Complex Muliplication</a>
     */
    public Complex multiply(Complex factor) {
        return multiply(real, imaginary, factor.real, factor.imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is:
     * <pre>
     *  (a + i b)(c + i d) = (ac - bd) + i (ad + bc)</pre>
     *
     * <p>Recalculates to recover infinities as specified in C99 standard G.5.1.
     *
     * @param re1 Real component of first number.
     * @param im1 Imaginary component of first number.
     * @param re2 Real component of second number.
     * @param im2 Imaginary component of second number.
     * @return (a + b i)(c + d i).
     */
    private static Complex multiply(double re1, double im1, double re2, double im2) {
        double a = re1;
        double b = im1;
        double c = re2;
        double d = im2;
        final double ac = a * c;
        final double bd = b * d;
        final double ad = a * d;
        final double bc = b * c;
        double x = ac - bd;
        double y = ad + bc;

        // --------------
        // NaN can occur if:
        // - any of (a,b,c,d) are NaN (for NaN or Infinite complex numbers)
        // - a multiplication of infinity by zero (ac,bd,ad,bc).
        // - a subtraction of infinity from infinity (e.g. ac - bd)
        //   Note that (ac,bd,ad,bc) can be infinite due to overflow.
        //
        // Detect a NaN result and perform correction.
        //
        // Modification from the listing in ISO C99 G.5.1 (6)
        // Do not correct infinity multiplied by zero. This is left as NaN.
        // --------------

        if (Double.isNaN(x) && Double.isNaN(y)) {
            // Recover infinities that computed as NaN+iNaN ...
            boolean recalc = false;
            if ((Double.isInfinite(a) || Double.isInfinite(b)) &&
                isNotZero(c, d)) {
                // This complex is infinite.
                // "Box" the infinity and change NaNs in the other factor to 0.
                a = boxInfinity(a);
                b = boxInfinity(b);
                c = changeNaNtoZero(c);
                d = changeNaNtoZero(d);
                recalc = true;
            }
            // (c, d) may have been corrected so do not use factor.isInfinite().
            if ((Double.isInfinite(c) || Double.isInfinite(d)) &&
                isNotZero(a, b)) {
                // This other complex is infinite.
                // "Box" the infinity and change NaNs in the other factor to 0.
                c = boxInfinity(c);
                d = boxInfinity(d);
                a = changeNaNtoZero(a);
                b = changeNaNtoZero(b);
                recalc = true;
            }
            if (!recalc && (Double.isInfinite(ac) || Double.isInfinite(bd) ||
                            Double.isInfinite(ad) || Double.isInfinite(bc))) {
                // The result overflowed to infinity.
                // Recover infinities from overflow by changing NaNs to 0 ...
                a = changeNaNtoZero(a);
                b = changeNaNtoZero(b);
                c = changeNaNtoZero(c);
                d = changeNaNtoZero(d);
                recalc = true;
            }
            if (recalc) {
                x = Double.POSITIVE_INFINITY * (a * c - b * d);
                y = Double.POSITIVE_INFINITY * (a * d + b * c);
            }
        }
        return new Complex(x, y);
    }

    /**
     * Box values for the real or imaginary component of an infinite complex number.
     * Any infinite value will be returned as one. Non-infinite values will be returned as zero.
     * The sign is maintained.
     *
     * <pre>
     *  inf  =  1
     * -inf  = -1
     *  x    =  0
     * -x    = -0
     * </pre>
     *
     * @param component the component
     * @return The boxed value
     */
    private static double boxInfinity(double component) {
        return Math.copySign(Double.isInfinite(component) ? 1.0 : 0.0, component);
    }

    /**
     * Checks if the complex number is not zero.
     *
     * @param real the real component
     * @param imaginary the imaginary component
     * @return true if the complex is not zero
     */
    private static boolean isNotZero(double real, double imaginary) {
        // The use of equals is deliberate.
        // This method must distinguish NaN from zero thus ruling out:
        // (real != 0.0 || imaginary != 0.0)
        return !(real == 0.0 && imaginary == 0.0);
    }

    /**
     * Change NaN to zero preserving the sign; otherwise return the value.
     *
     * @param value the value
     * @return The new value
     */
    private static double changeNaNtoZero(double value) {
        return Double.isNaN(value) ? Math.copySign(0.0, value) : value;
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a real number.
     * Implements the formula:
     *
     * <p>\[ (a + i b) c =  (ac) + i (bc) \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * real-only and complex numbers.</p>
     *
     * <p>Note: This method should be preferred over using
     * {@link #multiply(Complex) multiply(Complex.ofCartesian(factor, 0))}. Multiplication
     * can generate signed zeros if either {@code this} complex has zeros for the real
     * and/or imaginary component, or if the factor is zero. The summation of signed zeros
     * in {@link #multiply(Complex)} may create zeros in the result that differ in sign
     * from the equivalent call to multiply by a real-only number.
     *
     * @param  factor Value to be multiplied by this complex number.
     * @return {@code this * factor}.
     * @see #multiply(Complex)
     */
    public Complex multiply(double factor) {
        return new Complex(real * factor, imaginary * factor);
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as an imaginary number.
     * Implements the formula:
     *
     * <p>\[ (a + i b) id = (-bd) + i (ad) \]
     *
     * <p>This method can be used to compute the multiplication of this complex number \( z \)
     * by \( i \). This should be used in preference to
     * {@link #multiply(Complex) multiply(Complex.I)} with or without {@link #negate() negation}:</p>
     *
     * \[ iz = (-b + i a) \\
     *   -iz = (b - i a) \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * imaginary-only and complex numbers.</p>
     *
     * <p>Note: This method should be preferred over using
     * {@link #multiply(Complex) multiply(Complex.ofCartesian(0, factor))}. Multiplication
     * can generate signed zeros if either {@code this} complex has zeros for the real
     * and/or imaginary component, or if the factor is zero. The summation of signed zeros
     * in {@link #multiply(Complex)} may create zeros in the result that differ in sign
     * from the equivalent call to multiply by an imaginary-only number.
     *
     * @param  factor Value to be multiplied by this complex number.
     * @return {@code this * factor}.
     * @see #multiply(Complex)
     */
    public Complex multiplyImaginary(double factor) {
        return new Complex(-imaginary * factor, real * factor);
    }

    /**
     * Returns a {@code Complex} whose value is the negation of both the real and imaginary parts
     * of complex number \( z \).
     *
     * @return \( -z \).
     */
    public Complex negate() {
        return new Complex(-real, -imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this - subtrahend)}.
     * Implements the formula:
     *
     * <p>\[ (a + i b) - (c + i d) = (a - c) + i (b - d) \]
     *
     * @param  subtrahend Value to be subtracted from this complex number.
     * @return {@code this - subtrahend}.
     * @see <a href="http://mathworld.wolfram.com/ComplexSubtraction.html">Complex Subtraction</a>
     */
    public Complex subtract(Complex subtrahend) {
        return new Complex(real - subtrahend.real,
                           imaginary - subtrahend.imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this - subtrahend)},
     * with {@code subtrahend} interpreted as a real number.
     * Implements the formula:
     *
     * <p>\[ (a + i b) - c = (a - c) + i b \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * real-only and complex numbers.</p>
     *
     * @param  subtrahend Value to be subtracted from this complex number.
     * @return {@code this - subtrahend}.
     * @see #subtract(Complex)
     */
    public Complex subtract(double subtrahend) {
        return new Complex(real - subtrahend, imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this - subtrahend)},
     * with {@code subtrahend} interpreted as an imaginary number.
     * Implements the formula:
     *
     * <p>\[ (a + i b) - i d = a + i (b - d) \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * imaginary-only and complex numbers.</p>
     *
     * @param  subtrahend Value to be subtracted from this complex number.
     * @return {@code this - subtrahend}.
     * @see #subtract(Complex)
     */
    public Complex subtractImaginary(double subtrahend) {
        return new Complex(real, imaginary - subtrahend);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (minuend - this)},
     * with {@code minuend} interpreted as a real number.
     * Implements the formula:
     * \[ c - (a + i b) = (c - a) - i b \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * real-only and complex numbers.</p>
     *
     * <p>Note: This method inverts the sign of the imaginary component \( b \) if it is {@code 0.0}.
     * The sign would not be inverted if subtracting from \( c + i 0 \) using
     * {@link #subtract(Complex) Complex.ofCartesian(minuend, 0).subtract(this))} since
     * {@code 0.0 - 0.0 = 0.0}.
     *
     * @param  minuend Value this complex number is to be subtracted from.
     * @return {@code minuend - this}.
     * @see #subtract(Complex)
     * @see #ofCartesian(double, double)
     */
    public Complex subtractFrom(double minuend) {
        return new Complex(minuend - real, -imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this - subtrahend)},
     * with {@code minuend} interpreted as an imaginary number.
     * Implements the formula:
     * \[ i d - (a + i b) = -a + i (d - b) \]
     *
     * <p>This method is included for compatibility with ISO C99 which defines arithmetic between
     * imaginary-only and complex numbers.</p>
     *
     * <p>Note: This method inverts the sign of the real component \( a \) if it is {@code 0.0}.
     * The sign would not be inverted if subtracting from \( 0 + i d \) using
     * {@link #subtract(Complex) Complex.ofCartesian(0, minuend).subtract(this))} since
     * {@code 0.0 - 0.0 = 0.0}.
     *
     * @param  minuend Value this complex number is to be subtracted from.
     * @return {@code this - subtrahend}.
     * @see #subtract(Complex)
     * @see #ofCartesian(double, double)
     */
    public Complex subtractFromImaginary(double minuend) {
        return new Complex(-real, minuend - imaginary);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/InverseCosine.html">
     * inverse cosine</a> of this complex number.
     *
     * <p>\[ \cos^{-1}(z) = \frac{\pi}{2} + i \left(\ln{iz + \sqrt{1 - z^2}}\right) \]
     *
     * <p>The inverse cosine of \( z \) is in the range \( [0, \pi) \) along the real axis and
     * unbounded along the imaginary axis. Special cases:
     *
     * <ul>
     * <li>{@code z.conj().acos() == z.acos().conj()}.
     * <li>If {@code z} is ±0 + i0, returns π/2 − i0.
     * <li>If {@code z} is ±0 + iNaN, returns π/2 + iNaN.
     * <li>If {@code z} is x + i∞ for finite x, returns π/2 − i∞.
     * <li>If {@code z} is x + iNaN, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is −∞ + iy for positive-signed finite y, returns π − i∞.
     * <li>If {@code z} is +∞ + iy for positive-signed finite y, returns +0 − i∞.
     * <li>If {@code z} is −∞ + i∞, returns 3π/4 − i∞.
     * <li>If {@code z} is +∞ + i∞, returns π/4 − i∞.
     * <li>If {@code z} is ±∞ + iNaN, returns NaN ± i∞ where the sign of the imaginary part of the result is unspecified.
     * <li>If {@code z} is NaN + iy for finite y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + i∞, returns NaN − i∞.
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>The inverse cosine is a multivalued function and requires a branch cut in
     * the complex plane; the cut is conventionally placed at the line segments
     * \( (-\infty,-1) \) and \( (1,\infty) \) of the real axis.
     *
     * <p>This function is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \cos^{-1}(z) = \cos^{-1}(B) - i\ \text{sgn}(y) \ln\left(A + \sqrt{A^2-1}\right) \\
     *   A = \frac{1}{2} \left[ \sqrt{(x+1)^2+y^2} + \sqrt{(x-1)^2+y^2} \right] \\
     *   B = \frac{1}{2} \left[ \sqrt{(x+1)^2+y^2} - \sqrt{(x-1)^2+y^2} \right] \]
     *
     * <p>where \( \text{sgn}(y) \) is the sign function implemented using
     * {@link Math#copySign(double,double) copySign(1.0, y)}.
     *
     * <p>The implementation is based on the method described in:</p>
     * <blockquote>
     * T E Hull, Thomas F Fairgrieve and Ping Tak Peter Tang (1997)
     * Implementing the complex Arcsine and Arccosine Functions using Exception Handling.
     * ACM Transactions on Mathematical Software, Vol 23, No 3, pp 299-335.
     * </blockquote>
     *
     * <p>The code has been adapted from the <a href="https://www.boost.org/">Boost</a>
     * {@code c++} implementation {@code <boost/math/complex/acos.hpp>}. The function is well
     * defined over the entire complex number range, and produces accurate values even at the
     * extremes due to special handling of overflow and underflow conditions.</p>
     *
     * @return The inverse cosine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/ArcCos/">ArcCos</a>
     */
    public Complex acos() {
        return acos(real, imaginary, Complex::ofCartesian);
    }

    /**
     * Returns the inverse cosine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code acosh(z) = +-i acos(z)}.<p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param constructor Constructor.
     * @return The inverse cosine of the complex number.
     */
    private static Complex acos(final double real, final double imaginary,
                                final ComplexConstructor constructor) {
        // Compute with positive values and determine sign at the end
        final double x = Math.abs(real);
        final double y = Math.abs(imaginary);
        // The result (without sign correction)
        double re;
        double im;

        // Handle C99 special cases
        if (isPosInfinite(x)) {
            if (isPosInfinite(y)) {
                re = PI_OVER_4;
                im = y;
            } else if (Double.isNaN(y)) {
                // sign of the imaginary part of the result is unspecified
                return constructor.create(imaginary, real);
            } else {
                re = 0;
                im = Double.POSITIVE_INFINITY;
            }
        } else if (Double.isNaN(x)) {
            if (isPosInfinite(y)) {
                return constructor.create(x, -imaginary);
            }
            // No-use of the input constructor
            return NAN;
        } else if (isPosInfinite(y)) {
            re = PI_OVER_2;
            im = y;
        } else if (Double.isNaN(y)) {
            return constructor.create(x == 0 ? PI_OVER_2 : y, y);
        } else {
            // Special case for real numbers:
            if (y == 0 && x <= 1) {
                return constructor.create(x == 0 ? PI_OVER_2 : Math.acos(real), -imaginary);
            }

            final double xp1 = x + 1;
            final double xm1 = x - 1;

            if (inRegion(x, y, SAFE_MIN, SAFE_MAX)) {
                final double yy = y * y;
                final double r = Math.sqrt(xp1 * xp1 + yy);
                final double s = Math.sqrt(xm1 * xm1 + yy);
                final double a = 0.5 * (r + s);
                final double b = x / a;

                if (b <= B_CROSSOVER) {
                    re = Math.acos(b);
                } else {
                    final double apx = a + x;
                    if (x <= 1) {
                        re = Math.atan(Math.sqrt(0.5 * apx * (yy / (r + xp1) + (s - xm1))) / x);
                    } else {
                        re = Math.atan((y * Math.sqrt(0.5 * (apx / (r + xp1) + apx / (s + xm1)))) / x);
                    }
                }

                if (a <= A_CROSSOVER) {
                    double am1;
                    if (x < 1) {
                        am1 = 0.5 * (yy / (r + xp1) + yy / (s - xm1));
                    } else {
                        am1 = 0.5 * (yy / (r + xp1) + (s + xm1));
                    }
                    im = Math.log1p(am1 + Math.sqrt(am1 * (a + 1)));
                } else {
                    im = Math.log(a + Math.sqrt(a * a - 1));
                }
            } else {
                // Hull et al: Exception handling code from figure 6
                if (y <= (EPSILON * Math.abs(xm1))) {
                    if (x < 1) {
                        re = Math.acos(x);
                        im = y / Math.sqrt(xp1 * (1 - x));
                    } else {
                        // This deviates from Hull et al's paper as per
                        // https://svn.boost.org/trac/boost/ticket/7290
                        if ((Double.MAX_VALUE / xp1) > xm1) {
                            // xp1 * xm1 won't overflow:
                            re = y / Math.sqrt(xm1 * xp1);
                            im = Math.log1p(xm1 + Math.sqrt(xp1 * xm1));
                        } else {
                            re = y / x;
                            im = LN_2 + Math.log(x);
                        }
                    }
                } else if (y <= SAFE_MIN) {
                    // Hull et al: Assume x == 1.
                    // True if:
                    // E^2 > 8*sqrt(u)
                    //
                    // E = Machine epsilon: (1 + epsilon) = 1
                    // u = Double.MIN_NORMAL
                    re = Math.sqrt(y);
                    im = Math.sqrt(y);
                } else if (EPSILON * y - 1 >= x) {
                    re = PI_OVER_2;
                    im = LN_2 + Math.log(y);
                } else if (x > 1) {
                    re = Math.atan(y / x);
                    final double xoy = x / y;
                    im = LN_2 + Math.log(y) + 0.5 * Math.log1p(xoy * xoy);
                } else {
                    re = PI_OVER_2;
                    final double a = Math.sqrt(1 + y * y);
                    im = 0.5 * Math.log1p(2 * y * (y + a));
                }
            }
        }

        return constructor.create(negative(real) ? Math.PI - re : re,
                                  negative(imaginary) ? im : -im);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/InverseSine.html">
     * inverse sine</a> of this complex number.
     *
     * <p>\[ \sin^{-1}(z) = - i \left(\ln{iz + \sqrt{1 - z^2}}\right) \]
     *
     * <p>The inverse sine of \( z \) is unbounded along the imaginary axis and
     * in the range \( [-\pi, \pi] \) along the real axis. Special cases are handled
     * as if the operation is implemented using \( \sin^{-1}(z) = -i \sinh^{-1}(iz) \).
     *
     * <p>The inverse sine is a multivalued function and requires a branch cut in
     * the complex plane; the cut is conventionally placed at the line segments
     * \( (\infty,-1) \) and \( (1,\infty) \) of the real axis.
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \sin^{-1}(z) = \sin^{-1}(B) + i\ \text{sgn}(y)\ln \left(A + \sqrt{A^2-1} \right) \\
     *   A = \frac{1}{2} \left[ \sqrt{(x+1)^2+y^2} + \sqrt{(x-1)^2+y^2} \right] \\
     *   B = \frac{1}{2} \left[ \sqrt{(x+1)^2+y^2} - \sqrt{(x-1)^2+y^2} \right] \]
     *
     * <p>where \( \text{sgn}(y) \) is the sign function implemented using
     * {@link Math#copySign(double,double) copySign(1.0, y)}.
     *
     * <p>The implementation is based on the method described in:</p>
     * <blockquote>
     * T E Hull, Thomas F Fairgrieve and Ping Tak Peter Tang (1997)
     * Implementing the complex Arcsine and Arccosine Functions using Exception Handling.
     * ACM Transactions on Mathematical Software, Vol 23, No 3, pp 299-335.
     * </blockquote>
     *
     * <p>The code has been adapted from the <a href="https://www.boost.org/">Boost</a>
     * {@code c++} implementation {@code <boost/math/complex/asin.hpp>}. The function is well
     * defined over the entire complex number range, and produces accurate values even at the
     * extremes due to special handling of overflow and underflow conditions.</p>
     *
     * @return The inverse sine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/ArcSin/">ArcSin</a>
     */
    public Complex asin() {
        return asin(real, imaginary, Complex::ofCartesian);
    }

    /**
     * Returns the inverse sine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code asinh(z) = -i asin(iz)}.<p>
     *
     * <p>The code has been adapted from the <a href="https://www.boost.org/">Boost</a>
     * {@code c++} implementation {@code <boost/math/complex/asin.hpp>}.</p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param constructor Constructor.
     * @return The inverse sine of this complex number.
     */
    private static Complex asin(final double real, final double imaginary,
                                final ComplexConstructor constructor) {
        // Compute with positive values and determine sign at the end
        final double x = Math.abs(real);
        final double y = Math.abs(imaginary);
        // The result (without sign correction)
        double re;
        double im;

        // Handle C99 special cases
        if (Double.isNaN(x)) {
            if (isPosInfinite(y)) {
                re = x;
                im = y;
            } else {
                // No-use of the input constructor
                return NAN;
            }
        } else if (Double.isNaN(y)) {
            if (x == 0) {
                re = 0;
                im = y;
            } else if (isPosInfinite(x)) {
                re = y;
                im = x;
            } else {
                // No-use of the input constructor
                return NAN;
            }
        } else if (isPosInfinite(x)) {
            re = isPosInfinite(y) ? PI_OVER_4 : PI_OVER_2;
            im = x;
        } else if (isPosInfinite(y)) {
            re = 0;
            im = y;
        } else {
            // Special case for real numbers:
            if (y == 0 && x <= 1) {
                return constructor.create(Math.asin(real), imaginary);
            }

            final double xp1 = x + 1;
            final double xm1 = x - 1;

            if (inRegion(x, y, SAFE_MIN, SAFE_MAX)) {
                final double yy = y * y;
                final double r = Math.sqrt(xp1 * xp1 + yy);
                final double s = Math.sqrt(xm1 * xm1 + yy);
                final double a = 0.5 * (r + s);
                final double b = x / a;

                if (b <= B_CROSSOVER) {
                    re = Math.asin(b);
                } else {
                    final double apx = a + x;
                    if (x <= 1) {
                        re = Math.atan(x / Math.sqrt(0.5 * apx * (yy / (r + xp1) + (s - xm1))));
                    } else {
                        re = Math.atan(x / (y * Math.sqrt(0.5 * (apx / (r + xp1) + apx / (s + xm1)))));
                    }
                }

                if (a <= A_CROSSOVER) {
                    double am1;
                    if (x < 1) {
                        am1 = 0.5 * (yy / (r + xp1) + yy / (s - xm1));
                    } else {
                        am1 = 0.5 * (yy / (r + xp1) + (s + xm1));
                    }
                    im = Math.log1p(am1 + Math.sqrt(am1 * (a + 1)));
                } else {
                    im = Math.log(a + Math.sqrt(a * a - 1));
                }
            } else {
                // Hull et al: Exception handling code from figure 4
                if (y <= (EPSILON * Math.abs(xm1))) {
                    if (x < 1) {
                        re = Math.asin(x);
                        im = y / Math.sqrt(xp1 * (1 - x));
                    } else {
                        re = PI_OVER_2;
                        if ((Double.MAX_VALUE / xp1) > xm1) {
                            // xp1 * xm1 won't overflow:
                            im = Math.log1p(xm1 + Math.sqrt(xp1 * xm1));
                        } else {
                            im = LN_2 + Math.log(x);
                        }
                    }
                } else if (y <= SAFE_MIN) {
                    // Hull et al: Assume x == 1.
                    // True if:
                    // E^2 > 8*sqrt(u)
                    //
                    // E = Machine epsilon: (1 + epsilon) = 1
                    // u = Double.MIN_NORMAL
                    re = PI_OVER_2 - Math.sqrt(y);
                    im = Math.sqrt(y);
                } else if (EPSILON * y - 1 >= x) {
                    // Possible underflow:
                    re = x / y;
                    im = LN_2 + Math.log(y);
                } else if (x > 1) {
                    re = Math.atan(x / y);
                    final double xoy = x / y;
                    im = LN_2 + Math.log(y) + 0.5 * Math.log1p(xoy * xoy);
                } else {
                    final double a = Math.sqrt(1 + y * y);
                    // Possible underflow:
                    re = x / a;
                    im = 0.5 * Math.log1p(2 * y * (y + a));
                }
            }
        }

        return constructor.create(changeSign(re, real),
                                  changeSign(im, imaginary));
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/InverseTangent.html">
     * inverse tangent</a> of this complex number.
     *
     * <p>\[ \tan^{-1}(z) = \frac{i}{2} \ln \left( \frac{i + z}{i - z} \right) \]
     *
     * <p>The inverse hyperbolic tangent of \( z \) is unbounded along the imaginary axis and
     * in the range \( [-\pi/2, \pi/2] \) along the real axis.
     *
     * <p>The inverse tangent is a multivalued function and requires a branch cut in
     * the complex plane; the cut is conventionally placed at the line segments
     * \( (i \infty,-i] \) and \( [i,i \infty) \) of the imaginary axis.
     *
     * <p>As per the C99 standard this function is computed using the trigonomic identity:
     * \[ \tan^{-1}(z) = -i \tanh^{-1}(iz) \]
     *
     * @return The inverse tangent of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/ArcTan/">ArcTan</a>
     */
    public Complex atan() {
        // Define in terms of atanh
        // atan(z) = -i atanh(iz)
        // Multiply this number by I, compute atanh, then multiply by back
        return atanh(-imaginary, real, Complex::multiplyNegativeI);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/InverseHyperbolicSine.html">
     * inverse hyperbolic sine</a> of this complex number.
     *
     * <p>\[ \sinh^{-1}(z) = \ln \left(z + \sqrt{1 + z^2} \right) \]
     *
     * <p>The inverse hyperbolic sine of \( z \) is unbounded along the real axis and
     * in the range \( [-\pi, \pi] \) along the imaginary axis. Special cases:
     *
     * <ul>
     * <li>{@code z.conj().asinh() == z.asinh().conj()}.
     * <li>This is an odd function: \( \sinh^{-1}(z) = -\sinh^{-1}(-z) \).
     * <li>If {@code z} is +0 + i0, returns 0 + i0.
     * <li>If {@code z} is x + i∞ for positive-signed finite x, returns +∞ + iπ/2.
     * <li>If {@code z} is x + iNaN for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + iy for positive-signed finite y, returns +∞ + i0.
     * <li>If {@code z} is +∞ + i∞, returns +∞ + iπ/4.
     * <li>If {@code z} is +∞ + iNaN, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + i0, returns NaN + i0.
     * <li>If {@code z} is NaN + iy for finite nonzero y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + i∞, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>The inverse hyperbolic sine is a multivalued function and requires a branch cut in
     * the complex plane; the cut is conventionally placed at the line segments
     * \( (-i \infty,-i) \) and \( (i,i \infty) \) of the imaginary axis.
     *
     * <p>This function is computed using the trigonomic identity:
     *
     * <p>\[ \sinh^{-1}(z) = -i \sin^{-1}(iz) \]
     *
     * @return The inverse hyperbolic sine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/ArcSinh/">ArcSinh</a>
     */
    public Complex asinh() {
        // Define in terms of asin
        // asinh(z) = -i asin(iz)
        // Note: This is the opposite to the identity defined in the C99 standard:
        // asin(z) = -i asinh(iz)
        // Multiply this number by I, compute asin, then multiply by back
        return asin(-imaginary, real, Complex::multiplyNegativeI);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/InverseHyperbolicTangent.html">
     * inverse hyperbolic tangent</a> of this complex number.
     *
     * <p>\[ \tanh^{-1}(z) = \frac{1}{2} \ln \left( \frac{1 + z}{1 - z} \right) \]
     *
     * <p>The inverse hyperbolic tangent of \( z \) is unbounded along the real axis and
     * in the range \( [-\pi/2, \pi/2] \) along the imaginary axis. Special cases:
     *
     * <ul>
     * <li>{@code z.conj().atanh() == z.atanh().conj()}.
     * <li>This is an odd function: \( \tanh^{-1}(z) = -\tanh^{-1}(-z) \).
     * <li>If {@code z} is +0 + i0, returns +0 + i0.
     * <li>If {@code z} is +0 + iNaN, returns +0 + iNaN.
     * <li>If {@code z} is +1 + i0, returns +∞ + i0 ("divide-by-zero" floating-point operation).
     * <li>If {@code z} is x + i∞ for finite positive-signed x, returns +0 + iπ/2.
     * <li>If {@code z} is x+iNaN for nonzero finite x, returns NaN+iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + iy for finite positive-signed y, returns +0 + iπ/2.
     * <li>If {@code z} is +∞ + i∞, returns +0 + iπ/2.
     * <li>If {@code z} is +∞ + iNaN, returns +0 + iNaN.
     * <li>If {@code z} is NaN+iy for finite y, returns NaN+iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + i∞, returns ±0 + iπ/2 (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>The inverse hyperbolic tangent is a multivalued function and requires a branch cut in
     * the complex plane; the cut is conventionally placed at the line segments
     * \( (\infty,-1] \) and \( [1,\infty) \) of the real axis.
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \tanh^{-1}(z) = \frac{1}{4} \ln \left(1 + \frac{4x}{(1-x)^2+y^2} \right) + \\
     *                     i \frac{1}{2} \left( \tan^{-1} \left(\frac{2y}{1-x^2-y^2} \right) + \frac{\pi}{2} \left(\text{sgn}(x^2+y^2-1)+1 \right) \text{sgn}(y) \right) \]
     *
     * <p>The imaginary part is computed using {@link Math#atan2(double, double)} to ensure the
     * correct quadrant is returned from \( \tan^{-1} \left(\frac{2y}{1-x^2-y^2} \right) \).
     *
     * <p>The code has been adapted from the <a href="https://www.boost.org/">Boost</a>
     * {@code c++} implementation {@code <boost/math/complex/atanh.hpp>}. The function is well
     * defined over the entire complex number range, and produces accurate values even at the
     * extremes due to special handling of overflow and underflow conditions.
     *
     * @return The inverse hyperbolic tangent of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/ArcTanh/">ArcTanh</a>
     */
    public Complex atanh() {
        return atanh(real, imaginary, Complex::ofCartesian);
    }

    /**
     * Returns the inverse hyperbolic tangent of this complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code atan(z) = -i atanh(iz)}.<p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param constructor Constructor.
     * @return The inverse hyperbolic tangent of the complex number.
     */
    private static Complex atanh(final double real, final double imaginary,
                                 final ComplexConstructor constructor) {
        // Compute with positive values and determine sign at the end
        double x = Math.abs(real);
        double y = Math.abs(imaginary);
        // The result (without sign correction)
        double re;
        double im;

        // Handle C99 special cases
        if (Double.isNaN(x)) {
            if (isPosInfinite(y)) {
                // The sign of the real part of the result is unspecified
                return constructor.create(0, Math.copySign(PI_OVER_2, imaginary));
            }
            // No-use of the input constructor.
            // Optionally raises the ‘‘invalid’’ floating-point exception, for finite y.
            return NAN;
        } else if (Double.isNaN(y)) {
            if (isPosInfinite(x)) {
                return constructor.create(Math.copySign(0, real), Double.NaN);
            }
            if (x == 0) {
                return constructor.create(real, Double.NaN);
            }
            // No-use of the input constructor
            return NAN;
        } else {
            // x && y are finite or infinite.

            // Check the safe region.
            // The lower and upper bounds have been copied from boost::math::atanh.
            // They are different from the safe region for asin and acos.
            // x >= SAFE_UPPER: (1-x) == -x
            // x <= SAFE_LOWER: 1 - x^2 = 1

            if (inRegion(x, y, SAFE_LOWER, SAFE_UPPER)) {
                // Normal computation within a safe region.

                // minus x plus 1: (-x+1)
                final double mxp1 = 1 - x;
                final double yy = y * y;
                // The definition of real component is:
                // real = log( ((x+1)^2+y^2) / ((1-x)^2+y^2) ) / 4
                // This simplifies by adding 1 and subtracting 1 as a fraction:
                //      = log(1 + ((x+1)^2+y^2) / ((1-x)^2+y^2) - ((1-x)^2+y^2)/((1-x)^2+y^2) ) / 4
                //
                // real(atanh(z)) == log(1 + 4*x / ((1-x)^2+y^2)) / 4
                // imag(atanh(z)) == tan^-1 (2y, (1-x)(1+x) - y^2) / 2
                // imag(atanh(z)) == tan^-1 (2y, (1 - x^2 - y^2) / 2
                // The division is done at the end of the function.
                re = Math.log1p(4 * x / (mxp1 * mxp1 + yy));
                // Modified from boost which does not switch the magnitude of x and y.
                // The denominator for atan2 is 1 - x^2 - y^2.
                // This can be made more precise if |x| > |y|.
                final double numerator = 2 * y;
                double denominator;
                if (x < y) {
                    final double tmp = x;
                    x = y;
                    y = tmp;
                }
                // 1 - x is precise if |x| >= 1
                if (x >= 1) {
                    denominator = (1 - x) * (1 + x) - y * y;
                } else {
                    // |x| < 1: Use high precision if possible:
                    // 1 - x^2 - y^2 = -(x^2 + y^2 - 1)
                    denominator = -x2y2m1(x, y);
                }
                im = Math.atan2(numerator, denominator);
            } else {
                // This section handles exception cases that would normally cause
                // underflow or overflow in the main formulas.

                // C99. G.7: Special case for imaginary only numbers
                if (x == 0) {
                    if (imaginary == 0) {
                        return constructor.create(real, imaginary);
                    }
                    // atanh(iy) = i atan(y)
                    return constructor.create(real, Math.atan(imaginary));
                }

                // Real part:
                // real = Math.log1p(4x / ((1-x)^2 + y^2))
                // real = Math.log1p(4x / (1 - 2x + x^2 + y^2))
                // real = Math.log1p(4x / (1 + x(x-2) + y^2))
                // without either overflow or underflow in the squared terms.
                if (x >= SAFE_UPPER) {
                    // (1-x) = -x to machine precision:
                    // log1p(4x / (x^2 + y^2))
                    if (isPosInfinite(x) || isPosInfinite(y)) {
                        re = 0;
                    } else if (y >= SAFE_UPPER) {
                        // Big x and y: divide by x*y
                        re = Math.log1p((4 / y) / (x / y + y / x));
                    } else if (y > 1) {
                        // Big x: divide through by x:
                        re = Math.log1p(4 / (x + y * y / x));
                    } else {
                        // Big x small y, as above but neglect y^2/x:
                        re = Math.log1p(4 / x);
                    }
                } else if (y >= SAFE_UPPER) {
                    if (x > 1) {
                        // Big y, medium x, divide through by y:
                        final double mxp1 = 1 - x;
                        re = Math.log1p((4 * x / y) / (mxp1 * mxp1 / y + y));
                    } else {
                        // Big y, small x, as above but neglect (1-x)^2/y:
                        // Note: log1p(v) == v - v^2/2 + v^3/3 ... Taylor series when v is small.
                        // Here v is so small only the first term matters.
                        re = 4 * x / y / y;
                    }
                } else if (x == 1) {
                    // x = 1, small y:
                    // Special case when x == 1 as (1-x) is invalid.
                    // Simplify the following formula:
                    // real = log( sqrt((x+1)^2+y^2) ) / 2 - log( sqrt((1-x)^2+y^2) ) / 2
                    //      = log( sqrt(4+y^2) ) / 2 - log(y) / 2
                    // if: 4+y^2 -> 4
                    //      = log( 2 ) / 2 - log(y) / 2
                    //      = (log(2) - log(y)) / 2
                    // Multiply by 2 as it will be divided by 4 at the end.
                    // C99: if y=0 raises the ‘‘divide-by-zero’’ floating-point exception.
                    re = 2 * (LN_2 - Math.log(y));
                } else {
                    // Modified from boost which checks y > SAFE_LOWER.
                    // if y*y -> 0 it will be ignored so always include it.
                    final double mxp1 = 1 - x;
                    re = Math.log1p((4 * x) / (mxp1 * mxp1 + y * y));
                }

                // Imaginary part:
                // imag = atan2(2y, (1-x)(1+x) - y^2)
                // if x or y are large, then the formula:
                //   atan2(2y, (1-x)(1+x) - y^2)
                // evaluates to +(PI - theta) where theta is negligible compared to PI.
                if ((x >= SAFE_UPPER) || (y >= SAFE_UPPER)) {
                    im = Math.PI;
                } else if (x <= SAFE_LOWER) {
                    // (1-x)^2 -> 1
                    if (y <= SAFE_LOWER) {
                        // 1 - y^2 -> 1
                        im = Math.atan2(2 * y, 1);
                    } else {
                        im = Math.atan2(2 * y, 1 - y * y);
                    }
                } else {
                    // Medium x, small y.
                    // Modified from boost which checks (y == 0) && (x == 1) and sets re = 0.
                    // This is same as the result from calling atan2(0, 0) so exclude this case.
                    // 1 - y^2 = 1 so ignore subtracting y^2
                    im = Math.atan2(2 * y, (1 - x) * (1 + x));
                }
            }
        }

        re /= 4;
        im /= 2;
        return constructor.create(changeSign(re, real),
                                  changeSign(im, imaginary));
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/InverseHyperbolicCosine.html">
     * inverse hyperbolic cosine</a> of this complex number.
     *
     * <p>\[ \cosh^{-1}(z) = \ln \left(z + \sqrt{z + 1} \sqrt{z - 1} \right) \]
     *
     * <p>The inverse hyperbolic cosine of \( z \) is in the range \( [0, \infty) \) along the
     * real axis and in the range \( [-\pi, \pi] \) along the imaginary axis. Special cases:
     *
     * <ul>
     * <li>{@code z.conj().acosh() == z.acosh().conj()}.
     * <li>If {@code z} is ±0 + i0, returns +0 + iπ/2.
     * <li>If {@code z} is x + i∞ for finite x, returns +∞ + iπ/2.
     * <li>If {@code z} is 0 + iNaN, returns NaN + iπ/2 <sup>[1]</sup>.
     * <li>If {@code z} is x + iNaN for finite non-zero x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is −∞ + iy for positive-signed finite y, returns +∞ + iπ.
     * <li>If {@code z} is +∞ + iy for positive-signed finite y, returns +∞ + i0.
     * <li>If {@code z} is −∞ + i∞, returns +∞ + i3π/4.
     * <li>If {@code z} is +∞ + i∞, returns +∞ + iπ/4.
     * <li>If {@code z} is ±∞ + iNaN, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + iy for finite y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + i∞, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>[1] This has been updated as per
     * <a href="http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1892.htm#dr_471">
     * DR 471: Complex math functions cacosh and ctanh</a>.
     *
     * <p>The inverse hyperbolic cosine is a multivalued function and requires a branch cut in
     * the complex plane; the cut is conventionally placed at the line segment
     * \( (-\infty,-1) \) of the real axis.
     *
     * <p>This function is computed using the trigonomic identity:
     *
     * <p>\[ \cosh^{-1}(z) = \pm i \cos^{-1}(z) \]
     *
     * <p>The sign of the multiplier is chosen to give {@code z.acosh().real() >= 0}
     * and compatibility with the C99 standard.
     *
     * @return The inverse hyperbolic cosine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/ArcCosh/">ArcCosh</a>
     */
    public Complex acosh() {
        // Define in terms of acos
        // acosh(z) = +-i acos(z)
        // Note the special case:
        // acos(+-0 + iNaN) = π/2 + iNaN
        // acosh(0 + iNaN) = NaN + iπ/2
        // will not appropriately multiply by I to maintain positive imaginary if
        // acos() imaginary computes as NaN. So do this explicitly.
        if (Double.isNaN(imaginary) && real == 0) {
            return new Complex(Double.NaN, PI_OVER_2);
        }
        return acos(real, imaginary, (re, im) ->
            // Set the sign appropriately for real >= 0
            (negative(im)) ?
                // Multiply by I
                new Complex(-im, re) :
                // Multiply by -I
                new Complex(im, -re)
        );
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/Cosine.html">
     * cosine</a> of this complex number.
     *
     * <p>\[ \cos(z) = \frac{1}{2} \left( e^{iz} + e^{-iz} \right) \]
     *
     * <p>This is an even function: \( \cos(z) = \cos(-z) \).
     * The cosine is an entire function and requires no branch cuts.
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \cos(x + iy) = \cos(x)\cosh(y) - i \sin(x)\sinh(y) \]
     *
     * <p>As per the C99 standard this function is computed using the trigonomic identity:
     *
     * <p>\[ cos(z) = cosh(iz) \]
     *
     * @return The cosine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Cos/">Cos</a>
     */
    public Complex cos() {
        // Define in terms of cosh
        // cos(z) = cosh(iz)
        // Multiply this number by I and compute cosh.
        return cosh(-imaginary, real, Complex::ofCartesian);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/HyperbolicCosine.html">
     * hyperbolic cosine</a> of this complex number.
     *
     * <p>\[ \cosh(z) = \frac{1}{2} \left( e^{z} + e^{-z} \right) \]
     *
     * <p>The hyperbolic cosine of \( z \) is an entire function in the complex plane
     * and is periodic with respect to the imaginary component with period \( 2\pi i \).
     *
     * <ul>
     * <li>{@code z.conj().cosh() == z.cosh().conj()}.
     * <li>This is an even function: \( \cosh(z) = \cosh(-z) \).
     * <li>If {@code z} is +0 + i0, returns 1 + i0.
     * <li>If {@code z} is +0 + i∞, returns NaN ± i0 (where the sign of the imaginary part of the result is unspecified; "invalid" floating-point operation).
     * <li>If {@code z} is +0 + iNaN, returns NaN ± i0 (where the sign of the imaginary part of the result is unspecified).
     * <li>If {@code z} is x + i∞ for finite nonzero x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is x + iNaN for finite nonzero x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + i0, returns +∞ + i0.
     * <li>If {@code z} is +∞ + iy for finite nonzero y, returns +∞ cis(y) (see {@link #ofCis(double)}).
     * <li>If {@code z} is +∞ + i∞, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is +∞ + iNaN, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + i0, returns NaN ± i0 (where the sign of the imaginary part of the result is unspecified).
     * <li>If {@code z} is NaN + iy for all nonzero numbers y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \cosh(x + iy) = \cosh(x)\cos(y) + i \sinh(x)\sin(y) \]
     *
     * @return The hyperbolic cosine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Cosh/">Cosh</a>
     */
    public Complex cosh() {
        return cosh(real, imaginary, Complex::ofCartesian);
    }

    /**
     * Returns the hyperbolic cosine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code cos(z) = cosh(iz)}.<p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param constructor Constructor.
     * @return The hyperbolic cosine of the complex number.
     */
    private static Complex cosh(double real, double imaginary, ComplexConstructor constructor) {
        // ISO C99: Preserve the even function by mapping to positive
        // f(z) = f(-z)
        if (Double.isInfinite(real) && !Double.isFinite(imaginary)) {
            return constructor.create(Math.abs(real), Double.NaN);
        }
        if (real == 0) {
            // Imaginary-only cosh(iy) = cos(y).
            if (Double.isFinite(imaginary)) {
                // Maintain periodic property with respect to the imaginary component.
                // sinh(+/-0.0) * sin(+/-x) = +/-0 * sin(x)
                return constructor.create(Math.cos(imaginary),
                                          changeSign(real, Math.sin(imaginary)));
            }
            // If imaginary is inf/NaN the sign of the imaginary part is unspecified.
            // Although not required by C99 changing the sign maintains the conjugate equality.
            // It is not possible to also maintain the even function (hence the unspecified sign).
            return constructor.create(Double.NaN, changeSign(real, imaginary));
        }
        if (imaginary == 0) {
            // Real-only cosh(x).
            // Change sign to preserve conjugate equality and even function.
            // sin(+/-0) * sinh(+/-x) = +/-0 * +/-a (sinh is monotonic and same sign)
            // => change the sign of imaginary using real. Handles special case of infinite real.
            // If real is NaN the sign of the imaginary part is unspecified.
            return constructor.create(Math.cosh(real), changeSign(imaginary, real));
        }
        final double x = Math.abs(real);
        if (x > SAFE_EXP) {
            // Approximate sinh/cosh(x) using exp^|x| / 2
            return coshsinh(x, real, imaginary, false, constructor);
        }
        // No overflow of sinh/cosh
        return constructor.create(Math.cosh(real) * Math.cos(imaginary),
                                  Math.sinh(real) * Math.sin(imaginary));
    }

    /**
     * Compute cosh or sinh when the absolute real component |x| is large. In this case
     * cosh(x) and sinh(x) can be approximated by exp(|x|) / 2:
     *
     * <pre>
     * cosh(x+iy) real = (e^|x| / 2) * cos(y)
     * cosh(x+iy) imag = (e^|x| / 2) * sin(y) * sign(x)
     * sinh(x+iy) real = (e^|x| / 2) * cos(y) * sign(x)
     * sinh(x+iy) imag = (e^|x| / 2) * sin(y)
     * </pre>
     *
     * @param x Absolute real component |x|.
     * @param real Real part (x).
     * @param imaginary Imaginary part (y).
     * @param sinh Set to true to compute sinh, otherwise cosh.
     * @param constructor Constructor.
     * @return The hyperbolic sine/cosine of the complex number.
     */
    private static Complex coshsinh(double x, double real, double imaginary, boolean sinh,
                                    ComplexConstructor constructor) {
        // Always require the cos and sin.
        double re = Math.cos(imaginary);
        double im = Math.sin(imaginary);
        // Compute the correct function
        if (sinh) {
            re = changeSign(re, real);
        } else {
            im = changeSign(im, real);
        }
        // Multiply by (e^|x| / 2).
        // Overflow safe computation since sin/cos can be very small allowing a result
        // when e^x overflows: e^x / 2 = (e^m / 2) * e^m * e^(x-2m)
        if (x > SAFE_EXP * 3) {
            // e^x > e^m * e^m * e^m
            // y * (e^m / 2) * e^m * e^m will overflow when starting with Double.MIN_VALUE.
            // Note: Do not multiply by +inf to safeguard against sin(y)=0.0 which
            // will create 0 * inf = nan.
            re *= Double.MAX_VALUE * Double.MAX_VALUE * Double.MAX_VALUE;
            im *= Double.MAX_VALUE * Double.MAX_VALUE * Double.MAX_VALUE;
        } else {
            // Initial part of (e^x / 2) using (e^m / 2)
            re *= EXP_M / 2;
            im *= EXP_M / 2;
            double xm;
            if (x > SAFE_EXP * 2) {
                // e^x = e^m * e^m * e^(x-2m)
                re *= EXP_M;
                im *= EXP_M;
                xm = x - SAFE_EXP * 2;
            } else {
                // e^x = e^m * e^(x-m)
                xm = x - SAFE_EXP;
            }
            final double exp = Math.exp(xm);
            re *= exp;
            im *= exp;
        }
        return constructor.create(re, im);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/ExponentialFunction.html">
     * exponential function</a> of this complex number.
     *
     * <p>\[ \exp(z) = e^z \]
     *
     * <p>The exponential function of \( z \) is an entire function in the complex plane.
     * Special cases:
     *
     * <ul>
     * <li>{@code z.conj().exp() == z.exp().conj()}.
     * <li>If {@code z} is ±0 + i0, returns 1 + i0.
     * <li>If {@code z} is x + i∞ for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is x + iNaN for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + i0, returns +∞ + i0.
     * <li>If {@code z} is −∞ + iy for finite y, returns +0 cis(y) (see {@link #ofCis(double)}).
     * <li>If {@code z} is +∞ + iy for finite nonzero y, returns +∞ cis(y).
     * <li>If {@code z} is −∞ + i∞, returns ±0 ± i0 (where the signs of the real and imaginary parts of the result are unspecified).
     * <li>If {@code z} is +∞ + i∞, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified; "invalid" floating-point operation).
     * <li>If {@code z} is −∞ + iNaN, returns ±0 ± i0 (where the signs of the real and imaginary parts of the result are unspecified).
     * <li>If {@code z} is +∞ + iNaN, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is NaN + i0, returns NaN + i0.
     * <li>If {@code z} is NaN + iy for all nonzero numbers y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>Implements the formula:
     *
     * <p>\[ \exp(x + iy) = e^x (\cos(y) + i \sin(y)) \]
     *
     * @return <code>e<sup>this</sup></code>.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Exp/">Exp</a>
     */
    public Complex exp() {
        if (Double.isInfinite(real)) {
            // Set the scale factor applied to cis(y)
            double zeroOrInf;
            if (real < 0) {
                if (!Double.isFinite(imaginary)) {
                    // (−∞ + i∞) or (−∞ + iNaN) returns (±0 ± i0) (where the signs of the
                    // real and imaginary parts of the result are unspecified).
                    // Here we preserve the conjugate equality.
                    return new Complex(0, Math.copySign(0, imaginary));
                }
                // (−∞ + iy) returns +0 cis(y), for finite y
                zeroOrInf = 0;
            } else {
                // (+∞ + i0) returns +∞ + i0.
                if (imaginary == 0) {
                    return this;
                }
                // (+∞ + i∞) or (+∞ + iNaN) returns (±∞ + iNaN) and raises the invalid
                // floating-point exception (where the sign of the real part of the
                // result is unspecified).
                if (!Double.isFinite(imaginary)) {
                    return new Complex(real, Double.NaN);
                }
                // (+∞ + iy) returns (+∞ cis(y)), for finite nonzero y.
                zeroOrInf = real;
            }
            return new Complex(zeroOrInf * Math.cos(imaginary),
                               zeroOrInf * Math.sin(imaginary));
        } else if (Double.isNaN(real)) {
            // (NaN + i0) returns (NaN + i0)
            // (NaN + iy) returns (NaN + iNaN) and optionally raises the invalid floating-point exception
            // (NaN + iNaN) returns (NaN + iNaN)
            return imaginary == 0 ? this : NAN;
        } else if (!Double.isFinite(imaginary)) {
            // (x + i∞) or (x + iNaN) returns (NaN + iNaN) and raises the invalid
            // floating-point exception, for finite x.
            return NAN;
        }
        // real and imaginary are finite.
        // Compute e^a * (cos(b) + i sin(b)).

        // Special case:
        // (±0 + i0) returns (1 + i0)
        final double exp = Math.exp(real);
        if (imaginary == 0) {
            return new Complex(exp, imaginary);
        }
        return new Complex(exp * Math.cos(imaginary),
                           exp * Math.sin(imaginary));
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/NaturalLogarithm.html">
     * natural logarithm</a> of this complex number.
     *
     * <p>The natural logarithm of \( z \) is unbounded along the real axis and
     * in the range \( [-\pi, \pi] \) along the imaginary axis. The imaginary part of the
     * natural logarithm has a branch cut along the negative real axis \( (-infty,0] \).
     * Special cases:
     *
     * <ul>
     * <li>{@code z.conj().log() == z.log().conj()}.
     * <li>If {@code z} is −0 + i0, returns −∞ + iπ ("divide-by-zero" floating-point operation).
     * <li>If {@code z} is +0 + i0, returns −∞ + i0 ("divide-by-zero" floating-point operation).
     * <li>If {@code z} is x + i∞ for finite x, returns +∞ + iπ/2.
     * <li>If {@code z} is x + iNaN for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is −∞ + iy for finite positive-signed y, returns +∞ + iπ.
     * <li>If {@code z} is +∞ + iy for finite positive-signed y, returns +∞ + i0.
     * <li>If {@code z} is −∞ + i∞, returns +∞ + i3π/4.
     * <li>If {@code z} is +∞ + i∞, returns +∞ + iπ/4.
     * <li>If {@code z} is ±∞ + iNaN, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + iy for finite y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + i∞, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>Implements the formula:
     *
     * <p>\[ \ln(z) = \ln |z| + i \arg(z) \]
     *
     * <p>where \( |z| \) is the absolute and \( \arg(z) \) is the argument.
     *
     * <p>The implementation is based on the method described in:</p>
     * <blockquote>
     * T E Hull, Thomas F Fairgrieve and Ping Tak Peter Tang (1994)
     * Implementing complex elementary functions using exception handling.
     * ACM Transactions on Mathematical Software, Vol 20, No 2, pp 215-244.
     * </blockquote>
     *
     * @return The natural logarithm of this complex number.
     * @see Math#log(double)
     * @see #abs()
     * @see #arg()
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Log/">Log</a>
     */
    public Complex log() {
        return log(Math::log, HALF, LN_2, Complex::ofCartesian);
    }

    /**
     * Returns the base 10
     * <a href="http://mathworld.wolfram.com/CommonLogarithm.html">
     * common logarithm</a> of this complex number.
     *
     * <p>The common logarithm of \( z \) is unbounded along the real axis and
     * in the range \( [-\pi, \pi] \) along the imaginary axis. The imaginary part of the
     * common logarithm has a branch cut along the negative real axis \( (-infty,0] \).
     * Special cases are as defined in the {@link #log() natural logarithm}:
     *
     * <p>Implements the formula:
     *
     * <p>\[ \log_{10}(z) = \log_{10} |z| + i \arg(z) \]
     *
     * <p>where \( |z| \) is the absolute and \( \arg(z) \) is the argument.
     *
     * @return The base 10 logarithm of this complex number.
     * @see Math#log10(double)
     * @see #abs()
     * @see #arg()
     */
    public Complex log10() {
        return log(Math::log10, LOG_10E_O_2, LOG10_2, Complex::ofCartesian);
    }

    /**
     * Returns the logarithm of this complex number using the provided function.
     * Implements the formula:
     *
     * <pre>
     *   log(x + i y) = log(|x + i y|) + i arg(x + i y)</pre>
     *
     * <p>Warning: The argument {@code logOf2} must be equal to {@code log(2)} using the
     * provided log function otherwise scaling using powers of 2 in the case of overflow
     * will be incorrect. This is provided as an internal optimisation.
     *
     * @param log Log function.
     * @param logOfeOver2 The log function applied to e, then divided by 2.
     * @param logOf2 The log function applied to 2.
     * @param constructor Constructor for the returned complex.
     * @return The logarithm of this complex number.
     * @see #abs()
     * @see #arg()
     */
    private Complex log(UnaryOperation log, double logOfeOver2, double logOf2, ComplexConstructor constructor) {
        // Handle NaN
        if (Double.isNaN(real) || Double.isNaN(imaginary)) {
            // Return NaN unless infinite
            if (isInfinite()) {
                return constructor.create(Double.POSITIVE_INFINITY, Double.NaN);
            }
            // No-use of the input constructor
            return NAN;
        }

        // Returns the real part:
        // log(sqrt(x^2 + y^2))
        // log(x^2 + y^2) / 2

        // Compute with positive values
        double x = Math.abs(real);
        double y = Math.abs(imaginary);

        // Find the larger magnitude.
        if (x < y) {
            final double tmp = x;
            x = y;
            y = tmp;
        }

        if (x == 0) {
            // Handle zero: raises the ‘‘divide-by-zero’’ floating-point exception.
            return constructor.create(Double.NEGATIVE_INFINITY,
                                      negative(real) ? Math.copySign(Math.PI, imaginary) : imaginary);
        }

        double re;

        if (x > HALF && x < ROOT2) {
            // x^2+y^2 close to 1. Use log1p(x^2+y^2 - 1) / 2.
            re = Math.log1p(x2y2m1(x, y)) * logOfeOver2;
        } else if (y == 0) {
            // Handle real only number
            re = log.apply(x);
        } else if (x > SAFE_MAX || y < SAFE_MIN) {
            // Over/underflow of sqrt(x^2+y^2)
            // Note: Since y<x no check for y > SAFE_MAX or x < SAFE_MIN.
            if (isPosInfinite(x)) {
                // Handle infinity
                re = x;
            } else {
                // Do scaling
                final int expx = Math.getExponent(x);
                final int expy = Math.getExponent(y);
                // Hull et al: 2 * (expx - expy) > precision + 1
                // Modified to use the pre-computed half precision
                if ((expx - expy) > HALF_PRECISION) {
                    // y can be ignored
                    re = log.apply(x);
                } else {
                    // Hull et al:
                    // "It is important that the scaling be chosen so
                    // that there is no possibility of cancellation in this addition"
                    // i.e. sx^2 + sy^2 should not be close to 1.
                    // Their paper uses expx + 2 for underflow but expx for overflow.
                    // It has been modified here to use expx - 2.
                    int scale;
                    if (x > SAFE_MAX) {
                        // overflow
                        scale = expx - 2;
                    } else {
                        // underflow
                        scale = expx + 2;
                    }
                    final double sx = Math.scalb(x, -scale);
                    final double sy = Math.scalb(y, -scale);
                    re = scale * logOf2 + 0.5 * log.apply(sx * sx + sy * sy);
                }
            }
        } else {
            // Safe region that avoids under/overflow
            re = 0.5 * log.apply(x * x + y * y);
        }

        // All ISO C99 edge cases for the imaginary are satisfied by the Math library.
        return constructor.create(re, arg());
    }

    /**
     * Compute {@code x^2 + y^2 - 1} in high precision.
     * Assumes that the values x and y can be multiplied without overflow; that
     * {@code x >= y}; and both values are positive.
     *
     * @param x the x value
     * @param y the y value
     * @return {@code x^2 + y^2 - 1}.
     */
    private static double x2y2m1(double x, double y) {
        // Hull et al used (x-1)*(x+1)+y*y.
        // From the paper on page 236:

        // If x == 1 there is no cancellation.

        // If x > 1, there is also no cancellation, but the argument is now accurate
        // only to within a factor of 1 + 3 EPSILSON (note that x – 1 is exact),
        // so that error = 3 EPSILON.

        // If x < 1, there can be serious cancellation:

        // If 4 y^2 < |x^2 – 1| the cancellation is not serious ... the argument is accurate
        // only to within a factor of 1 + 4 EPSILSON so that error = 4 EPSILON.

        // Otherwise there can be serious cancellation and the relative error in the real part
        // could be enormous.

        final double xx = x * x;
        final double yy = y * y;
        // Modify to use high precision before the threshold set by Hull et al.
        // This is to preserve the monotonic output of the computation at the switch.
        // Set the threshold when x^2 + y^2 is above 0.5 thus subtracting 1 results in a number
        // that can be expressed with a higher precision than any number in the range 0.5-1.0
        // due to the variable exponent used below 0.5.
        if (x < 1 && xx + yy > 0.5) {
            // Large relative error.
            // This does not use o.a.c.numbers.LinearCombination.value(x, x, y, y, 1, -1).
            // It is optimised knowing that:
            // - the products are squares
            // - the final term is -1 (which does not require split multiplication and addition)
            // - The answer will not be NaN as the terms are not NaN components
            // - The order is known to be 1 > |x| >= |y|
            // The squares are computed using a split multiply algorithm and
            // the summation using an extended precision summation algorithm.

            // Split x and y as one 26 bits number and one 27 bits number
            final double xHigh = splitHigh(x);
            final double xLow  = x - xHigh;
            final double yHigh = splitHigh(y);
            final double yLow  = y - yHigh;

            // Accurate split multiplication x * x and y * y
            final double x2Low = squareLow(xLow, xHigh, xx);
            final double y2Low = squareLow(yLow, yHigh, yy);

            return sumx2y2m1(xx, x2Low, yy, y2Low);
        }
        return (x - 1) * (x + 1) + yy;
    }

    /**
     * Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) create
     * a big value from which to derive the two split parts.
     * <pre>
     * c = (2^s + 1) * a
     * a_big = c - a
     * a_hi = c - a_big
     * a_lo = a - a_hi
     * a = a_hi + a_lo
     * </pre>
     *
     * <p>The multiplicand must be odd allowing a p-bit value to be split into
     * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
     * Combined they have (p􏰔-1) bits of significand but the sign bit of {@code a_lo}
     * contains a bit of information.
     *
     * @param a Value.
     * @return the high part of the value.
     * @see <a href="https://doi.org/10.1007/BF01397083">
     * Dekker (1971) A floating-point technique for extending the available precision</a>
     */
    private static double splitHigh(double a) {
        final double c = MULTIPLIER * a;
        return c - (c - a);
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
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Dekker's two-sum algorithm. The values are required to be ordered by magnitude:
     * {@code |a| >= |b|}.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return <code>b - (x - a)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 6</a>
     */
    private static double fastSumLow(double a, double b, double x) {
        // x = a + b
        // bVirtual = x - a
        // y = b - bVirtual
        return b - (x - a);
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Knuth's two-sum algorithm. The values are not required to be ordered by magnitude.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return <code>(a - (x - (x - a))) + (b - (x - a))</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 7</a>
     */
    private static double sumLow(double a, double b, double x) {
        // x = a + b
        // bVirtual = x - a
        // aVirtual = x - bVirtual
        // bRoundoff = b - bVirtual
        // aRoundoff = a - aVirtual
        // y = aRoundoff + bRoundoff
        final double bVirtual = x - a;
        return (a - (x - bVirtual)) + (b - bVirtual);
    }

    /**
     * Sum x^2 + y^2 - 1. It is assumed that {@code y <= x < 1}.
     *
     * <p>Implement Shewchuk's expansion-sum algorithm: [x2Low, x2High, -1] + [y2Low, y2High].
     *
     * @param x2High High part of x^2.
     * @param x2Low Low part of x^2.
     * @param y2High High part of y^2.
     * @param y2Low Low part of y^2.
     * @return x^2 + y^2 - 1
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 12</a>
     */
    private static double sumx2y2m1(double x2High, double x2Low, double y2High, double y2Low) {
        // Let e and f be non-overlapping expansions of components of length m and n.
        // The following algorithm will produce a non-overlapping expansion h where the
        // sum h_i = e + f and components of h are in increasing order of magnitude.

        // Expansion sum proceeds by a grow-expansion of the first part from one expansion
        // into the other, extending its length by 1. The process repeats for the next part
        // but the grow expansion starts at the previous merge position + 1.
        // Thus expansion sum requires mn two-sum operations to merge length m into length n
        // resulting in length m+n-1.

        // Variables numbered from 1 as per Figure 7 (p.12). The output expansion h is placed
        // into e increasing its length for each grow expansion.

        // We have two expansions for x^2 and y^2 and the whole number -1.
        // Expecting (x^2 + y^2) close to 1 we generate first the intermediate expansion
        // (x^2 - 1) moving the result away from 1 where there are sparse floating point
        // representations. This is then added to a similar magnitude y^2.

        // q=running sum
        double q = x2Low - 1;
        double e1 = fastSumLow(-1, x2Low, q);
        double e3 = q + x2High;
        double e2 = sumLow(q, x2High, e3);

        final double f1 = y2Low;
        final double f2 = y2High;

        // Grow expansion of f1 into e
        q = f1 + e1;
        e1 = sumLow(f1, e1, q);
        double p = q;
        q += e2;
        e2 = sumLow(p, e2, q);
        double e4 = q + e3;
        e3 = sumLow(q, e3, e4);

        // Grow expansion of f2 into e (only required to start at e2)
        q = f2 + e2;
        e2 = sumLow(f2, e2, q);
        p = q;
        q += e3;
        e3 = sumLow(p, e3, q);
        final double e5 = q + e4;
        e4 = sumLow(q, e4, e5);

        // Final summation
        return e1 + e2 + e3 + e4 + e5;
    }

    /**
     * Returns the complex power of this complex number raised to the power of \( x \).
     * Implements the formula:
     *
     * <p>\[ z^x = e^{x \ln(z)} \]
     *
     * <p>If this complex number is zero then this method returns zero if \( x \) is positive
     * in the real component and zero in the imaginary component;
     * otherwise it returns NaN + iNaN.
     *
     * @param  x The exponent to which this complex number is to be raised.
     * @return <code>this<sup>x</sup></code>.
     * @see #log()
     * @see #multiply(Complex)
     * @see #exp()
     * @see <a href="http://mathworld.wolfram.com/ComplexExponentiation.html">Complex exponentiation</a>
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Power/">Power</a>
     */
    public Complex pow(Complex x) {
        if (real == 0 &&
            imaginary == 0) {
            // This value is zero. Test the other.
            if (x.real > 0 &&
                x.imaginary == 0) {
                // 0 raised to positive number is 0
                return ZERO;
            }
            // 0 raised to anything else is NaN
            return NAN;
        }
        return log().multiply(x).exp();
    }

    /**
     * Returns the complex power of this complex number raised to the power of \( x \).
     * Implements the formula:
     *
     * <p>\[ z^x = e^{x \ln(z)} \]
     *
     * <p>If this complex number is zero then this method returns zero if \( x \) is positive;
     * otherwise it returns NaN + iNaN.
     *
     * @param  x The exponent to which this complex number is to be raised.
     * @return <code>this<sup>x</sup></code>.
     * @see #log()
     * @see #multiply(double)
     * @see #exp()
     * @see #pow(Complex)
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Power/">Power</a>
     */
    public Complex pow(double x) {
        if (real == 0 &&
            imaginary == 0) {
            // This value is zero. Test the other.
            if (x > 0) {
                // 0 raised to positive number is 0
                return ZERO;
            }
            // 0 raised to anything else is NaN
            return NAN;
        }
        return log().multiply(x).exp();
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/Sine.html">
     * sine</a> of this complex number.
     *
     * <p>\[ \sin(z) = \frac{1}{2} i \left( e^{-iz} - e^{iz} \right) \]
     *
     * <p>This is an odd function: \( \sin(z) = -\sin(-z) \).
     * The sine is an entire function and requires no branch cuts.
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \sin(x + iy) = \sin(x)\cosh(y) + i \cos(x)\sinh(y) \]
     *
     * <p>As per the C99 standard this function is computed using the trigonomic identity:
     *
     * <p>\[ \sin(z) = -i \sinh(iz) \]
     *
     * @return The sine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Sin/">Sin</a>
     */
    public Complex sin() {
        // Define in terms of sinh
        // sin(z) = -i sinh(iz)
        // Multiply this number by I, compute sinh, then multiply by back
        return sinh(-imaginary, real, Complex::multiplyNegativeI);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/HyperbolicSine.html">
     * hyperbolic sine</a> of this complex number.
     *
     * <p>\[ \sinh(z) = \frac{1}{2} \left( e^{z} - e^{-z} \right) \]
     *
     * <p>The hyperbolic sine of \( z \) is an entire function in the complex plane
     * and is periodic with respect to the imaginary component with period \( 2\pi i \).
     *
     * <ul>
     * <li>{@code z.conj().sinh() == z.sinh().conj()}.
     * <li>This is an odd function: \( \sinh(z) = -\sinh(-z) \).
     * <li>If {@code z} is +0 + i0, returns +0 + i0.
     * <li>If {@code z} is +0 + i∞, returns ±0 + iNaN (where the sign of the real part of the result is unspecified; "invalid" floating-point operation).
     * <li>If {@code z} is +0 + iNaN, returns ±0 + iNaN (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is x + i∞ for positive finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is x + iNaN for finite nonzero x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + i0, returns +∞ + i0.
     * <li>If {@code z} is +∞ + iy for positive finite y, returns +∞ cis(y) (see {@link #ofCis(double)}.
     * <li>If {@code z} is +∞ + i∞, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified; "invalid" floating-point operation).
     * <li>If {@code z} is +∞ + iNaN, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is NaN + i0, returns NaN + i0.
     * <li>If {@code z} is NaN + iy for all nonzero numbers y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \sinh(x + iy) = \sinh(x)\cos(y) + i \cosh(x)\sin(y) \]
     *
     * @return The hyperbolic sine of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Sinh/">Sinh</a>
     */
    public Complex sinh() {
        return sinh(real, imaginary, Complex::ofCartesian);
    }

    /**
     * Returns the hyperbolic sine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code sin(z) = -i sinh(iz)}.<p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param constructor Constructor.
     * @return The hyperbolic sine of the complex number.
     */
    private static Complex sinh(double real, double imaginary, ComplexConstructor constructor) {
        if (Double.isInfinite(real) && !Double.isFinite(imaginary)) {
            return constructor.create(real, Double.NaN);
        }
        if (real == 0) {
            // Imaginary-only sinh(iy) = i sin(y).
            if (Double.isFinite(imaginary)) {
                // Maintain periodic property with respect to the imaginary component.
                // sinh(+/-0.0) * cos(+/-x) = +/-0 * cos(x)
                return constructor.create(changeSign(real, Math.cos(imaginary)),
                                          Math.sin(imaginary));
            }
            // If imaginary is inf/NaN the sign of the real part is unspecified.
            // Returning the same real value maintains the conjugate equality.
            // It is not possible to also maintain the odd function (hence the unspecified sign).
            return constructor.create(real, Double.NaN);
        }
        if (imaginary == 0) {
            // Real-only sinh(x).
            return constructor.create(Math.sinh(real), imaginary);
        }
        final double x = Math.abs(real);
        if (x > SAFE_EXP) {
            // Approximate sinh/cosh(x) using exp^|x| / 2
            return coshsinh(x, real, imaginary, true, constructor);
        }
        // No overflow of sinh/cosh
        return constructor.create(Math.sinh(real) * Math.cos(imaginary),
                                  Math.cosh(real) * Math.sin(imaginary));
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html">
     * square root</a> of this complex number.
     *
     * <p>\[ \sqrt{x + iy} = \frac{1}{2} \sqrt{2} \left( \sqrt{ \sqrt{x^2 + y^2} + x } + i\ \text{sgn}(y) \sqrt{ \sqrt{x^2 + y^2} - x } \right) \]
     *
     * <p>The square root of \( z \) is in the range \( [0, +\infty) \) along the real axis and
     * is unbounded along the imaginary axis. The imaginary part of the square root has a
     * branch cut along the negative real axis \( (-infty,0) \). Special cases:
     *
     * <ul>
     * <li>{@code z.conj().sqrt() == z.sqrt().conj()}.
     * <li>If {@code z} is ±0 + i0, returns +0 + i0.
     * <li>If {@code z} is x + i∞ for all x (including NaN), returns +∞ + i∞.
     * <li>If {@code z} is x + iNaN for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is −∞ + iy for finite positive-signed y, returns +0 + i∞.
     * <li>If {@code z} is +∞ + iy for finite positive-signed y, returns +∞ + i0.
     * <li>If {@code z} is −∞ + iNaN, returns NaN ± i∞ (where the sign of the imaginary part of the result is unspecified).
     * <li>If {@code z} is +∞ + iNaN, returns +∞ + iNaN.
     * <li>If {@code z} is NaN + iy for finite y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>Implements the following algorithm to compute \( \sqrt{x + iy} \):
     * <ol>
     * <li>Let \( t = \sqrt{2 (|x| + |x + iy|)} \)
     * <li>if \( x \geq 0 \) return \( \frac{t}{2} + i \frac{y}{t} \)
     * <li>else return \( \frac{|y|}{t} + i\ \text{sgn}(y) \frac{t}{2} \)
     * </ol>
     * where:
     * <ul>
     * <li>\( |x| =\ \){@link Math#abs(double) abs}(x)
     * <li>\( |x + y i| =\ \){@link Complex#abs}
     * <li>\( \text{sgn}(y) =\ \){@link Math#copySign(double,double) copySign}(1.0, y)
     * </ul>
     *
     * <p>The implementation is overflow and underflow safe based on the method described in:</p>
     * <blockquote>
     * T E Hull, Thomas F Fairgrieve and Ping Tak Peter Tang (1994)
     * Implementing complex elementary functions using exception handling.
     * ACM Transactions on Mathematical Software, Vol 20, No 2, pp 215-244.
     * </blockquote>
     *
     * @return The square root of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Sqrt/">Sqrt</a>
     */
    public Complex sqrt() {
        return sqrt(real, imaginary);
    }

    /**
     * Returns the square root of the complex number {@code sqrt(x + i y)}.
     *
     * @param real Real component.
     * @param imaginary Imaginary component.
     * @return The square root of the complex number.
     */
    private static Complex sqrt(double real, double imaginary) {
        // Handle NaN
        if (Double.isNaN(real) || Double.isNaN(imaginary)) {
            // Check for infinite
            if (Double.isInfinite(imaginary)) {
                return new Complex(Double.POSITIVE_INFINITY, imaginary);
            }
            if (Double.isInfinite(real)) {
                if (real == Double.NEGATIVE_INFINITY) {
                    return new Complex(Double.NaN, Math.copySign(Double.POSITIVE_INFINITY, imaginary));
                }
                return new Complex(Double.POSITIVE_INFINITY, Double.NaN);
            }
            return NAN;
        }

        // Compute with positive values and determine sign at the end
        final double x = Math.abs(real);
        final double y = Math.abs(imaginary);

        // Compute
        double t;

        if (inRegion(x, y, SAFE_MIN, SAFE_MAX)) {
            // No over/underflow of x^2 + y^2
            t = Math.sqrt(2 * (Math.sqrt(x * x + y * y) + x));
        } else {
            // Potential over/underflow. First check infinites and real/imaginary only.

            // Check for infinite
            if (isPosInfinite(y)) {
                return new Complex(Double.POSITIVE_INFINITY, imaginary);
            } else if (isPosInfinite(x)) {
                if (real == Double.NEGATIVE_INFINITY) {
                    return new Complex(0, Math.copySign(Double.POSITIVE_INFINITY, imaginary));
                }
                return new Complex(Double.POSITIVE_INFINITY, Math.copySign(0, imaginary));
            } else if (y == 0) {
                // Real only
                final double sqrtAbs = Math.sqrt(x);
                if (real < 0) {
                    return new Complex(0, Math.copySign(sqrtAbs, imaginary));
                }
                return new Complex(sqrtAbs, imaginary);
            } else if (x == 0) {
                // Imaginary only
                final double sqrtAbs = Math.sqrt(y) / ROOT2;
                return new Complex(sqrtAbs, Math.copySign(sqrtAbs, imaginary));
            } else {
                // Over/underflow
                // scale so that abs(x) is near 1, with even exponent.
                final int scale = getScale(x, y) & MASK_INT_TO_EVEN;
                final double sx = Math.scalb(x, -scale);
                final double sy = Math.scalb(y, -scale);
                final double st = Math.sqrt(2 * (Math.sqrt(sx * sx + sy * sy) + sx));
                // Rescale. This works if exponent is even:
                // st * sqrt(2^scale) = st * (2^scale)^0.5 = st * 2^(scale*0.5)
                t = Math.scalb(st, scale / 2);
            }
        }

        if (real >= 0) {
            return new Complex(t / 2, imaginary / t);
        }
        return new Complex(y / t, Math.copySign(t / 2, imaginary));
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/Tangent.html">
     * tangent</a> of this complex number.
     *
     * <p>\[ \tan(z) = \frac{i(e^{-iz} - e^{iz})}{e^{-iz} + e^{iz}} \]
     *
     * <p>This is an odd function: \( \tan(z) = -\tan(-z) \).
     * The tangent is an entire function and requires no branch cuts.
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:</p>
     * \[ \tan(x + iy) = \frac{\sin(2x)}{\cos(2x)+\cosh(2y)} + i \frac{\sinh(2y)}{\cos(2x)+\cosh(2y)} \]
     *
     * <p>As per the C99 standard this function is computed using the trigonomic identity:</p>
     * \[ \tan(z) = -i \tanh(iz) \]
     *
     * @return The tangent of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Tan/">Tangent</a>
     */
    public Complex tan() {
        // Define in terms of tanh
        // tan(z) = -i tanh(iz)
        // Multiply this number by I, compute tanh, then multiply by back
        return tanh(-imaginary, real, Complex::multiplyNegativeI);
    }

    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/HyperbolicTangent.html">
     * hyperbolic tangent</a> of this complex number.
     *
     * <p>\[ \tanh(z) = \frac{e^z - e^{-z}}{e^z + e^{-z}} \]
     *
     * <p>The hyperbolic tangent of \( z \) is an entire function in the complex plane
     * and is periodic with respect to the imaginary component with period \( \pi i \)
     * and has poles of the first order along the imaginary line, at coordinates
     * \( (0, \pi(\frac{1}{2} + n)) \).
     * Note that the {@code double} floating-point representation is unable to exactly represent
     * \( \pi/2 \) and there is no value for which a pole error occurs.
     *
     * <ul>
     * <li>{@code z.conj().tanh() == z.tanh().conj()}.
     * <li>This is an odd function: \( \tanh(z) = -\tanh(-z) \).
     * <li>If {@code z} is +0 + i0, returns +0 + i0.
     * <li>If {@code z} is 0 + i∞, returns 0 + iNaN.
     * <li>If {@code z} is x + i∞ for finite non-zero x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is 0 + iNaN, returns 0 + iNAN.
     * <li>If {@code z} is x + iNaN for finite non-zero x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + iy for positive-signed finite y, returns 1 + i0 sin(2y).
     * <li>If {@code z} is +∞ + i∞, returns 1 ± i0 (where the sign of the imaginary part of the result is unspecified).
     * <li>If {@code z} is +∞ + iNaN, returns 1 ± i0 (where the sign of the imaginary part of the result is unspecified).
     * <li>If {@code z} is NaN + i0, returns NaN + i0.
     * <li>If {@code z} is NaN + iy for all nonzero numbers y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>[1] This has been updated as per
     * <a href="http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1892.htm#dr_471">
     * DR 471: Complex math functions cacosh and ctanh</a>.
     *
     * <p>This is implemented using real \( x \) and imaginary \( y \) parts:
     *
     * <p>\[ \tan(x + iy) = \frac{\sinh(2x)}{\cosh(2x)+\cos(2y)} + i \frac{\sin(2y)}{\cosh(2x)+\cos(2y)} \]
     *
     * @return The hyperbolic tangent of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Tanh/">Tanh</a>
     */
    public Complex tanh() {
        return tanh(real, imaginary, Complex::ofCartesian);
    }

    /**
     * Returns the hyperbolic tangent of this complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code tan(z) = -i tanh(iz)}.<p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param constructor Constructor.
     * @return The hyperbolic tangent of the complex number.
     */
    private static Complex tanh(double real, double imaginary, ComplexConstructor constructor) {
        // Cache the absolute real value
        final double x = Math.abs(real);

        // Handle inf or nan.
        // Deliberate logic inversion using x to match !Double.isFinite(x) knowing x is absolute.
        if (!(x <= Double.MAX_VALUE) || !Double.isFinite(imaginary)) {
            if (isPosInfinite(x)) {
                if (Double.isFinite(imaginary)) {
                    return constructor.create(Math.copySign(1, real),
                                              Math.copySign(0, sin2(imaginary)));
                }
                // imaginary is infinite or NaN
                return constructor.create(Math.copySign(1, real), Math.copySign(0, imaginary));
            }
            // Remaining cases:
            // (0 + i inf), returns (0 + i NaN)
            // (0 + i NaN), returns (0 + i NaN)
            // (x + i inf), returns (NaN + i NaN) for non-zero x (including infinite)
            // (x + i NaN), returns (NaN + i NaN) for non-zero x (including infinite)
            // (NaN + i 0), returns (NaN + i 0)
            // (NaN + i y), returns (NaN + i NaN) for non-zero y (including infinite)
            // (NaN + i NaN), returns (NaN + i NaN)
            return constructor.create(real == 0 ? real : Double.NaN,
                                      imaginary == 0 ? imaginary : Double.NaN);
        }

        // Finite components
        // tanh(x+iy) = (sinh(2x) + i sin(2y)) / (cosh(2x) + cos(2y))

        if (real == 0) {
            // Imaginary-only tanh(iy) = i tan(y)
            // Identity: sin 2y / (1 + cos 2y) = tan(y)
            return constructor.create(real, Math.tan(imaginary));
        }
        if (imaginary == 0) {
            // Identity: sinh 2x / (1 + cosh 2x) = tanh(x)
            return constructor.create(Math.tanh(real), imaginary);
        }

        // The double angles can be avoided using the identities:
        // sinh(2x) = 2 sinh(x) cosh(x)
        // sin(2y) = 2 sin(y) cos(y)
        // cosh(2x) = 2 sinh^2(x) + 1
        // cos(2y) = 2 cos^2(y) - 1
        // tanh(x+iy) = (sinh(x)cosh(x) + i sin(y)cos(y)) / (sinh^2(x) + cos^2(y))
        // tanh(x+iy) = (sinh(x)cosh(x) + i 0.5 sin(2y)) / (sinh^2(x) + cos^2(y))

        if (x > SAFE_EXP / 2) {
            // Potential overflow in sinh/cosh(2x).
            // Approximate sinh/cosh using exp^x.
            // Ignore cos^2(y) in the divisor as it is insignificant.
            // real = sinh(x)cosh(x) / sinh^2(x) = +/-1
            final double re = Math.copySign(1, real);
            // imag = sin(2y) / 2 sinh^2(x)
            // sinh(x) -> sign(x) * e^|x| / 2 when x is large.
            // sinh^2(x) -> e^2|x| / 4 when x is large.
            // imag = sin(2y) / 2 (e^2|x| / 4) = 2 sin(2y) / e^2|x|
            // Underflow safe divide as e^2|x| may overflow:
            // imag = 2 sin(2y) / e^m / e^(2|x| - m)
            double im = sin2(imaginary);
            if (x > SAFE_EXP) {
                // e^2|x| > e^m * e^m
                // This will underflow 2.0 / e^m / e^m
                im = Math.copySign(0.0, im);
            } else {
                // e^2|x| = e^m * e^(2|x| - m)
                im = 2 * im / EXP_M / Math.exp(2 * x - SAFE_EXP);
            }
            return constructor.create(re, im);
        }

        // No overflow of sinh(2x) and cosh(2x)

        // Note: This does not use the double angle identities and returns the
        // definitional formula when 2y is finite. The transition when 2y overflows
        // and cos(2y) and sin(2y) switch to identities is monotonic at the junction
        // but the function is not smooth due to the sampling
        // of the 2 pi period at very large jumps of x. Thus this returns a value
        // but the usefulness as y -> inf may be limited.

        final double real2 = 2 * real;
        final double divisor = Math.cosh(real2) + cos2(imaginary);
        return constructor.create(Math.sinh(real2) / divisor,
                                  sin2(imaginary) / divisor);
    }

    /**
     * Safely compute {@code cos(2*a)} when {@code a} is finite.
     * Note that {@link Math#cos(double)} returns NaN when the input is infinite.
     * If {@code 2*a} is finite use {@code Math.cos(2*a)}; otherwise use the identity:
     *
     * <pre>
     *  <code>cos(2a) = 2 cos<sup>2</sup>(a) - 1</code></pre>
     *
     * @param a Angle a.
     * @return The cosine of 2a.
     * @see Math#cos(double)
     */
    private static double cos2(double a) {
        final double twoA = 2 * a;
        if (Double.isFinite(twoA)) {
            return Math.cos(twoA);
        }
        final double cosA = Math.cos(a);
        return 2 * cosA * cosA - 1;
    }

    /**
     * Safely compute {@code sin(2*a)} when {@code a} is finite.
     * Note that {@link Math#sin(double)} returns NaN when the input is infinite.
     * If {@code 2*a} is finite use {@code Math.sin(2*a)}; otherwise use the identity:
     *
     * <pre>
     *  <code>sin(2a) = 2 sin(a) cos(a)</code></pre>
     *
     * @param a Angle a.
     * @return The sine of 2a.
     * @see Math#sin(double)
     */
    private static double sin2(double a) {
        final double twoA = 2 * a;
        if (Double.isFinite(twoA)) {
            return Math.sin(twoA);
        }
        return 2 * Math.sin(a) * Math.cos(a);
    }

    /**
     * Returns the argument of this complex number.
     *
     * <p>The argument is the angle phi between the positive real axis and
     * the point representing this number in the complex plane.
     * The value returned is between \( -\pi \) (not inclusive)
     * and \( \pi \) (inclusive), with negative values returned for numbers with
     * negative imaginary parts.
     *
     * <p>If either real or imaginary part (or both) is NaN, NaN is returned.
     * Infinite parts are handled as {@linkplain Math#atan2} handles them,
     * essentially treating finite parts as zero in the presence of an
     * infinite coordinate and returning a multiple of \( \frac{\pi}{4} \) depending on
     * the signs of the infinite parts.
     *
     * <p>This code follows the
     * <a href="http://www.iso-9899.info/wiki/The_Standard">ISO C Standard</a>, Annex G,
     * in calculating the returned value using the {@code atan2(y, x)} method for complex
     * \( x + iy \).
     *
     * @return The argument of this complex number.
     * @see Math#atan2(double, double)
     */
    public double arg() {
        // Delegate
        return Math.atan2(imaginary, real);
    }

    /**
     * Returns the n-th roots of this complex number.
     * The nth roots are defined by the formula:
     *
     * <p>\[ z_k = |z|^{\frac{1}{n}} \left( \cos \left(\phi + \frac{2\pi k}{n} \right) + i \sin \left(\phi + \frac{2\pi k}{n} \right) \right) \]
     *
     * <p>for \( k=0, 1, \ldots, n-1 \), where \( |z| \) and \( \phi \)
     * are respectively the {@link #abs() modulus} and
     * {@link #arg() argument} of this complex number.
     *
     * <p>If one or both parts of this complex number is NaN, a list with all
     * all elements set to {@code NaN + i NaN} is returned.</p>
     *
     * @param n Degree of root.
     * @return A list of all {@code n}-th roots of this complex number.
     * @throws IllegalArgumentException if {@code n} is zero.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Root/">Root</a>
     */
    public List<Complex> nthRoot(int n) {
        if (n == 0) {
            throw new IllegalArgumentException("cannot compute zeroth root");
        }

        final List<Complex> result = new ArrayList<>();

        // nth root of abs -- faster / more accurate to use a solver here?
        final double nthRootOfAbs = Math.pow(abs(), 1.0 / n);

        // Compute nth roots of complex number with k = 0, 1, ... n-1
        final double nthPhi = arg() / n;
        final double slice = 2 * Math.PI / n;
        double innerPart = nthPhi;
        for (int k = 0; k < Math.abs(n); k++) {
            // inner part
            final double realPart = nthRootOfAbs *  Math.cos(innerPart);
            final double imaginaryPart = nthRootOfAbs *  Math.sin(innerPart);
            result.add(new Complex(realPart, imaginaryPart));
            innerPart += slice;
        }

        return result;
    }

    /**
     * Returns a string representation of the complex number.
     *
     * <p>The string will represent the numeric values of the real and imaginary parts.
     * The values are split by a separator and surrounded by parentheses.
     * The string can be {@link #parse(String) parsed} to obtain an instance with the same value.
     *
     * <p>The format for complex number \( x + i y \) is {@code "(x,y)"}, with \( x \) and
     * \( y \) converted as if using {@link Double#toString(double)}.
     *
     * @return A string representation of the complex number.
     * @see #parse(String)
     * @see Double#toString(double)
     */
    @Override
    public String toString() {
        return new StringBuilder(TO_STRING_SIZE)
            .append(FORMAT_START)
            .append(real).append(FORMAT_SEP)
            .append(imaginary)
            .append(FORMAT_END)
            .toString();
    }

    /**
     * Returns {@code true} if the values are equal according to semantics of
     * {@link Double#equals(Object)}.
     *
     * @param x Value
     * @param y Value
     * @return {@code Double.valueof(x).equals(Double.valueOf(y))}.
     */
    private static boolean equals(double x, double y) {
        return Double.doubleToLongBits(x) == Double.doubleToLongBits(y);
    }

    /**
     * Check that a value is negative. It must meet all the following conditions:
     * <ul>
     *  <li>it is not {@code NaN},</li>
     *  <li>it is negative signed,</li>
     * </ul>
     *
     * <p>Note: This is true for negative zero.</p>
     *
     * @param d Value.
     * @return {@code true} if {@code d} is negative.
     */
    private static boolean negative(double d) {
        return d < 0 || Double.doubleToLongBits(d) == NEGATIVE_ZERO_LONG_BITS;
    }

    /**
     * Check that a value is positive infinity. Used to replace {@link Double#isInfinite()}
     * when the input value is known to be positive (i.e. in the case where it has been
     * set using {@link Math#abs(double)}).
     *
     * @param d Value.
     * @return {@code true} if {@code d} is +inf.
     */
    private static boolean isPosInfinite(double d) {
        return d == Double.POSITIVE_INFINITY;
    }

    /**
     * Create a complex number given the real and imaginary parts, then multiply by {@code -i}.
     * This is used in functions that implement trigonomic identities. It is the functional
     * equivalent of:
     *
     * <pre>
     *  z = new Complex(real, imaginary).multiplyImaginary(-1);</pre>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @return {@code Complex} object.
     */
    private static Complex multiplyNegativeI(double real, double imaginary) {
        return new Complex(imaginary, -real);
    }

    /**
     * Change the sign of the magnitude based on the signed value.
     *
     * <p>If the signed value is negative then the result is {@code -magnitude}; otherwise
     * return {@code magnitude}.
     *
     * <p>A signed value of {@code -0.0} is treated as negative. A signed value of {@code NaN}
     * is treated as positive.
     *
     * <p>This is not the same as {@link Math#copySign(double, double)} as this method
     * will change the sign based on the signed value rather than copy the sign.
     *
     * @param magnitude the magnitude
     * @param signedValue the signed value
     * @return magnitude or -magnitude.
     * @see #negative(double)
     */
    private static double changeSign(double magnitude, double signedValue) {
        return negative(signedValue) ? -magnitude : magnitude;
    }

    /**
     * Returns a scale suitable for use with {@link Math#scalb(double, int)} to normalise
     * the number to the interval {@code [1, 2)}.
     *
     * <p>The scale is typically the largest unbiased exponent used in the representation of the
     * two numbers. In contrast to {@link Math#getExponent(double)} this handles
     * sub-normal numbers by computing the number of leading zeros in the mantissa
     * and shifting the unbiased exponent. The result is that for all finite, non-zero,
     * numbers {@code a, b}, the magnitude of {@code scalb(x, -getScale(x))} is
     * always in the range {@code [1, 2)}, where {@code x = max(|a|, |b|)}.
     *
     * <p>This method is a functional equivalent of the c function ilogb(double) adapted for
     * two input arguments.
     *
     * <p>The result is to be used to scale a complex number using {@link Math#scalb(double, int)}.
     * Hence the special case of both zero arguments is handled using the return value for NaN
     * as zero cannot be scaled. This is different from {@link Math#getExponent(double)}
     * or {@link #getMaxExponent(double, double)}.
     *
     * <p>Special cases:
     *
     * <ul>
     * <li>If either argument is NaN or infinite, then the result is
     * {@link Double#MAX_EXPONENT} + 1.
     * <li>If both arguments are zero, then the result is
     * {@link Double#MAX_EXPONENT} + 1.
     * </ul>
     *
     * @param a the first value
     * @param b the second value
     * @return The maximum unbiased exponent of the values to be used for scaling
     * @see Math#getExponent(double)
     * @see Math#scalb(double, int)
     * @see <a href="http://www.cplusplus.com/reference/cmath/ilogb/">ilogb</a>
     */
    private static int getScale(double a, double b) {
        // Only interested in the exponent and mantissa so remove the sign bit
        final long x = Double.doubleToRawLongBits(a) & UNSIGN_MASK;
        final long y = Double.doubleToRawLongBits(b) & UNSIGN_MASK;
        // Only interested in the maximum
        final long max = Math.max(x, y);
        // Get the unbiased exponent
        int exp = (int) ((max >>> 52) - EXPONENT_OFFSET);

        // Do not distinguish nan/inf
        if (exp == Double.MAX_EXPONENT + 1) {
            return exp;
        }
        // Handle sub-normal numbers
        if (exp == Double.MIN_EXPONENT - 1) {
            // Special case for zero, return as nan/inf to indicate scaling is not possible
            if (max == 0) {
                return Double.MAX_EXPONENT + 1;
            }
            // A sub-normal number has an exponent below -1022. The amount below
            // is defined by the number of shifts of the most significant bit in
            // the mantissa that is required to get a 1 at position 53 (i.e. as
            // if it were a normal number with assumed leading bit)
            final long mantissa = max & MANTISSA_MASK;
            exp -= Long.numberOfLeadingZeros(mantissa << 12);
        }
        return exp;
    }

    /**
     * Returns the largest unbiased exponent used in the representation of the
     * two numbers. Special cases:
     *
     * <ul>
     * <li>If either argument is NaN or infinite, then the result is
     * {@link Double#MAX_EXPONENT} + 1.
     * <li>If both arguments are zero or subnormal, then the result is
     * {@link Double#MIN_EXPONENT} -1.
     * </ul>
     *
     * <p>This is used by {@link #divide(double, double, double, double)} as
     * a simple detection that a number may overflow if multiplied
     * by a value in the interval [1, 2).
     *
     * @param a the first value
     * @param b the second value
     * @return The maximum unbiased exponent of the values.
     * @see Math#getExponent(double)
     * @see #divide(double, double, double, double)
     */
    private static int getMaxExponent(double a, double b) {
        // This could return:
        // Math.getExponent(Math.max(Math.abs(a), Math.abs(b)))
        // A speed test is required to determine performance.
        return Math.max(Math.getExponent(a), Math.getExponent(b));
    }

    /**
     * Checks if both x and y are in the region defined by the minimum and maximum.
     *
     * @param x x value.
     * @param y y value.
     * @param min the minimum (exclusive).
     * @param max the maximum (exclusive).
     * @return true if inside the region.
     */
    private static boolean inRegion(double x, double y, double min, double max) {
        return (x < max) && (x > min) && (y < max) && (y > min);
    }

    /**
     * Creates an exception.
     *
     * @param message Message prefix.
     * @param error Input that caused the error.
     * @param cause Underlying exception (if any).
     * @return A new instance.
     */
    private static NumberFormatException parsingException(String message,
                                                          Object error,
                                                          Throwable cause) {
        // Not called with a null message or error
        final StringBuilder sb = new StringBuilder(100)
            .append(message)
            .append(" '").append(error).append('\'');
        if (cause != null) {
            sb.append(": ").append(cause.getMessage());
        }

        return new NumberFormatException(sb.toString());
    }

    /**
     * Returns {@code sqrt(x^2 + y^2)} without intermediate overflow or underflow.
     *
     * <p>Special cases:
     * <ul>
     * <li>If either argument is infinite, then the result is positive infinity.
     * <li>If either argument is NaN and neither argument is infinite, then the result is NaN.
     * </ul>
     *
     * <p>The computed result is expected to be within 1 ulp of the exact result.
     *
     * <p>This method is a replacement for {@link Math#hypot(double, double)}. In some cases
     * there may be differences between this method and {@code Math.hypot(double, double)} due
     * to the use of a different split algorithm to compute the high precision sum of
     * {@code x^2 + y^2}. Differences are expected to be 1 ULP indicating a rounding
     * change in the sum.
     *
     * <p>JDK9 ported the hypot function to Java for bug JDK-7130085 due to the slow performance
     * of the method as a native function. Benchmarks of the Complex class for functions that
     * use hypot confirm this is slow pre-Java 9. This implementation outperforms
     * {@code Math.hypot(double, double)} on JDK 9 and 11. See the Commons numbers examples JMH
     * module for benchmarks. Differences to Math.hypot have been observed at 1 ULP. Note
     * that the standard precision result {@code Math.sqrt(x * x + y * y)} is typically the
     * same result or within a few ULP. Comparisons with alternative implementations indicate
     * performance gains are related to edge case handling. This implementation handles most
     * infinite/finite and nan/finite combinations without testing for inf or nan as the test
     * on the exponent covers these cases.
     *
     * <p>This port was adapted from the c source code for the
     * "Freely Distributable Math Library" hypot function.
     * The original notice is shown below:
     * <pre>
     * ====================================================
     * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
     *
     * Developed at SunSoft, a Sun Microsystems, Inc. business.
     * Permission to use, copy, modify, and distribute this
     * software is freely granted, provided that this notice
     * is preserved.
     * ====================================================
     * </pre>
     *
     * <p>Note: The fdlibm c code makes use of the language ability to read and write directly
     * to the upper and lower 32-bits of the 64-double. The function performs
     * checking on the upper 32-bits for the magnitude of the two numbers by accessing
     * the exponent and 20 most significant bits of the mantissa. These upper bits
     * are manipulated during scaling and then used to perform extended precision
     * computation of the sum {@code x^2 + y^2} where the high part of the number has 20-bit
     * precision. Manipulation of direct bits has no equivalent in Java
     * other than use of {@link Double#doubleToLongBits(double)} and
     * {@link Double#longBitsToDouble(long)}. To avoid conversion to and from long and double
     * representations this implementation only scales the double representation. The high
     * and low parts of a double for the extended precision computation are extracted
     * using the method of Dekker (1971) to create two 26-bit numbers.
     *
     * <p>The computation of the high precision sum {@code x^2 + y^2} has been moved to a
     * separate method for use in different complex number algorithms that require the
     * absolute or absolute square. The algorithms thus function as if they all use the same
     * absolute value of the complex number as defined by this hypot(x, y) method.
     *
     * @param x Value x
     * @param y Value y
     * @return sqrt(x^2 + y^2)
     * @see Math#hypot(double, double)
     * @see <a href="https://www.netlib.org/fdlibm/e_hypot.c">fdlibm e_hypot.c</a>
     * @see <a href="https://bugs.java.com/bugdatabase/view_bug.do?bug_id=7130085">JDK-7130085 : Port fdlibm hypot to Java</a>
     */
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
        // The original computed scaling by directly writing to the exponent bits.
        // The high part was maintained during scaling for use in the high
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
        // a binary float value of 1.0 to create powers of 2, e.g. 0x1.0p+600 is 2^600.

        /* High word of x & y */
        // The mask is used to remove the sign.
        int ha = ((int) (Double.doubleToRawLongBits(x) >>> 32)) & UNSIGN_INT_MASK;
        int hb = ((int) (Double.doubleToRawLongBits(y) >>> 32)) & UNSIGN_INT_MASK;

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

        // Check if the smaller part is significant.
        // Do not replace this with 27 since the product x^2 is computed in
        // extended precision for an effective mantissa of 105-bits. Potentially it could be
        // replaced with 54 where y^2 will not overlap extended precision x^2.
        if ((ha - hb) > EXP_60) {
            /* x/y > 2**60 */
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        /* a <- |a| */
        /* b <- |b| */
        // No equivalent to directly writing back the high bits.
        // Just use Math.abs(). It is a hotspot intrinsic in Java 8+.
        a = Math.abs(a);
        b = Math.abs(b);

        // Second re-order in the rare event the upper 32-bits are the same.
        // This could be moved into the sub-normals section where it is critical.
        // It is left here to ensure a > b in all cases.
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
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return b == Double.POSITIVE_INFINITY ? b : a;
            }
            /* scale a and b by 2^-600 */
            a *= 0x1.0p-600;
            b *= 0x1.0p-600;
            rescale = 0x1.0p+600;
        } else if (hb < EXP_NEG_500) {
            /* b < 2^-500 */
            if (hb < EXP_NEG_1022) {
                /* sub-normal or 0 */
                // Intentional comparison with zero.
                if (b == 0.0) {
                    return a;
                }
                /* scale a and b by 2^1022 */
                a *= 0x1.0p+1022;
                b *= 0x1.0p+1022;
                rescale = 0x1.0p-1022;
            } else {
                /* scale a and b by 2^600 */
                a *= 0x1.0p+600;
                b *= 0x1.0p+600;
                rescale = 0x1.0p-600;
            }
        }

        return Math.sqrt(x2y2(a, b)) * rescale;
    }

    /**
     * Return {@code x^2 + y^2} with high accuracy.
     *
     * <p>It is assumed that {@code 2^500 > x >= y > 2^-500}. Thus there will be no
     * overflow or underflow of the result.
     *
     * <p>Translated from "Freely Distributable Math Library" e_hypot.c
     * to minimize rounding errors. The code has been modified to use the split
     * method from Dekker. See the method {@link #hypot(double, double)} for the
     * original license.
     *
     * <p>Note that fdlibm conditioned {@code x > y} using only the upper 20 bits of the
     * mantissa and the 11-bit exponent and thus x could be slightly smaller than y
     * for cases where the upper 31-bits of the unsigned 64-bit representation are identical.
     * There are cases for {@code 2x > y > x} where the result is less exact than if the
     * arguments were reversed so the equality should be conditioned on the entire magnitude
     * of the two floating point values. This is especially relevant if the numbers are
     * sub-normal where the upper 20-bits of the mantissa are zero and the entire information
     * on the magnitude is in the lower 32-bits.
     *
     * @param x Value x.
     * @param y Value y
     * @return x^2 + y^2
     */
    private static double x2y2(double x, double y) {
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
            final double t1 = splitHigh(x);
            final double t2 = x - t1;
            return t1 * t1 - (y * (-y) - t2 * (x + t1));
        }
        // 2y > x > y
        final double t = x + x;
        final double y1 = splitHigh(y);
        final double y2 = y - y1;
        final double t1 = splitHigh(t);
        final double t2 = t - t1;
        return t1 * y1 - (w * (-w) - (t1 * y2 + t2 * y));
    }
}
