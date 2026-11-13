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
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.SplittableRandom;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.BiFunction;
import java.util.function.BooleanSupplier;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoublePredicate;
import java.util.function.DoubleSupplier;
import java.util.function.DoubleUnaryOperator;
import java.util.function.IntSupplier;
import java.util.stream.Stream;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.numbers.core.DD;
import org.apache.commons.numbers.core.DDMath;
import org.apache.commons.numbers.fraction.BigFraction;
import org.apache.commons.numbers.rootfinder.BrentSolver;
import org.apache.commons.rng.sampling.distribution.AliasMethodDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateDiscreteSampler;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Test the {@link HurwitzZeta} function.
 *
 * <p>This class contains a test implementation for the Euler-Maclaurin formula which can
 * be configured using an evaluation context. This changes the number of terms in the
 * series summation (N); the power function used for the terms; and methods to compute the
 * tail series. Test fixtures exist to run reference data for different contexts to output
 * the precision. These are disabled but can be used to benchmark any changes to the zeta
 * implementation. The final implementation has been chosen using options that achieve a
 * consistent low maximum and root mean squared error across various test data.
 *
 * <p>Test fixtures used to generate the data are also provided. The evaluation of
 * zeta(s, a) on the test data was performed using variable precision implementations in
 * MATLAB 2026a; or Python's mpmath 1.4.1. Note that mpmath was preferred as the software is
 * freely available. In some cases mpmath does not converge (it computes different results
 * in the first 15 digits when the number of digits of precision is changed between 30 and
 * 100). When the java implementation failed on test cases the data was either recomputed
 * with higher digits of precision; the reference implementation changed; and/or cases
 * cross verified using WolframAlpha HurwitzZeta online evaluator and cases updated to the
 * consensus.
 *
 * <p>The Euler-Maclaurin formula integrals are valid for (a+N) above zero, and evaluation
 * is more stable when (a+N) is many integers above 0. Evaluation of the zeta function
 * with negative argument {@code a} requires summation of (a+k)^-s until (a+k) is
 * positive. This sum can be added to the zeta evaluation for the positive (a+k). When the
 * power term is negative the sum of terms for (a+k) below zero is negative and
 * cancellation occurs. This class contains a BigDecimal and double-double (DD) test
 * implementation that sums the terms in extended precision to reduce cancellation. This
 * is possible as the power term s for negative a must be an integer, otherwise the zeta
 * result is complex which this implementation does not support. Both BigDecimal and DD
 * support the power function with an integer exponent, and can compute values that
 * overflow or underflow a double.
 *
 * <p>Cases of extreme cancellation occurs at the roots of the function. A test fixture
 * is provided that searches for roots and saves the (s, a) cases for reference evaluation.
 * These cases require a double-double implementation to evaluate. Using the double implementation
 * the cancellation can be so large that the result may have no correct bits.
 *
 * <p>Relation to polygamma
 *
 * <p>The Hurwitz zeta function is related to the polygamma function (psi) for positive integer
 * n != 0 (which is the digamma function). A reflection formula exists for polygamma:
 *
 * <pre>
 * psi(n, x) = (-1)^(n+1) n! zeta(n+1, x)
 *
 *                                           d^n cot(pi x)
 * psi(n, 1-x) = (-1)^(n) psi(n, x) + (-1)^n -------------
 *                                              d x^n
 * </pre>
 *
 * <p>Where [ d^n cot(pi x) / d x^n ] is the n-th partial derivative of cot(pi x)
 * (cos(pi x) / sin(pi x)).
 * This derivative can be tabulated for small n and has a general form for large n (see
 * <a href="https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/sf_gamma/polygamma.html">
 * Boost Polygamma documentation</a>). No attempt was made to use a polygamma reflection
 * for negative arguments. The derivative of the cot function would be increasingly
 * imprecise with large n when implemented using double precision trigonometric functions.
 */
class HurwitzZetaTest {
    /** Seed for data generation. */
    private static final long SEED = 6516587940839803692L;
    /** Table used to create a histogram of the number of steps to converge the sum series.
     * Used in the test zeta implementation. */
    private static final int[] N = new int[107];
    /** Table used to create a histogram of the number of steps to converge the tail series.
     * Used in the test zeta implementation. */
    private static final int[] M = new int[107];
    /** Retain the number of terms used from the tail series in the test implementation. */
    private static int lastM = 0;
    /** Filenames of resources used for the test zeta function. */
    private static final String[] TEST_RESOURCES = {
        "hzeta_s1_4_a1_8.csv",
        "hzeta_s1_4_a8_32.csv",
        "hzeta_s1_4_a32_2147483648.csv",
        "hzeta_s1_4_a0_1.csv",
        "hzeta_s1_4_a1e-16_1e-14.csv",
        "hzeta_s4_32_a1_8.csv",
        "hzeta_s1_19_a2147483648_10000000000000000.csv",
        "hzeta_domain_2p-1022.csv",
        "hzeta_domain_2p-1075.csv",
    };
    /** Filenames of resources used for the extended precision zeta function using integer s
     * and a in the expected range for the negative a use case. This requires {@code 0 < a < 1}
     * and {@code a ~ 2N}. Uses small s. When s is large the power series is strongly converging
     * for {@code 0 < a < 1}, and the term very small for {@code a ~ 2N} so it does not impact
     * the final result. */
    private static final String[] INT_TEST_RESOURCES = {
        // a -> 1 so that the initial term approaches the smallest magnitude (1.0), and is
        // the closest magnitude to remaining terms below 1.0.
        "hzeta_is2_11_a1_1.csv",
        // a ~ 2N
        "hzeta_is2_11_a40_41.csv",
    };
    /** Filenames of resources used for the extended precision zeta function using integer s.
     * Includes {@link #INT_TEST_RESOURCES} and contains data where s ~ a to test convergence
     * of the tail series in extended precision. */
    private static final String[] EXT_INT_TEST_RESOURCES = {
        "hzeta_is2_11_a1_1.csv",
        "hzeta_is2_11_a40_41.csv",
        // Data with s ~ a to test tail convergence
        "hzeta_s9_1075_tail.csv",
        "hzeta_s9_1837_tail.csv",
        "hzeta_s15_1128_tail.csv",
        "hzeta_is_2p-1075.csv",
    };
    /** Filenames of resources used for the roots of the zeta function using integer s. */
    private static final String[] ROOT_TEST_RESOURCES = {
        "hzeta_root_s3_21_na2b0_2b6.csv",
        "hzeta_root_s23_1067_na2b0_2b6.csv",
        "hzeta_root_s3_25_na2b7_2b11.csv",
        "hzeta_root_s3_21_na2b12_2b14.csv",
        "hzeta_root_s3_3_na1001_16384.csv",
    };
    /** Filenames of resources used for zeta function using sub-normal terms. */
    private static final String[] SUB_NORMAL_RESOURCES = {
        "hzeta_s20_1024_ap9_subnormal1023.csv",
        "hzeta_s20_1024_ap1_subnormal1040.csv",
        "hzeta_s20_1024_ap1_subnormal1053.csv",
        "hzeta_s20_1024_ap0.1_subnormal1073.csv",
        // This tests the tail is used to at least 7 terms when N=7
        "hzeta_2p-1022.csv",
        // This tests scaling of terms S and I. The tail is not used
        // as the result values only have 10-bits.
        "hzeta_2p-1065.csv",
    };
    /** Filenames of resources used for zeta function using the regions of the
     * domain (s, a). */
    private static final String[] DOMAIN_RESOURCES = {
        "hzeta_domain_2p-1022.csv",
        "hzeta_domain_2p-1075.csv",
    };

    // For the configurable double implementation these constants are from the main class

    /** Asymptotic threshold for large {@code a}. Used when {@code a + N} is not accurate.
     * Equal to 2^53. Configurable N should be subtracted from this value. */
    private static final double LARGE_APN = 1L << 53;
    /** Threshold used when {@code a < 1} where {@code a^-s} == zeta(s, a) with the smallest supported
     * {@code s = 1 + 2^-52}. Set using zeta(s, 1) * 2^54. No possible evaluation of
     * all remaining terms for any s could be added to {@code a^-s} when {@code a < 1}.
     *
     * <p>Approximately equal to 2^106 thus remaining terms (all less than 1) cannot be added in
     * double-double precision. */
    private static final double LARGE_T0 = 8.11296384146066920939820185959e+31;
    /** Downscale to use on the initial rising factorial tail term (always above 1).
     * Equal to 2^-1022 * 2^106 = 2^-916. */
    private static final double DOWNSCALE_916 = 0x1.0p-916;
    /** Downscale to use on the final sum.
     * Equal to 2^-106 */
    private static final double DOWNSCALE_106 = 0x1.0p-106;
    /** Upscale to use on the initial probability tail term (always below 1).
     * Multiplicative inverse of Double.MIN_NORMAL: 2^1022. */
    private static final double UPSCALE_1022 = 0x1.0p1022;
    /** Upscale to use on asymptotic approximation. Equal to 2^916. */
    private static final double UPSCALE_916 = 0x1.0p916;
    /** Upscale to use on the sum before scaled tail computation. Allows convergence
     * on a normal number using a sum of normal terms.
     * Multiplicative inverse of DOWNSCALE2: 2^106. */
    private static final double UPSCALE_106 = 0x1.0p106;
    /** Threshold for z in the Taylor series approximation of (1 + t)^z.
     * Set where |z * t| < 2^-26  =>  z < 2^-26 / 2^-52 = 2^26. */
    private static final double Z_LIMIT = 0x1p26;

    /**
     * Flag set when reporting to the console. Used for testing.
     * If negative no output is printed.
     *
     * <p>Some test fixtures will exercise the testing code. When no reporting
     * is configured then the test will still run but with a subset of the
     * cases. This ensures test code is asserted but without excessive runtimes.
     */
    private static int reporting = 0;

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses {@code 2k!} divided by Bernoulli number {@code B_2k}.
     * Precomputed for use in test implementations.
     * Note: Prefix R is for reciprocal and this must be used with division.
     */
    private static final double[] RF = {
        12.0,
        -720.0,
        30240.0,
        -1209600.0,
        4.790016E7,
        -1.8924375803183792E9,
        7.47242496E10,
        -2.950130727918164E12,
        1.1646782814350067E14,
        -4.597978722407473E15,
        1.81521054019435456E17,
        -7.1661652561756672E18,
        2.82908877253043E20,
        -1.1168794925000445E22,
        4.4092635141854666E23,
        -1.7407074646225822E25,
        6.872037622739274E26,
        -2.7129717107520044E28,
        1.0710383014704457E30,
        -4.228289733582729E31,
        1.669261878547101E33,
        -6.589981753232787E34,
        2.6016205165923062E36,
        -1.0270786120209628E38,
        4.054743835786749E39,
        -1.6007487042788347E41,
        6.319502582715391E42,
        -2.4948396201225357E44,
        9.849232037909393E45,
        -3.888320954748035E47,
        1.535047584313167E49,
        -6.060124957607529E50,
        2.3924414381101893E52,
        -9.444980218768351E53,
        3.728728733414322E55,
        -1.4720431007109737E57,
        5.811393226148101E58,
        -2.2942460864500874E60,
        9.057320508803927E61,
        -3.575686814230726E63,
        1.4116245727459505E65,
        -5.5728704383437275E66,
        2.2000810641991213E68,
        -8.685571901589203E69,
        3.428926346636115E71,
        -1.3536858624708422E73,
        5.344137578373867E74,
        -2.10978095054183E76,
        8.329081341920855E77,
        -3.288189514770133E79,
        1.2981251882636475E81,
        -5.124792828500739E82,
        2.023187114193683E84,
        -7.9872225785893E85,
        3.1532290845650105E87,
        -1.2448449460266318E89,
        4.914450863191316E90,
        -1.9401474347316165E92,
        7.659395064235756E93,
        -3.023807969426533E95,
        1.1937515377240477E97,
        -4.71274217221138E98,
        1.8605160353622745E100,
        -7.345022900363532E101,
        2.899698813741201E103,
        -1.14475520695735E105,
        4.519312411502489E106,
        -1.7841530266585092E108,
        7.0435538256502885E109,
        -2.780683593477915E111,
        1.097769881289064E113,
        -4.333821780701556E114,
        1.71092426081396E116,
        -6.754458245784006E117,
        2.6665532331825656E119,
        -1.0527130210383084E121,
        4.155944426209505E122,
        -1.640701095984005E124,
        6.4772283031183474E125,
        -2.557107238692696E127,
        1.0095054742823555E129,
        -3.985367868760375E130,
        1.5733601702991004E132,
        -6.211376984513082E133,
        2.4521533449270183E135,
        -9.680713378095076E136,
        3.8217924544852726E138,
        -1.5087831851535177E140,
        5.956437265792312E141,
        -2.351507178131059E143,
        9.283378237790206E144,
        -3.664930828506856E146,
        1.4468566973887733E148,
        -5.7119612913175405E149,
        2.2549919319855863E151,
        -8.902351318538374E152,
        3.5145074301396E154,
        -1.387471920006682E156,
        5.47751958723434E157,
        -2.1624380570088485E159,
        8.53696326581506E160,
        -3.37025800880906E162,
        1.3305245310619423E164,
        -5.252700307010521E165,
        2.0736829627069785E167,
        -8.186572198078725E168,
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * Precomputed for use in test implementations.
     * This is the inverse of table {@link #M}.
     */
    private static final double[] F = {
        0.08333333333333333,
        -0.001388888888888889,
        3.306878306878307E-5,
        -8.267195767195768E-7,
        2.08767569878681E-8,
        -5.284190138687493E-10,
        1.3382536530684679E-11,
        -3.3896802963225827E-13,
        8.586062056277845E-15,
        -2.174868698558062E-16,
        5.5090028283602295E-18,
        -1.3954464685812522E-19,
        3.534707039629467E-21,
        -8.953517427037546E-23,
        2.267952452337683E-24,
        -5.744790668872202E-26,
        1.455172475614865E-27,
        -3.6859949406653103E-29,
        9.336734257095045E-31,
        -2.36502241570063E-32,
        5.990671762482134E-34,
        -1.5174548844682903E-35,
        3.843758125454189E-37,
        -9.736353072646691E-39,
        2.466247044200681E-40,
        -6.247076741820743E-42,
        1.5824030244644914E-43,
        -4.008273685948936E-45,
        1.0153075855569557E-46,
        -2.5718041582418717E-48,
        6.514456035233815E-50,
        -1.6501309906896525E-51,
        4.179830628539476E-53,
        -1.058763466770291E-54,
        2.6818791912607708E-56,
        -6.793279351107421E-58,
        1.7207577616681404E-59,
        -4.358730329348894E-61,
        1.1040792903684666E-62,
        -2.7966655133781345E-64,
        7.084036501679471E-66,
        -1.794407408289224E-67,
        4.545287063611096E-69,
        -1.1513346631982051E-70,
        2.9163647710923614E-72,
        -7.387238263497337E-74,
        1.8712093117637953E-75,
        -4.739828557761799E-77,
        1.2006125993354507E-78,
        -3.0411872415142924E-80,
        7.703417274705106E-82,
        -1.951298390909883E-83,
        4.942696565159462E-85,
        -1.2519996659171848E-86,
        3.1713522017635153E-88,
        -8.033128970735334E-90,
        2.0348153391661465E-91,
        -5.154247466447474E-93,
        1.3055861352149468E-94,
        -3.307088314175091E-96,
        8.376952560049092E-98,
        -2.1219068717497136E-99,
        5.37485289561228E-101,
        -1.3614661432172069E-102,
        3.448634027993399E-104,
        -8.735492041638356E-106,
        2.2127259833925496E-107,
        -5.604900392837224E-109,
        1.419737854999179E-110,
        -3.596237998258763E-112,
        9.109377266078232E-114,
        -2.307432217109123E-115,
        5.844794085299002E-117,
        -1.4805036371705745E-118,
        3.7501595226227194E-120,
        -9.499265041992958E-122,
        2.40619194446752E-123,
        -6.094955397102685E-125,
        1.5438702377042473E-126,
        -3.9106689968592925E-128,
        9.90584028987943E-130,
        -2.5091786578563554E-131,
        6.35582378960246E-133,
        -1.609948973461625E-134,
        4.0780483898724625E-136,
        -1.032981724531519E-137,
        2.6165732752609202E-139,
        -6.627857533408623E-141,
        1.6788559257443675E-142,
        -4.2525917390343004E-144,
        1.0771940713664573E-145,
        -2.728564458083958E-147,
        6.911534513437014E-149,
        -1.7507121442157682E-150,
        4.43460566672392E-152,
        -1.123298737848715E-153,
        2.8453489425693972E-155,
        -7.207353068415137E-157,
        1.8256438595501416E-158,
        -4.624409918974655E-160,
        1.1713767165946985E-161,
        -2.9671318854112526E-163,
        7.515832866319735E-165,
        -1.9037827051837493E-166,
        4.822337927175732E-168,
        -1.221512466761957E-169,
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * Precomputed for use in test implementations.
     * This is the same as table {@link #F} in 70-digits of precision.
     */
    private static final BigDecimal[] FBD = {
        new BigDecimal("0.08333333333333333333333333333333333333333333333333333333333333333333333"),
        new BigDecimal("-0.001388888888888888888888888888888888888888888888888888888888888888888889"),
        new BigDecimal("0.00003306878306878306878306878306878306878306878306878306878306878306878307"),
        new BigDecimal("-8.267195767195767195767195767195767195767195767195767195767195767195767E-7"),
        new BigDecimal("2.087675698786809897921009032120143231254342365453476564587675698786810E-8"),
        new BigDecimal("-5.284190138687493184847682202179556676911174265671620168974666329163684E-10"),
        new BigDecimal("1.338253653068467883282698097512912327727142541957356772171586986401801E-11"),
        new BigDecimal("-3.389680296322582866830195391249442499572181074411596249961225581103162E-13"),
        new BigDecimal("8.586062056277844564135905450425627133953956127486574969878449381578335E-15"),
        new BigDecimal("-2.174868698558061873041516423865917899851791600154166083524811874481630E-16"),
        new BigDecimal("5.509002828360229515202652608902254877861582704468840789377736849075333E-18"),
        new BigDecimal("-1.395446468581252334070768626406354976391763666911099523783965934267938E-19"),
        new BigDecimal("3.534707039629467471693229977803799214724594564714920359506997731309122E-21"),
        new BigDecimal("-8.953517427037546850402611318112741051627139242784962516444340865535696E-23"),
        new BigDecimal("2.267952452337683060310950738868166063220354329744171261064703326867830E-24"),
        new BigDecimal("-5.744790668872202445263881987607018399624776691789095171764503043537801E-26"),
        new BigDecimal("1.455172475614864901866264867271329335720888955829047351288633402121478E-27"),
        new BigDecimal("-3.685994940665310178181782479908660374446298206435960844182904657554629E-29"),
        new BigDecimal("9.336734257095044672032555152785623295443688712038371425860508861190119E-31"),
        new BigDecimal("-2.365022415700629934559635196369838240069656250349772795125530143735169E-32"),
        new BigDecimal("5.990671762482134304659912396819657826449369990395618132966642413079763E-34"),
        new BigDecimal("-1.517454884468290261710813135864718931540884301245442686973614561814395E-35"),
        new BigDecimal("3.843758125454188232229445290990232105901809047372886539813079625237514E-37"),
        new BigDecimal("-9.736353072646691035267621279250454180955109079541241688773435358308703E-39"),
        new BigDecimal("2.466247044200680957106400280288842885924177338404622199793367088926577E-40"),
        new BigDecimal("-6.247076741820743693148756794723368692576573980221262151470368932374574E-42"),
        new BigDecimal("1.582403024464491429751081706828763940328602762417916616214025303901512E-43"),
        new BigDecimal("-4.008273685948935968530012190521982662681125480849604027743896411659203E-45"),
        new BigDecimal("1.015307585556955631163071394537876232706779463883299635682782384116126E-46"),
        new BigDecimal("-2.571804158241871749924819409764454885573157750221760160818715882306806E-48"),
        new BigDecimal("6.514456035233814931558434858641858023142039648061506239697847672837607E-50"),
        new BigDecimal("-1.650130990689652455506098780479323009187918309124390360514110851564175E-51"),
        new BigDecimal("4.179830628539475894850187234709407032931287099533911224921631428065113E-53"),
        new BigDecimal("-1.058763466770290877027042024279117287335207124617836065073387616228668E-54"),
        new BigDecimal("2.681879191260770666140984858841510339769094268609011851819904574290494E-56"),
        new BigDecimal("-6.793279351107421209527180299533894611894542453231341127791219156100086E-58"),
        new BigDecimal("1.720757761668140490536349940758230664281587517995397800783755191115580E-59"),
        new BigDecimal("-4.358730329348893843400199849773161109127526522657475973184850396137227E-61"),
        new BigDecimal("1.104079290368466675083839597644427323087385806923686190398281791181660E-62"),
        new BigDecimal("-2.796665513378134507204793753118626553864489833892558755670206138679735E-64"),
        new BigDecimal("7.084036501679470198509388422380349334323387828252242164209997075994523E-66"),
        new BigDecimal("-1.794407408289224066605257309336755986704921973905220023359133870266054E-67"),
        new BigDecimal("4.545287063611096107085079124639448185863419826897374715989665643907358E-69"),
        new BigDecimal("-1.151334663198205181273002900786243377460268013686045202405101967894282E-70"),
        new BigDecimal("2.916364771092361354703368980052970742786734157547911896568980746396497E-72"),
        new BigDecimal("-7.387238263497337562573375394709736697786504856753883327773382790549657E-74"),
        new BigDecimal("1.871209311763795306225261871408074659301018846222338710574930531250668E-75"),
        new BigDecimal("-4.739828557761799405499563441210887079197749225080190616793274013578142E-77"),
        new BigDecimal("1.200612599335450651981710022114129123651624711698687507489471047133328E-78"),
        new BigDecimal("-3.041187241514292383037122074838895532798552557535210001089401413468719E-80"),
        new BigDecimal("7.703417274705106272876503394597218882452448065616497738522362835029019E-82"),
        new BigDecimal("-1.951298390909883071112323768145617747486777932628731259552205113296608E-83"),
        new BigDecimal("4.942696565159461474896400153888327970275086487601097361422020921854949E-85"),
        new BigDecimal("-1.251999665917184792168185768825538225661310147841267501859258817882884E-86"),
        new BigDecimal("3.171352201763515460645395167591419654594408306209802797871271530980165E-88"),
        new BigDecimal("-8.033128970735334461394957124140604320908863896357851007376393644320001E-90"),
        new BigDecimal("2.034815339166146570781512075871595570547503242380639808183083362513841E-91"),
        new BigDecimal("-5.154247466447473859069509136779375774820347553341655041460275446975888E-93"),
        new BigDecimal("1.305586135214946724576575477100187893779561015049220316951635323254731E-94"),
        new BigDecimal("-3.307088314175091248451753981253057414970484831126237049698444217345473E-96"),
        new BigDecimal("8.376952560049091303035075843520726698500871006499291574198603410342316E-98"),
        new BigDecimal("-2.121906871749713769528898847038058622493036075204024715103185826773423E-99"),
        new BigDecimal("5.374852895612280256280556992243202193671362471266005683176269739853926E-101"),
        new BigDecimal("-1.361466143217206939250409531337985056551186635590368193154033073280195E-102"),
        new BigDecimal("3.448634027993399034277374049573816904552346268902603942292514314919319E-104"),
        new BigDecimal("-8.735492041638355060271261209502790508674540352967083633115097064950438E-106"),
        new BigDecimal("2.212725983392549708431115456213317404319773840076968700330947158249593E-107"),
        new BigDecimal("-5.604900392837224170822067787299057694158160355283377216951480870287638E-109"),
        new BigDecimal("1.419737854999178767309682770362704493752286624636918355949881058878973E-110"),
        new BigDecimal("-3.596237998258762663674623971191536517120563267758289607457687703620313E-112"),
        new BigDecimal("9.109377266078231864576892767139751023249607479194325282636316738898584E-114"),
        new BigDecimal("-2.307432217109123288503643570464850300232976414806213340886269469255711E-115"),
        new BigDecimal("5.844794085299001994493088929867304932408032086274505089273608956599896E-117"),
        new BigDecimal("-1.480503637170574495251820316298791125081247803550335046176209349324214E-118"),
        new BigDecimal("3.750159522622719689081176129191400562349688962715649262138170984615923E-120"),
        new BigDecimal("-9.499265041992958297390439583396763016764691865952385722782710815110294E-122"),
        new BigDecimal("2.406191944467519874474430124120706183373632626965386180864448667151536E-123"),
        new BigDecimal("-6.094955397102684750831061254616638145129437344495943436935650807455717E-125"),
        new BigDecimal("1.543870237704247145637618175952270849803149674487976961421141600579795E-126"),
        new BigDecimal("-3.910668996859292320146383100735052562051738589560318147853630366021334E-128"),
        new BigDecimal("9.905840289879429767932893840368913666671560408193911745315523946847522E-130"),
        new BigDecimal("-2.509178657856355318831305287294069779076443074810090295827624189594996E-131"),
        new BigDecimal("6.355823789602459805358145464286720693075403061385930768584203812506772E-133"),
        new BigDecimal("-1.609948973461625039730184241507041411509020052499606137340285253441347E-134"),
        new BigDecimal("4.078048389872462246989648408768057411410825453086287743426045724368654E-136"),
        new BigDecimal("-1.032981724531519021447019174458280023312613407565355929241614077244873E-137"),
        new BigDecimal("2.616573275260920172946141886267130378865746840754578100904165069822076E-139"),
        new BigDecimal("-6.627857533408622801054493484333313814259947931062379736028850847212972E-141"),
        new BigDecimal("1.678855925744367305152093934168830190524066196163748280275319889626768E-142"),
        new BigDecimal("-4.252591739034300621563005751491573858481678065504069572450932695353594E-144"),
        new BigDecimal("1.077194071366457292665858859870733628695784553089312112918448661341321E-145"),
        new BigDecimal("-2.728564458083957965649286777541512367443093131510097703167206577666725E-147"),
        new BigDecimal("6.911534513437013691212132922676748830088259405954869110009975085805334E-149"),
        new BigDecimal("-1.750712144215768249498497911921168010135817605446930257283751951348322E-150"),
        new BigDecimal("4.434605666723919586207981310089151331192950428073856854356039822189180E-152"),
        new BigDecimal("-1.123298737848714970168262281651073691524745415149383192677583791635122E-153"),
        new BigDecimal("2.845348942569397093568104834095780921404892435257915257819150653602841E-155"),
        new BigDecimal("-7.207353068415136796682689400677207295096178887722574917666061808587523E-157"),
        new BigDecimal("1.825643859550141758304265123776272180908856385216002502109492295542019E-158"),
        new BigDecimal("-4.624409918974655474068228810977201137320753264948261557312628309361878E-160"),
        new BigDecimal("1.171376716594698503117925402612594894238068184862096636384363155836369E-161"),
        new BigDecimal("-2.967131885411252338073120134667304797306190472043207035461074138316972E-163"),
        new BigDecimal("7.515832866319735281505615656771271060714655242703701767051043339757299E-165"),
        new BigDecimal("-1.903782705183749422163676479489486484370791872044812155268219988500410E-166"),
        new BigDecimal("4.822337927175731599757997875780745381026933501142317572734300292543235E-168"),
        new BigDecimal("-1.221512466761956929673298475247707655443475938093769491269909895545982E-169"),
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * Precomputed for use in test implementations.
     * This is the same as table {@link #F} in double-double precision.
     */
    private static final DD[] FDD = {
        DD.ofSum(0.08333333333333333, 4.625929269271485E-18),
        DD.ofSum(-0.001388888888888889, 5.300543954373577E-20),
        DD.ofSum(3.306878306878307E-5, -2.2300719288557665E-21),
        DD.ofSum(-8.267195767195768E-7, 3.457597454003665E-23),
        DD.ofSum(2.08767569878681E-8, -1.2073450591132599E-24),
        DD.ofSum(-5.284190138687493E-10, 3.517096671929869E-27),
        DD.ofSum(1.3382536530684679E-11, -2.828354019907999E-29),
        DD.ofSum(-3.3896802963225827E-13, -1.4986928409964295E-29),
        DD.ofSum(8.586062056277845E-15, -6.05252374381974E-31),
        DD.ofSum(-2.174868698558062E-16, 4.961617782549996E-33),
        DD.ofSum(5.5090028283602295E-18, -1.49827152194499E-35),
        DD.ofSum(-1.3954464685812522E-19, -1.0350590497256251E-35),
        DD.ofSum(3.534707039629467E-21, 1.894231142684204E-37),
        DD.ofSum(-8.953517427037546E-23, -5.728752743153026E-39),
        DD.ofSum(2.267952452337683E-24, 1.3043458462619563E-40),
        DD.ofSum(-5.744790668872202E-26, 1.663242973708004E-43),
        DD.ofSum(1.455172475614865E-27, -5.613265715443096E-44),
        DD.ofSum(-3.6859949406653103E-29, 1.0778256413554197E-45),
        DD.ofSum(9.336734257095045E-31, -3.9347970210731877E-47),
        DD.ofSum(-2.36502241570063E-32, 2.0347170931532494E-49),
        DD.ofSum(5.990671762482134E-34, 1.6265467158179092E-50),
        DD.ofSum(-1.5174548844682903E-35, 5.493014407946745E-52),
        DD.ofSum(3.843758125454189E-37, -3.685053096067968E-53),
        DD.ofSum(-9.736353072646691E-39, 2.258059165188444E-55),
        DD.ofSum(2.466247044200681E-40, -1.505641802268162E-56),
        DD.ofSum(-6.247076741820743E-42, -2.7106815859687654E-58),
        DD.ofSum(1.5824030244644914E-43, 2.545428531496969E-60),
        DD.ofSum(-4.008273685948936E-45, -2.2124211668946826E-61),
        DD.ofSum(1.0153075855569557E-46, -9.404269751258486E-63),
        DD.ofSum(-2.5718041582418717E-48, -6.537655454012542E-65),
        DD.ofSum(6.514456035233815E-50, -2.763626172529861E-66),
        DD.ofSum(-1.6501309906896525E-51, 3.1794529475063687E-68),
        DD.ofSum(4.179830628539476E-53, 2.617556823159939E-69),
        DD.ofSum(-1.058763466770291E-54, 6.6915528436035195E-71),
        DD.ofSum(2.6818791912607708E-56, -8.70695425146146E-73),
        DD.ofSum(-6.793279351107421E-58, 2.795667911354165E-74),
        DD.ofSum(1.7207577616681404E-59, 4.65433497191727E-76),
        DD.ofSum(-4.358730329348894E-61, 2.8840522874209336E-77),
        DD.ofSum(1.1040792903684666E-62, 6.624841731022409E-79),
        DD.ofSum(-2.7966655133781345E-64, 2.628041826403209E-81),
        DD.ofSum(7.084036501679471E-66, -5.026235239023924E-82),
        DD.ofSum(-1.794407408289224E-67, 1.5372719769275798E-84),
        DD.ofSum(4.545287063611096E-69, 9.87696151726261E-87),
        DD.ofSum(-1.1513346631982051E-70, -7.192856523313341E-87),
        DD.ofSum(2.9163647710923614E-72, -3.8911087510195904E-89),
        DD.ofSum(-7.387238263497337E-74, -6.923136687699924E-90),
        DD.ofSum(1.8712093117637953E-75, 1.5886680102062367E-92),
        DD.ofSum(-4.739828557761799E-77, -9.517121002177184E-94),
        DD.ofSum(1.2006125993354507E-78, -1.109850335891779E-95),
        DD.ofSum(-3.0411872415142924E-80, 5.125117133572647E-97),
        DD.ofSum(7.703417274705106E-82, 3.948211996024456E-99),
        DD.ofSum(-1.951298390909883E-83, -1.2533409416284754E-99),
        DD.ofSum(4.942696565159462E-85, -2.8094990080509668E-101),
        DD.ofSum(-1.2519996659171848E-86, -2.5859037541718075E-103),
        DD.ofSum(3.1713522017635153E-88, 1.808466797142321E-104),
        DD.ofSum(-8.033128970735334E-90, -5.501602254025203E-106),
        DD.ofSum(2.0348153391661465E-91, 4.956249311149982E-108),
        DD.ofSum(-5.154247466447474E-93, -2.4639060810936437E-109),
        DD.ofSum(1.3055861352149468E-94, -3.669056422783558E-111),
        DD.ofSum(-3.307088314175091E-96, -6.886640190175871E-113),
        DD.ofSum(8.376952560049092E-98, -6.082546927852001E-114),
        DD.ofSum(-2.1219068717497136E-99, -1.893178093634648E-115),
        DD.ofSum(5.37485289561228E-101, 1.5431494060360965E-118),
        DD.ofSum(-1.3614661432172069E-102, -4.0565948103494285E-119),
        DD.ofSum(3.448634027993399E-104, 2.9458047085233967E-120),
        DD.ofSum(-8.735492041638356E-106, 9.4189357054137E-122),
        DD.ofSum(2.2127259833925496E-107, 1.2061839538632616E-123),
        DD.ofSum(-5.604900392837224E-109, 1.1590480702591175E-125),
        DD.ofSum(1.419737854999179E-110, -1.0538661729107658E-126),
        DD.ofSum(-3.596237998258763E-112, 1.3118059698383575E-128),
        DD.ofSum(9.109377266078232E-114, -5.075775794901921E-131),
        DD.ofSum(-2.307432217109123E-115, -1.0179374029871243E-131),
        DD.ofSum(5.844794085299002E-117, 3.158983160261219E-133),
        DD.ofSum(-1.4805036371705745E-118, -1.8032521549667245E-135),
        DD.ofSum(3.7501595226227194E-120, 2.8585141457797E-136),
        DD.ofSum(-9.499265041992958E-122, -2.9137013725177297E-139),
        DD.ofSum(2.40619194446752E-123, -2.156821690478659E-140),
        DD.ofSum(-6.094955397102685E-125, -2.4131212075875234E-141),
        DD.ofSum(1.5438702377042473E-126, -1.4827854235963976E-142),
        DD.ofSum(-3.9106689968592925E-128, 1.9958189536928135E-144),
        DD.ofSum(9.90584028987943E-130, -5.990937108052624E-146),
        DD.ofSum(-2.5091786578563554E-131, 1.280123839457573E-147),
        DD.ofSum(6.35582378960246E-133, -3.3994016155804947E-149),
        DD.ofSum(-1.609948973461625E-134, 1.1085599508572997E-151),
        DD.ofSum(4.0780483898724625E-136, -2.2687099531588687E-152),
        DD.ofSum(-1.032981724531519E-137, -9.220601136512268E-155),
        DD.ofSum(2.6165732752609202E-139, -4.3351815619058493E-156),
        DD.ofSum(-6.627857533408623E-141, 5.518180594483479E-157),
        DD.ofSum(1.6788559257443675E-142, -1.7064658138350165E-158),
        DD.ofSum(-4.2525917390343004E-144, -2.3021245512754632E-160),
        DD.ofSum(1.0771940713664573E-145, -2.977662208872027E-162),
        DD.ofSum(-2.728564458083958E-147, 7.87584867255881E-164),
        DD.ofSum(6.911534513437014E-149, 2.3434612729446457E-166),
        DD.ofSum(-1.7507121442157682E-150, -3.610113360795167E-167),
        DD.ofSum(4.43460566672392E-152, -9.633244302445706E-169),
        DD.ofSum(-1.123298737848715E-153, 6.407992585647102E-170),
        DD.ofSum(2.8453489425693972E-155, -1.261702070020168E-171),
        DD.ofSum(-7.207353068415137E-157, 1.8745992502888998E-173),
        DD.ofSum(1.8256438595501416E-158, 1.1410294135768043E-174),
        DD.ofSum(-4.624409918974655E-160, -2.6688757360284575E-176),
        DD.ofSum(1.1713767165946985E-161, 5.7933128244158574E-179),
        DD.ofSum(-2.9671318854112526E-163, 2.6093750293274783E-179),
        DD.ofSum(7.515832866319735E-165, 5.996609111087652E-182),
        DD.ofSum(-1.9037827051837493E-166, -1.0499629634577001E-182),
        DD.ofSum(4.822337927175732E-168, -4.3170232110295926E-184),
        DD.ofSum(-1.221512466761957E-169, 2.050787849799304E-186),
    };

    /**
     * Compute the power function {@code (x+n)^s} in double precision.
     */
    @FunctionalInterface
    interface PowOperator {
        /**
         * Compute number {@code x+n} raised to the power {@code s}.
         *
         * @param x Number.
         * @param n Number.
         * @param s Power.
         * @return result
         */
        double apply(double x, int n, double s);
    }

    /**
     * Compute the power function {@code x^n} in double-double precision.
     */
    @FunctionalInterface
    interface ScaledPowOperator {
        /**
         * Compute number {@code x} raised to the power {@code n}.
         *
         * <p>The value is returned as fractional {@code f} and integral
         * {@code 2^exp} components.
         * <pre>
         * (x+xx)^n = (f+ff) * 2^exp
         * </pre>
         *
         * @param x Number.
         * @param n Power.
         * @param exp Result power of two scale factor (integral exponent).
         * @return Fraction part.
         * @see DD#frexp(int[])
         */
        DD apply(DD x, int n, long[] exp);
    }

    /** Context for the zeta implementation.
     * This is used to control how the test zeta implementation is computed.
     * Contains options for double, double-double and BigDecimal implementations. */
    private static final class Context {
        // Options flags may not be supported by all implementations

        /** Use extended precision power in the main power terms. */
        static final int EXT_POW_TERMS = 1;
        /** Use extended precision power in a series summation of power terms. */
        static final int EXT_POW_SERIES = 2;

        /** Use extended precision summation. */
        static final int EXT_SUM = 4;

        /** Return the result scaled by 2^106. This is the scale used for the
         * tail scaling option. This is only applied when tail scaling is enabled. */
        static final int SCALED_RESULT = 8;

        /** Use the power function to compute the power term in the tail series. */
        static final int TAIL_POW = 16;

        /** Use accurate add for the main terms in the DD implementation. */
        static final int AA_TERMS = 64;
        /** Use accurate add for the series summation in the DD implementation. */
        static final int AA_SERIES = 128;

        /** Option to use scaling in the tail series to prevent underflow. */
        static final int TAIL_SCALING = 256;

        /** Option to support summation of sub-normal terms. */
        static final int SUB_NORMAL_SUM = 512;

        // Default contexts optimised using the precision tests

        /** Context for double precision zeta implementation. */
        static final Context DOUBLE = Context.of(9, 15);
        /** Context for double precision zeta implementation for a reference
         * implementation where tail divergence will never occur.
         * The N and M are chosen using the test fixture
         * {@code testTailConvergenceLimit}.
         * Could use (50, 53) but this will be ~8% faster by using fewer
         * power operations. */
        static final Context DOUBLE_REF = Context.of(46, 39);

        /** Context for double precision using the DD zeta implementation.
         * Uses standard DD pow for all terms. Converge the tail beyond 53-bits. */
        static final Context DD_DOUBLE = Context.of(9, 15)
            .withEpsilon(0x1p-55).withOptions(1);
        /** Context for double-double precision using the DD zeta implementation.
         * Uses DDMath only for the most important terms. Note that using
         * DDMath for the series, or using accurate add do not improve the error
         * around the roots, or on the negative a test data. */
        static final Context DD_DOUBLE_DOUBLE = Context.of(15, 30)
            .withEpsilon(0x1p-106).withOptions(EXT_POW_TERMS);

        // BigDecimal implementation is very robust to different N.
        // Here N ~ M on the test data with precision +1 digit than the expected 17 per double.

        // TODO - re-optimise. Requires better test data

        /** Context for double precision using the BigDecimal zeta implementation. */
        static final Context BD_DOUBLE = Context.of(9, 15).withMathContext(new MathContext(18));
        /** Context for double-double precision using the BigDecimal zeta implementation. */
        static final Context BD_DOUBLE_DOUBLE = Context.of(17, 30).withMathContext(new MathContext(36));
        /** Context for quad-double precision using the BigDecimal zeta implementation. */
        static final Context BD_QUAD_DOUBLE = Context.of(37, 50).withMathContext(new MathContext(70));

        /** Default epsilon for the series summations. */
        private static final double EPS = 0x1p-53;
        /** Default math context for extended precision BigDecimal implementation. */
        private static final MathContext MC = MathContext.DECIMAL128;
        /** Default options. Configured for the double implementation:
         * use extended precision pow everywhere; extended precision sum;
         * tail scaling to reduce underflow; and support sub-normal summation. */
        private static final int DEFAULT_OPTIONS = EXT_POW_TERMS | EXT_POW_SERIES | EXT_SUM | TAIL_SCALING | SUB_NORMAL_SUM;

        /** N. */
        private final int n;
        /** M. */
        private final int m;
        /** Epsilon for convergence of series summations. */
        private final double eps;
        /** Math context for extended precision evaluations. */
        private final MathContext mc;
        /** Options flags. */
        private final int options;

        /**
         * Create an instance.
         *
         * @param n the number of terms for the series sum
         * @param m the number of terms for the tail sum
         * @param eps the epsilon for series summations
         * @param mc the math context for BigDecimal evaluations
         * @param options the options
         */
        private Context(int n, int m, double eps, MathContext mc,
            int options) {
            // Assume all tables for the tail series are the same size
            if (m > RF.length) {
                throw new IllegalArgumentException("Unsupported M: " + m);
            }
            this.n = n;
            this.m = m;
            this.eps = eps;
            this.mc = mc;
            this.options = options;
        }

        /**
         * Create a context.
         *
         * @param n the number of terms N.
         * @param m the number of terms M.
         * @return the context
         */
        static Context of(int n, int m) {
            return new Context(n, m, EPS, MC, DEFAULT_OPTIONS);
        }

        /**
         * Return a context that uses the provided number of terms N.
         *
         * @param value the value
         * @return the context
         */
        Context withN(int value) {
            return new Context(value, m, eps, mc, options);
        }

        /**
         * Return a context that uses the provided convergence epsilon in summations.
         *
         * @param value the value
         * @return the context
         */
        Context withEpsilon(double value) {
            return new Context(n, m, value, mc, options);
        }

        /**
         * Return a context that uses the provided MathContext.
         *
         * @param value the value
         * @return the context
         */
        Context withMathContext(MathContext value) {
            return new Context(n, m, eps, value, options);
        }

        /**
         * Return a context that uses the provided options.
         *
         * @param value the value
         * @return the context
         */
        Context withOptions(int value) {
            return new Context(n, m, eps, mc, value);
        }

        /**
         * Return a context that adds the provided options to the current options.
         *
         * @param value the value
         * @return the context
         */
        Context addOptions(int value) {
            return new Context(n, m, eps, mc, options | value);
        }

        /**
         * Gets N.
         * @return n
         */
        int getN() {
            return n;
        }

        /**
         * Gets M.
         * @return m
         */
        int getM() {
            return m;
        }

        /**
         * Gets the convergence epsilon for series summations.
         * @return the epsilon
         */
        double getEps() {
            return eps;
        }

        /**
         * Gets the math context for extended precision evaluations.
         * @return the math context
         */
        MathContext getMathContext() {
            return mc;
        }

        /**
         * Checks if the option is set. All bits in the option must be enabled
         * in the configured options.
         *
         * @param option the option
         * @return true if all the option bits are enabled
         */
        boolean isSet(int option) {
            return (options & option) == option;
        }

        /**
         * Extended precision power function {@code (x+y)^z}.
         *
         * <p>Warning: assumes x+y is finite; z is negative and finite.
         * This method is used for testing where the arguments should be {@code x + n}
         * with x finite and less than 2^53; and n a small integer.
         *
         * @param x the argument x ({@code < 2^53})
         * @param y the argument y
         * @param z the exponent z ({@code < 0})
         * @return the result
         */
        static double extPow(double x, int y, double z) {
            // Use a test function with support for (a+n) > 0.
            // This must handle possible overflow when (a+n) < 1 and z is -large.
            return PowXPNTest.powxpn(x, y, z);
        }

        /**
         * Standard precision {@code (x+y)^z}.
         *
         * @param x the argument x.
         * @param y the argument y.
         * @param z the exponent z.
         * @return the result
         */
        static double stdPow(double x, int y, double z) {
            return Math.pow(x + y, z);
        }
    }

    /**
     * Compute the sum of {@code (x, xx)} and {@code (y, yy)}.
     *
     * <p><strong>*** This code has been extracted from the DD class for testing***</strong>.
     *
     * <p>It is not part of the public API. The method has been adapted to use
     * DD arguments and DD.ofSum for all twoSum and fastTwoSum operations. This makes it
     * slower than an optimal implementation that has access to the private DD constructor
     * (and fastTwoSum to compute the round-off).
     *
     * <p>The high-part of the result is within 1 ulp of the true sum {@code e}. The
     * low-part of the result is within 1 ulp of the result of the high-part subtracted
     * from the true sum {@code e - hi}.
     *
     * <p>Note: This is an internal helper method used when accuracy is required. The
     * computed result is within 1 eps of the exact result where eps is 2<sup>-106</sup>.
     * The performance is approximately 2-fold slower than {@link DD#add(DD)} (perhaps
     * more given the lack of fastTwoSum).
     *
     * @param x (x, xx).
     * @param y (y, yy).
     * @return the sum
     */
    private static DD accurateAdd(DD x, DD y) {
        // Expansion sum (Schewchuk Fig 7): (x, xx) + (x, yy) -> (s0, s1, s2, s3)
        DD s = DD.ofSum(x.lo(), y.lo());
        double s3 = s.lo();
        s = DD.ofSum(x.hi(), s.hi());
        // (s0, s1, s2) == (s.hi(), s.lo(), s3)
        double s0 = s.hi();
        s = DD.ofSum(s.lo(), y.hi());
        double s2 = s.lo();
        s = DD.ofSum(s0, s.hi());
        // s1 = s.lo()
        s0 = s.hi();
        // Compress (Schewchuk Fig. 15) (s0, s1, s2, s3) -> (s0, s1)
        // fastTwoSum
        s = DD.ofSum(s.lo(), s2);
        final double s1 = s.hi();
        // fastTwoSum
        s = DD.ofSum(s.lo(), s3);
        // s2 = s.hi()
        s3 = s.lo();
        // fastTwoSum
        s = DD.ofSum(s1, s.hi());
        s2 = s.lo();
        // fastTwoSum
        s = DD.ofSum(s0, s.hi());
        // Here (s0, s1) = s
        // e = exact 212-bit result
        // |e - s0| <= ulp(s0)
        // |s1 + s2 + s3| <= ulp(e - s0)   (Sum magnitudes small to high)
        // fastTwoSum
        return DD.ofSum(s.hi(), s3 + s2 + s.lo());
    }

    /** Define the expected error for a test. */
    private interface TestError {
        /**
         * @return maximum allowed error
         */
        double getTolerance();

        /**
         * @return maximum allowed RMS error
         */
        double getRmsTolerance();
    }

    /** Define a test for an extended precision zeta function.
     * Results are computed as BigDecimal. Errors are reported using the
     * scaled relative error 2^b * (actual - expected) / expected. */
    private interface ExtendedPrecisionTestCase extends TestError {
        /**
         * @return function to test
         */
        BiFunction<Integer, Double, BigDecimal> getFunction();

        /**
         * @return Filenames of the test data
         */
        String[] getFilenames();

        /**
         * Get the scale for the relative error. Use 52 for the approximate equivalent
         * of 1 ULP in double precision; 105 for double-double precision.
         * Note that arbitrary precision computation using BigDecimal does
         * not have an equivalent to a ULP as the mantissa is not fixed. Scaling
         * the relative error is for convenience.
         *
         * @return Scale to apply to the relative error.
         * @see Math#scalb(double, int)
         */
        int scale();
    }

    /** Define a test for a double precision zeta function. */
    private interface DoublePrecisionTestCase extends TestError {
        /**
         * @return function to test
         */
        DoubleBinaryOperator getFunction();

        /**
         * @return Filenames of the test data
         */
        String[] getFilenames();
    }

    /**
     * Define the test cases for the zeta function.
     * This encapsulates the function to test, the expected maximum and RMS error, and
     * the resource file(s) containing the data.
     */
    private enum ZetaTestCase implements DoublePrecisionTestCase {
        // Prefix T_ for test code in this class.
        // Note: Testing code may skip tests so here the RMS tolerance is configured
        // to be more lenient. It may approach the max depending on the included cases.

        // Testing implementations for negative a.
        // Requires accurate evaluation with integer s (IS).

        // The double implementation is accurate until zeta values are sub-normal
        T_DOUBLE_ZETA_IS((s, a) -> HurwitzZetaTest.zeta(s, a, Context.DOUBLE), INT_TEST_RESOURCES, 1.5, 0.75),
        // Standard double precision sum and power function works
        // with sub-normal term scaling
        T_ZETA_2_POW_N1022_STD((s, a) -> HurwitzZetaTest.zeta(s, a,
            Context.DOUBLE.withOptions(Context.SUB_NORMAL_SUM + Context.TAIL_SCALING)
        ), "hzeta_2p-1022.csv", 40, 0.95),
        T_ZETA_DOMAIN_2_POW_N1022_STD((s, a) -> HurwitzZetaTest.zeta(s, a,
            Context.DOUBLE.withOptions(Context.SUB_NORMAL_SUM + Context.TAIL_SCALING)
        ), "hzeta_domain_2p-1022.csv", 65, 0.95),
        // These are within 0.5 ULP as a double-double but rounding put single errors at just over 0.5 ulp.
        // The zeta values are far below the min double value.
        T_BD_ZETA_IS((s, a) -> HurwitzZetaTest.zeta((int) s, new BigDecimal(a), null, Context.BD_DOUBLE).doubleValue(),
            EXT_INT_TEST_RESOURCES, 0.56, 0.2),
        T_DD_ZETA_IS((s, a) -> {
            final int[] exp = {0};
            final DD r = HurwitzZetaTest.zeta((int) s, DD.of(a), null, exp, Context.DD_DOUBLE);
            return Math.scalb(r.hi(), exp[0]);
        }, EXT_INT_TEST_RESOURCES, 0.52, 0.2),

        // No cancellation - all implementations work
        T_DOUBLE_ZETA_S2_4_N_A1_15((s, a) -> HurwitzZetaTest.zetaNegative((int) s, a, false), "hzeta_s2_4_na1_15_p0.5_p0x1p-1.csv", 1.5, 0.75),
        T_DOUBLE_ZETA_S2_4_N_A40_99((s, a) -> HurwitzZetaTest.zetaNegative((int) s, a, false), "hzeta_s2_4_na40_99_p0.5_p0x1p-1.csv", 2.0, 0.75),
        // These use a context to evaluate to double precision
        T_BD_ZETA_S2_4_N_A1_15((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_s2_4_na1_15_p0.5_p0x1p-1.csv", 0.57, 0.2),
        T_BD_ZETA_S2_4_N_A40_99((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_s2_4_na40_99_p0.5_p0x1p-1.csv", 0.57, 0.2),
        T_DD_ZETA_S2_4_N_A1_15((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_s2_4_na1_15_p0.5_p0x1p-1.csv", 0.57, 0.2),
        T_DD_ZETA_S2_4_N_A40_99((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_s2_4_na40_99_p0.5_p0x1p-1.csv", 0.57, 0.2),

        // Double arithmetic is max error ~24-bits when cancellation is expected
        // Note: Testing code may skip tests so here the RMS is configured to the same as the max
        T_DOUBLE_ZETA_S3_5_N_A1_15_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegative((int) s, a, false), "hzeta_s3_5_na1_15_p0.5_p0x1p-30.csv", 0x1p24, 0x1p24),
        // Computing the largest terms in extended precision error improves 6 bits
        T_DOUBLE_P_ZETA_S3_5_N_A1_15_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegative((int) s, a, true), "hzeta_s3_5_na1_15_p0.5_p0x1p-30.csv", 0x1p17, 0x1p17),

        // Extended precision implementations
        T_BD_ZETA_S3_5_N_A1_15_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_s3_5_na40_99_p0.5_p0x1p-30.csv", 0, 0),
        T_BD_ZETA_S3_5_N_A40_99_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_s3_5_na40_99_p0.5_p0x1p-30.csv", 0, 0),
        T_DD_ZETA_S3_5_N_A1_15_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_s3_5_na40_99_p0.5_p0x1p-30.csv", 0, 0),
        T_DD_ZETA_S3_5_N_A40_99_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_s3_5_na40_99_p0.5_p0x1p-30.csv", 0, 0),

        // Roots are the point of maximum cancellation
        // Double arithmetic has only a few bits of precision on average; and may have no bits correct
        T_DOUBLE_ZETA_ROOT_S3_21_N_A0_100((s, a) -> HurwitzZetaTest.zetaNegative((int) s, a, false), "hzeta_root_s3_21_na2b0_2b6.csv", 0x1p55, 0x1p52),
        // Computing the largest terms in extended precision allows a few bits of precision
        T_DOUBLE_P_ZETA_ROOT_S3_21_N_A0_100((s, a) -> HurwitzZetaTest.zetaNegative((int) s, a, true), "hzeta_root_s3_21_na2b0_2b6.csv", 0x1p48, 0x1p45),

        T_BD_ZETA_ROOTS((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), ROOT_TEST_RESOURCES, 0, 0),
        // 1 ULP on the case of total cancellation: 3.0, -2.4994443912584825
        T_DD_ZETA_ROOTS((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), ROOT_TEST_RESOURCES, 0.7, 0.7),

        // Final implementation
        ZETA_S1_4_A1_8(HurwitzZeta::value, "hzeta_s1_4_a1_8.csv", 1.5, 0.3),
        ZETA_S1_4_A8_32(HurwitzZeta::value, "hzeta_s1_4_a8_32.csv", 2, 0.49),
        ZETA_S1_4_A32_2147483648(HurwitzZeta::value, "hzeta_s1_4_a32_2147483648.csv", 1.8, 0.5),
        ZETA_S1_4_A0_1(HurwitzZeta::value, "hzeta_s1_4_a0_1.csv", 1.5, 0.35),
        ZETA_S1_4_A0(HurwitzZeta::value, "hzeta_s1_4_a1e-16_1e-14.csv", 1.25, 0.22),
        ZETA_S4_32_A1_8(HurwitzZeta::value, "hzeta_s4_32_a1_8.csv", 1.3, 0.4),
        ZETA_S1_19_A2147483648_1E16(HurwitzZeta::value, "hzeta_s1_19_a2147483648_10000000000000000.csv", 1.8, 0.5),
        // Data that tests the convergence of the tail series with large s and s ~ (a+N).
        // Most of the results are sub-normal. The final implementation should compute zero
        // and not create +/-infinity from overflow of the tail series Pochammer term.
        ZETA_CONVERGENCE(HurwitzZeta::value, EXT_INT_TEST_RESOURCES, 3, 0.25),
        // TODO: These may not be very useful.
        // We require a set of data where a^-s is normal and (a+1)^-s is sub-normal
        // Terms of the series (a+k)^-s are sub-normal
        // 1.33% of results are sub-normal
        ZETA_S20_1024_AP9_SUB_NORMAL_1023(HurwitzZeta::value, "hzeta_s20_1024_ap9_subnormal1023.csv", 1.95, 0.3),
        // 38.4% of results are sub-normal
        ZETA_S20_1024_AP1_SUB_NORMAL_1040(HurwitzZeta::value, "hzeta_s20_1024_ap1_subnormal1040.csv", 1.55, 0.2),
        // 48.3% of results are sub-normal
        ZETA_S20_1024_AP1_SUB_NORMAL_1053(HurwitzZeta::value, "hzeta_s20_1024_ap1_subnormal1053.csv", 1, 0.15),
        // 100% of results are sub-normal
        ZETA_S20_1024_AP0_1_SUB_NORMAL_1073(HurwitzZeta::value, "hzeta_s20_1024_ap0.1_subnormal1073.csv", 0.6, 0.1),
        // Full domain range that evaluates at 2^-1022 and 2^-1065
        ZETA_2_POW_N1022(HurwitzZeta::value, "hzeta_2p-1022.csv", -1.35, 0.35),
        ZETA_2_POW_N1065(HurwitzZeta::value, "hzeta_2p-1065.csv", 0.6, 0.1),
        ZETA_2_POW_N1075(HurwitzZeta::value, "hzeta_2p-1075.csv", 0.6, 0.1),
        // Random samples from the domain
        ZETA_DOMAIN(HurwitzZeta::value, DOMAIN_RESOURCES, 2.0, 0.5),
        // Negative a
        ZETA_S2_4_N_A1_15(HurwitzZeta::value, "hzeta_s2_4_na1_15_p0.5_p0x1p-1.csv", 0.55, 0.07),
        ZETA_S2_4_N_A40_99(HurwitzZeta::value, "hzeta_s2_4_na40_99_p0.5_p0x1p-1.csv", 0.55, 0.07),
        ZETA_S2_4_N_A1_15_B30(HurwitzZeta::value, "hzeta_s2_4_na1_15_p0x1p-30.csv", 0, 0),
        ZETA_S2_4_N_A1_15_HALF_B30(HurwitzZeta::value, "hzeta_s2_4_na1_15_p0.5_p0x1p-30.csv", 0.53, 0.07),
        ZETA_S3_5_N_A1_15(HurwitzZeta::value, "hzeta_s3_5_na1_15_p0.5_p0x1p-1.csv", 0, 0),
        ZETA_S3_5_N_A40_99(HurwitzZeta::value, "hzeta_s3_5_na40_99_p0.5_p0x1p-1.csv", 0, 0),
        ZETA_S3_5_N_A1_15_B30(HurwitzZeta::value, "hzeta_s3_5_na1_15_p0x1p-30.csv", 0, 0),
        ZETA_S3_5_N_A1_15_HALF_B30(HurwitzZeta::value, "hzeta_s3_5_na1_15_p0.5_p0x1p-30.csv", 0, 0),
        ZETA_S3_5_N_A40_99_HALF_B30(HurwitzZeta::value, "hzeta_s3_5_na40_99_p0.5_p0x1p-30.csv", 0, 0),
        ZETA_S6_50_N_A400_999_HALF_B3(HurwitzZeta::value, "hzeta_s6_50_na400_999_p0.5_p0x1p-3.csv", 0.52, 0.1),
        ZETA_S7_51_N_A400_999_HALF_B3(HurwitzZeta::value, "hzeta_s7_51_na400_999_p0.5_p0x1p-3.csv", 1, 0),
        ZETA_ROOTS(HurwitzZeta::value, ROOT_TEST_RESOURCES, 0.7, 0.05);

//        JDK Temurin 25.492-b09
//        T_DOUBLE_ZETA_IS                         max        1.05975   RMS       0.260675   mean    -0.00757331  n  6000  (392ms)
//        T_ZETA_2_POW_N1022_STD                   max        36.0615   RMS       0.881910   mean      0.0364887  n  1998  (116ms)
//        T_ZETA_DOMAIN_2_POW_N1022_STD            max        58.8128   RMS       0.883161   mean    -0.00500754  n 50000  (714ms)
//        T_BD_ZETA_IS                             max       0.552489   RMS      0.0369241   mean    0.000364038  n 16836  (1.33e+03ms)
//        T_DD_ZETA_IS                             max       0.509511   RMS      0.0149965   mean   -2.98738e-05  n 16836  (328ms)
//        T_DOUBLE_ZETA_S2_4_N_A1_15               max        1.43054   RMS       0.363239   mean    -0.00657329  n  3000  (59.0ms)
//        T_DOUBLE_ZETA_S2_4_N_A40_99              max        1.40007   RMS       0.398324   mean    -0.00538358  n  3000  (81.8ms)
//        T_BD_ZETA_S2_4_N_A1_15                   max       0.550507   RMS      0.0510756   mean    0.000694938  n  3000  (171ms)
//        T_BD_ZETA_S2_4_N_A40_99                  max       0.566598   RMS      0.0524151   mean   -0.000150341  n  3000  (285ms)
//        T_DD_ZETA_S2_4_N_A1_15                   max       0.504587   RMS      0.0159107   mean   -0.000503140  n  3000  (49.7ms)
//        T_DD_ZETA_S2_4_N_A40_99                  max       0.507904   RMS      0.0260235   mean   -0.000670253  n  3000  (69.6ms)
//        T_DOUBLE_ZETA_S3_5_N_A1_15_HALF_B30      max    8.40748e+06   RMS    1.04226e+06   mean        15829.3  n  3000  (61.8ms)
//        T_DOUBLE_P_ZETA_S3_5_N_A1_15_HALF_B30    max        87729.7   RMS        7926.41   mean        93.2080  n  3000  (55.0ms)
//        T_BD_ZETA_S3_5_N_A1_15_HALF_B30          max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (528ms)
//        T_BD_ZETA_S3_5_N_A40_99_HALF_B30         max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (428ms)
//        T_DD_ZETA_S3_5_N_A1_15_HALF_B30          max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (57.2ms)
//        T_DD_ZETA_S3_5_N_A40_99_HALF_B30         max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (55.4ms)
//        T_DOUBLE_ZETA_ROOT_S3_21_N_A0_100        max    1.30515e+16   RMS    6.71498e+14   mean    2.99768e+13  n   642  (37.5ms)
//        T_DOUBLE_P_ZETA_ROOT_S3_21_N_A0_100      max    9.57925e+13   RMS    6.03802e+12   mean   -2.72521e+11  n   642  (18.2ms)
//        T_BD_ZETA_ROOTS                          max        0.00000   RMS        0.00000   mean        0.00000  n  2640  (585ms)
//        T_DD_ZETA_ROOTS                          max       0.638297   RMS      0.0158942   mean   -4.88173e-05  n  2640  (72.6ms)
//        ZETA_S1_4_A1_8                           max        1.29125   RMS       0.287148   mean     0.00943939  n  3000  (89.3ms)
//        ZETA_S1_4_A8_32                          max        1.74869   RMS       0.361734   mean     0.00183242  n  3000  (48.0ms)
//        ZETA_S1_4_A32_2147483648                 max        1.77377   RMS       0.467224   mean     0.00262765  n  3000  (32.4ms)
//        ZETA_S1_4_A0_1                           max        1.18864   RMS       0.316203   mean     0.00322568  n  3000  (31.7ms)
//        ZETA_S1_4_A0                             max       0.952404   RMS       0.108517   mean    -0.00255755  n  3000  (43.5ms)
//        ZETA_S4_32_A1_8                          max        1.03142   RMS       0.331588   mean     0.00607164  n  3000  (41.4ms)
//        ZETA_S1_19_A2147483648_1E16              max        1.76616   RMS       0.473537   mean      0.0109014  n  3000  (39.5ms)
//        ZETA_CONVERGENCE                         max        1.04351   RMS       0.169635   mean    -0.00212817  n 16836  (194ms)
//        ZETA_S20_1024_AP9_SUB_NORMAL_1023        max        1.89978   RMS       0.250231   mean     0.00651038  n  3000  (37.7ms)
//        ZETA_S20_1024_AP1_SUB_NORMAL_1040        max        1.48305   RMS       0.159991   mean     0.00413054  n  3000  (28.0ms)
//        ZETA_S20_1024_AP1_SUB_NORMAL_1053        max       0.895140   RMS       0.103818   mean    7.07467e-05  n  3000  (26.7ms)
//        ZETA_S20_1024_AP0_1_SUB_NORMAL_1073      max       0.500519   RMS     0.00913819   mean    0.000166840  n  3000  (27.6ms)
//        ZETA_2_POW_N1022                         max        1.27858   RMS       0.270062   mean      0.0120844  n  1998  (23.8ms)
//        ZETA_2_POW_N1065                         max        0.00000   RMS        0.00000   mean        0.00000  n  1998  (27.8ms)
//        ZETA_DOMAIN                              max        1.96329   RMS       0.428365   mean     0.00334990  n 50000  (531ms)
//        ZETA_S2_4_N_A1_15                        max       0.518098   RMS      0.0402256   mean    -0.00218941  n  3000  (54.3ms)
//        ZETA_S2_4_N_A40_99                       max       0.527688   RMS      0.0523742   mean    -0.00135787  n  3000  (86.5ms)
//        ZETA_S2_4_N_A1_15_B30                    max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (28.5ms)
//        ZETA_S2_4_N_A1_15_HALF_B30               max       0.511209   RMS      0.0503863   mean   -0.000326508  n  3000  (36.4ms)
//        ZETA_S3_5_N_A1_15                        max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (60.7ms)
//        ZETA_S3_5_N_A40_99                       max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (94.8ms)
//        ZETA_S3_5_N_A1_15_B30                    max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (22.5ms)
//        ZETA_S3_5_N_A1_15_HALF_B30               max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (25.8ms)
//        ZETA_S3_5_N_A40_99_HALF_B30              max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (50.3ms)
//        ZETA_S6_50_N_A400_999_HALF_B3            max       0.517315   RMS      0.0207425   mean   -0.000509997  n  3000  (40.9ms)
//        ZETA_S7_51_N_A400_999_HALF_B3            max        0.00000   RMS        0.00000   mean        0.00000  n  3000  (56.9ms)
//        ZETA_ROOTS                               max       0.638297   RMS      0.0158942   mean   -4.88173e-05  n  2640  (56.7ms)

        /** The function. */
        private final DoubleBinaryOperator fun;

        /** The filenames containing the test data. */
        private final String[] filenames;

        /** The maximum allowed ulp. */
        private final double maxUlp;

        /** The maximum allowed RMS ulp. */
        private final double rmsUlp;

        /**
         * Create an instance.
         *
         * @param fun function to test
         * @param filename Filename of the test data
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         */
        ZetaTestCase(DoubleBinaryOperator fun, String filename, double maxUlp, double rmsUlp) {
            this(fun, new String[] {filename}, maxUlp, rmsUlp);
        }

        /**
         * Create an instance.
         *
         * @param fun function to test
         * @param filenames Filenames of the test data
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         */
        ZetaTestCase(DoubleBinaryOperator fun, String[] filenames, double maxUlp, double rmsUlp) {
            this.fun = fun;
            this.filenames = filenames;
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
        }

        @Override
        public DoubleBinaryOperator getFunction() {
            return fun;
        }

        @Override
        public String[] getFilenames() {
            return filenames;
        }

        @Override
        public double getTolerance() {
            return maxUlp;
        }

        @Override
        public double getRmsTolerance() {
            return rmsUlp;
        }
    }
    /**
     * Define the test cases for the extended precision zeta function.
     * This encapsulates the function to test, the expected maximum and RMS error, and
     * the resource file(s) containing the data.
     */
    private enum ExtZetaTestCase implements ExtendedPrecisionTestCase {
        // Prefix T_ for test code in this class.

        T_DD_ZETA_2_POW_N1075((s, a) -> HurwitzZetaTest.zeta(s, DD.of(a), Context.DD_DOUBLE_DOUBLE),
            "hzeta_is_2p-1075.csv", 20, 1.5, 106),
        T_DD_ZETA_2_POW_N1075_SERIES((s, a) -> HurwitzZetaTest.zeta(s, DD.of(a),
            Context.DD_DOUBLE_DOUBLE.addOptions(Context.EXT_POW_SERIES)),
            "hzeta_is_2p-1075.csv", 2, 0.25, 106),
        T_DD_ZETA_2_POW_N1075_AA((s, a) -> HurwitzZetaTest.zeta(s, DD.of(a),
            Context.DD_DOUBLE_DOUBLE.addOptions(Context.EXT_POW_SERIES | Context.AA_SERIES | Context.AA_TERMS)),
            "hzeta_is_2p-1075.csv", 0.7, 0.2, 106),
        T_BD_ZETA_2_POW_N1075((s, a) -> HurwitzZetaTest.zeta(s, new BigDecimal(a), null, Context.BD_DOUBLE_DOUBLE),
            "hzeta_is_2p-1075.csv", 0.001, 0.0005, 106);

//        JDK Temurin 25.492-b09
//        T_DD_ZETA_2_POW_N1075                    max        19.9135   RMS        1.25660   mean       0.161528  n  1836  (320ms)
//        T_DD_ZETA_2_POW_N1075_SERIES             max        1.51680   RMS       0.210924   mean    -0.00947771  n  1836  (112ms)
//        T_DD_ZETA_2_POW_N1075_AA                 max       0.655684   RMS       0.171493   mean    -0.00829713  n  1836  (97.6ms)
//        T_BD_ZETA_2_POW_N1075                    max    0.000575613   RMS    0.000118868   mean    3.64535e-06  n  1836  (642ms)

        /** The function. */
        private final BiFunction<Integer, Double, BigDecimal> fun;

        /** The filenames containing the test data. */
        private final String[] filenames;

        /** The maximum allowed ulp. */
        private final double maxUlp;

        /** The maximum allowed RMS ulp. */
        private final double rmsUlp;

        /** Scale to apply to the relative error. */
        private final int scale;

        /**
         * Create an instance.
         *
         * @param fun function to test
         * @param filename Filename of the test data
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         * @param scale scale to apply to the relative error
         */
        ExtZetaTestCase(BiFunction<Integer, Double, BigDecimal> fun, String filename,
            double maxUlp, double rmsUlp, int scale) {
            this(fun, new String[] {filename}, maxUlp, rmsUlp, scale);
        }

        /**
         * Create an instance.
         *
         * @param fun function to test
         * @param filenames Filenames of the test data
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         * @param scale scale to apply to the relative error
         */
        ExtZetaTestCase(BiFunction<Integer, Double, BigDecimal> fun, String[] filenames,
            double maxUlp, double rmsUlp, int scale) {
            this.fun = fun;
            this.filenames = filenames;
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
            this.scale = scale;
        }

        @Override
        public BiFunction<Integer, Double, BigDecimal> getFunction() {
            return fun;
        }

        @Override
        public String[] getFilenames() {
            return filenames;
        }

        @Override
        public int scale() {
            return scale;
        }

        @Override
        public double getTolerance() {
            return maxUlp;
        }

        @Override
        public double getRmsTolerance() {
            return rmsUlp;
        }
    }

    // Test zeta implementation
    //
    // This class contains a parameterized version of the final implementation.
    // A context is used to configure the N and M parameters and enable features.
    // With all features enabled the result should match the final version.
    //
    // Other features were tested during development and removed from this version
    // to simplify the code:
    //
    // - Use of the power function to compute the tail terms. This did not improve
    // precision as the tail term is not required to ULP accuracy. Recycling an
    // existing power term is faster. Use of power for the tail term is an option in
    // the DD version where the power function can never underflow.
    // - Use of an extended precision sum in the tail. The tail term is not
    // required to high accuracy.
    // - Different updates to the tail power term. Either divide by (a+n) twice
    // or multiply by (a+n)^-2.
    // - Use of log-sum-exp to support sub-normal terms. This identifies the maximum
    // of a^-s and (a+n)^(1-s)/(s-1). This term dominates the result and all other
    // power terms are smaller and can be sub-normal. Remap the power function. In
    // this example a^-s is the major term:
    // pow(a + k, -s) = exp(-s * log(a + k))
    //                = exp(-s * log(a)) * exp(-s * log(a + k) + s * log(a))  [scaled]
    // The result can control error but does not outperform an implementation that scales
    // sub-normal terms. It may have higher runtime cost due to use of log+exp vs pow.
    // Also note that exp(log(x) * y) is not as accurate as pow(x, y).
    // - Computation of I = (a+n)^-(1+s) as (a+n)^-s * (a+n). When this term is the
    // major term this small modification increases error.
    //
    // Note: The method is sensitive to the initial loop over N to create S.
    // Under certain conditions the N cannot be too high if using an ascending
    // sum of k as the sum does not converge and later terms are added with
    // low precision.
    //
    // Better results are obtained using descending k. However this prevents
    // an early exit if the series is rapidly converging and the term (a+k)^-s
    // drops below machine epsilon of the sum.
    //
    // The tail series may be divergent depending on s and (a+N). This
    // can be mitigated by using a small M for the tail series before it diverges,
    // or increasing N so that the series sum contains most of the result
    // when s is large and the term (a+k)^-s is strongly converging. The tail
    // then contributes a negligible magnitude and is ignored when multiplied
    // by the scaling factor (a+N)^-s.

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     * See {@link HurwitzZeta} for the formula details.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    static double zeta(double s, double a, Context c) {
        return zeta(s, a, Double.NaN, c);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     * See {@link HurwitzZeta} for the formula details.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @param a0 {@code a^-s} (can be nan)
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    static double zeta(double s, double a, double a0, Context c) {
        // First term may be provided
        double t0 = Double.isNaN(a0) ? Math.pow(a, -s) : a0;
        lastM = 0;

        // For testing the tail convergence this supports return of the result
        // using the scaling of the tail series
        final double resultScale = c.isSet(Context.SCALED_RESULT | Context.TAIL_SCALING) ?
            UPSCALE_106 : 1;

        // Can overflow if 0 <= a < 1
        if (t0 >= LARGE_T0) {
            // Nothing can be added in double-double precision
            return t0 * resultScale;
        }

        int n = c.getN();

        // Asymptotic behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0 (no series sum S): the I term and the first term of T.
        // Threshold so (a+N) <= 2^53.
        //
        // Tiny s trap using N=0.
        // zeta ~ a^(1-s) / (s-1) + a^-s / 2
        //      = a*a^-s / (s-1) + a^s / 2
        //      = a^s * (a / (s-1) + 0.5)
        //      = a^(1-s) * (1 / (s-1) + 1 / 2a)
        // When 0.5 cannot be added to a / (s-1) then N=0 is a good approximation.
        // a / (s - 1) > 2^52 : a > 2^52 * (s - 1)
        // As s -> 1 then a can reduce in magnitude until a ~ 1 and s ~ 1: a^-s -> 1 / a.
        // TODO: Test which option is better?
        // Relying on ti could be used before t0 is even computed.
        if (a > LARGE_APN - n || a * c.getEps() > s - 1) {
            // Note: If a^-s is sub-normal then a^(1-s)=(a^-s * a) may not be
            // sub-normal as a is large. Compute using pow.
            double ti;
            if (c.isSet(Context.SUB_NORMAL_SUM) && t0 < Double.MIN_NORMAL) {
                // Large a. The tiny s trap cannot reach here as t0 -> 1/a.
                // Avoid computing the division with a sub-normal term
                if (t0 * a < Double.MIN_NORMAL) {
                    ti = Math.pow(a, (1 - s) * 0.5);
                    ti = 2 * ti * UPSCALE_916 * ti / (s - 1);
                } else {
                    ti = 2 * UPSCALE_916 * Math.pow(a, 1 - s) / (s - 1);
                }
                // Recompute t0. Is this worth the cost? The maximum overlap
                // between ti and t0 * 0.5 is small.
                t0 = Math.pow(a, -s * 0.5);
                t0 = t0 * UPSCALE_916 * t0;
                //t0 *= UPSCALE_916;
                return (ti + t0) * 0.5 * (DOWNSCALE_916 * resultScale);

//                // Do not recompute t0 and rely on ti having double precision
//                if (t0 * a < Double.MIN_NORMAL) {
//                    ti = Math.pow(a, (1 - s) * 0.5);
//                    ti = ti * UPSCALE_1022 * ti;
//                } else {
//                    ti = UPSCALE_1022 * Math.pow(a, 1 - s);
//                }
//                return ti * (1 / (s - 1) + 1 / (2 * a)) * Double.MIN_NORMAL;
            }
            // Compute (2 * ti + t0) * 0.5 to avoid (ti + t0 * 0.5) losing a bit on sub-normal t0
            ti = 2 * Math.pow(a, 1 - s) / (s - 1);
            return (ti + t0) * 0.5 * resultScale;
        }

        // Now any (a+n)^-s cannot overflow and the sum cannot overflow

        // Methods are used to increase precision
        // 1. Extended precision summation (low cost)
        // 2. Extended precision power function using round-off from (a+k) (selectively applied)
        // 3. Computation of sub-normal power terms in parts (low cost)
        // 4. Scaling for sub-normal tail terms (no cost)
        // 5. Adding series terms in magnitude order (no cost)
        //
        // Options 1-4 must be enabled by the context.
        //
        // Option 5 requires reversal of the summation loop and is not supported.
        // Testing shows that forward summation loses precision unless performed in
        // extended precision. The sum must stop not at 2^-53 of the current sum but well
        // past it (e.g. another 10 bits) otherwise a slowly converging sum will have
        // errors without computing I and T.
        // TODO:
        // Create a test to sample within the zeta domain for each power of 2 below 2^53
        // and histogram the N required for convergence compared to the N used by the
        // current implementation.

        PowOperator pow = c.isSet(Context.EXT_POW_TERMS) ? Context::extPow : Context::stdPow;
        // Power function for the series (allows switching to a faster pow)
        PowOperator powSeries = c.isSet(Context.EXT_POW_SERIES) ? Context::extPow : Context::stdPow;

        // Check the extended precision power can be used.
        // If s is large or there is no round-off from (a+n) use standard precision.
        if (s >= Z_LIMIT || DD.ofSum(a, n).lo() == 0) {
            pow = powSeries = Context::stdPow;
        }

        // TODO:
        // Find a data line where the first term is normal and the next term is sub-normal.
        // Do we have this already?

        // Scaled power for the series sum to support sub-normal terms.
        // This includes early-exit traps for when s is too large and the value is zero or t0.
        PowOperator powS = powSeries;
        boolean scaledSum = false;
        final double originalT0 = t0;
        double t1;
        if (c.isSet(Context.SUB_NORMAL_SUM)) {
            // Check t0
            // Recompute scaled t0 if applicable.
            if (t0 < Double.MIN_NORMAL) {
                t0 = Math.pow(a, -s * 0.5);
                t0 = t0 * UPSCALE_1022 * t0;
                if (t0 == 0) {
                    // I = (a+n)^(1-s) / (s-1) = (a+n)^-s * (a+n) / (s-1)
                    // For any s > 1 the I term is not 2^1022 larger than a^-s.
                    // The max is ~ 2^53 with a~2^53 and s~1.
                    // No need to compute I
                    return 0;
                }
                scaledSum = true;
                // Redefine the series power function
                PowOperator originalPowSeries = powSeries;
                powS = (x, y, z) -> {
                    final double p = originalPowSeries.apply(x, y, z * 0.5);
                    return p * UPSCALE_1022 * p;
                };
            }

            // Compute t1 = (a+1)^-s. This may be scaled or unscaled.
            t1 = powS.apply(a, 1, -s);
            if (t1 < Double.MIN_VALUE) {
                // If t0 was sub-normal then return t0 as t1 cannot be added.
                // If t0 was normal check if a sub-normal term can be added
                // otherwise s is very large and the result is t0 to 1 ULP.
                if (scaledSum) {
                    // Ideally return originalT0. But we must support result scaling.
                    if (resultScale == 1) {
                        return originalT0;
                    }
                    // Here we must return the best we can compute for t0 * 2^106
                    return scaledSum ? t0 * DOWNSCALE_916 : originalT0 * UPSCALE_106;
                }
                t0 *= UPSCALE_1022;
                scaledSum = true;
                // Redefine the series power function
                PowOperator originalPowSeries = powSeries;
                powS = (x, y, z) -> {
                    final double p = originalPowSeries.apply(x, y, z * 0.5);
                    return p * UPSCALE_1022 * p;
                };
                t1 = powS.apply(a, 1, -s);
            }
        } else {
            t1 = powSeries.apply(a, 1, -s);
        }

//        // TODO - Remove original
//        // If t0 is sub-normal the summation S is of sub-normal terms.
//        // Avoid by using: scale * x^y = x^(y/2) * scale * x^(y/2).
//        final boolean scaledSum = c.isSet(Context.SUB_NORMAL_SUM) && t0 < Double.MIN_NORMAL;
//        if (scaledSum) {
//            // Re-compute t0 using scaling
//            t0 = Math.pow(a, -s * 0.5);
//            t0 = t0 * UPSCALE_1022 * t0;
//            if (t0 == 0) {
//                // For any s > 1 the I term is not 2^1022 larger than a^-s.
//                // No need to compute 2^-1022 * powS.apply(a, N, 1 - s) / (s - 1).
//                return 0;
//            }
//            // Only the power function for the series is scaled
//            PowOperator originalPowSeries = powSeries;
//            powS = (x, y, z) -> {
//                final double p = originalPowSeries.apply(x, y, z * 0.5);
//                return p * UPSCALE_1022 * p;
//            };
//        } else {
//            powS = powSeries;
//        }

        // Here we compute the entire evaluation of the series N terms.
        // If t1 cannot be added to t0 reduce n.
        if (t1 < t0 * c.getEps()) {
            // This still requires evaluation of the I term which may dominate with small s
            n = 2;
        }
        double p = powS.apply(a, n, -s);

        // We always use a DD sum. If not using extended precision
        // we add in double precision and create a new DD.
        BiFunction<DD, Double, DD> add = c.isSet(Context.EXT_SUM) ?
            DD::add :
            (x, y) -> DD.of(x.hi() + y);

        // Continue optional sub-normal support
        if (c.isSet(Context.SUB_NORMAL_SUM)) {
            if (p < Double.MIN_NORMAL && !scaledSum) {
                // The final term is sub-normal.
                // If this occurs when the scaled sum has been activated then:
                // - the series is very strongly converging and the result is t0
                // - the zeta result is dominated by the I term for small s
                //   after summing many sub-normal terms and the full
                //   evaluation is required.
                // Activate scaling when possible. Very strongly converging series use n=2.
                if (t0 < 1) {
                    t0 *= UPSCALE_1022;
                    t1 *= UPSCALE_1022;
                    scaledSum = true;
                    // Redefine the series power function
                    PowOperator originalPowSeries = powSeries;
                    powS = (x, y, z) -> {
                        final double pp = originalPowSeries.apply(x, y, z * 0.5);
                        return pp * UPSCALE_1022 * pp;
                    };
                    p = powS.apply(a, n, -s);
                }
            }
        }

        // Initialise sum with the first tail term
        DD sum = DD.of(0.5 * p);
        // S : k in [0, n-1]
        // Note k=0,1 already computed as t0,t1
        for (int k = n - 1; k > 1; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum = add.apply(sum, powS.apply(a, k, -s));
        }
        sum = add.apply(sum, t1);

        // I: (a+n)^(1-s) / (s-1) = (a+n)^-s * (a+n) / (s-1)
        // Note: This term represents the sum of the remaining terms to infinity.
        // For S ~ I and I to be sub-normal requires a sum of 9 sub-normal terms.
        // Only compute scaled when the sum is scaled: a^-s is sub-normal and a < 2^53
        // (asymptotic case) hence s > 2 and division by (s-1) is safe.
        double ti;
        if (scaledSum) {
            // Determine if (a+n)^(1-s) is sub-normal and compute in parts.
            if (p * (a + n) < 1) {
                // Note: This must use the pow for the main terms and apply scaling
                // as powS may be standard precision
                ti = pow.apply(a, n, (1 - s) * 0.5);
                ti = ti * UPSCALE_1022 * ti;
            } else {
                ti = UPSCALE_1022 * pow.apply(a, n, 1 - s);
            }
            ti /= s - 1;
            // Add in magnitude order. When a in [0, 1] it may be the dominant term.
            // Note: This makes a measurable difference to the RMS error when the summation
            // is in single precision. When in double-double precision the difference
            // cannot be measured with the test resources.
            if (t0 > ti) {
                sum = add.apply(sum, ti);
                sum = add.apply(sum, t0);
            } else {
                sum = add.apply(sum, t0);
                sum = add.apply(sum, ti);
            }
            // Partial downscale of the sum. The tail can use normal numbers.
            sum = DD.ofSum(sum.hi() * DOWNSCALE_916, sum.lo() * DOWNSCALE_916);
        } else {
            ti = pow.apply(a, n, 1 - s) / (s - 1);
            // Add in magnitude order
            if (t0 > ti) {
                sum = add.apply(sum, ti);
                sum = add.apply(sum, t0);
            } else {
                sum = add.apply(sum, t0);
                sum = add.apply(sum, ti);
            }
            // Upscale the sum so the tail can use normal numbers
            sum = DD.ofSum(sum.hi() * UPSCALE_106, sum.lo() * UPSCALE_106);
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s) using division by (a+N)^2.
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        //
        // If a or s are large then (a+n)^-(2k-1+s) will underflow.
        // We can downscale the rising factorial (always above 1)
        // and upscale p (always below 1) by the same amount.
        // Scaling could dynamically use the exponent of f and p. Here we use
        // fixed scaling which can be pre-applied to a possible sub-normal
        // power term. This allows a few terms of the tail to be computed
        // when the zeta result or terms are sub-normal and a significant
        // improvement in accuracy for sub-normal results.
        //
        // The loop will exit immediately if p = (a+n)^-(1-s) = 0:
        // Extend support for large s by computing sub-normal (a+n)^-(1+s)
        // in parts with fixed scaling:
        //   scale * (a+n)^-(1+s) = (a+n)^-(1+s)/2 * scale * (a+n)^-(1+s)/2
        // The limit on (1+s):
        //   9^-(1+s) = 2^-(1075+1022) : 1+s = log(2)/log(9) * (1075+1022) = 662
        //
        // Rising factorial term : (s)_{2k-1}
        // The maximum s is < 662.
        // (662)_{2k-1}; k=53 = Pochammer(662, 105) = 3.96e299
        // The rising factorial will not overflow for k <= 53 before (a+n)^-(1+s) is
        // zero. This is within the length of table F.
        // For larger s where p=0 the loop exits immediately.
        //
        // Note that downscaling of f does not match upscaling of p. The mismatch of 2^106
        // allows the tail to sum normal terms and converge on a normal number.
        // Due to a prior check on large a^-s the zeta result cannot be above 2^106 and
        // overflow cannot occur.
        double f = s * DOWNSCALE_916;
        // 2k - 1
        double k2 = 1;

        // Use to divide by (a+n)^2
        final double apn = DD.ofSum(a, n).square().hi();

        // Scaled power term (a+n)^-(1+s).
        // The tail is expected to be a small contribution to the result.
        // This term is not required to high accuracy so compute as
        // (a+n)^-s / (a+n).
        if (scaledSum) {
            p /= a + n;
        } else {
            // Power term was not scaled.
            // Note if p < MIN_NORMAL this is inexact. The smallest s for this to occur
            // is (2^-53)^-s < 2^-1022 : s = 1022/53 ~ 19.28. In this case the tail
            // is not required to high accuracy.
            p = UPSCALE_1022 * p / (a + n);
        }

        // The result is 2^106 too high.
        // Remove this scaling unless the scaled result has been requested.
        // The result scale is only changed when tail scaling is active.
        double downscale = resultScale == 1 ? DOWNSCALE_106 : 1;
        if (!c.isSet(Context.TAIL_SCALING)) {
            // Currently: p * 2^1022; f * 2^-916; sum * 2^106
            // The terms sum is 2^106 too large. This allows the tail sum to use
            // normal numbers and converge on a normal result.
            // Disable this scaling. The function will readily underflow when
            // the tail power term is sub-normal and there can be significant
            // loss of precision if the tail was required.
            p *= 0x1p-1022;
            f = s;
            sum = DD.ofSum(sum.hi() * DOWNSCALE_106, sum.lo() * DOWNSCALE_106);
            downscale = 1;
        }

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Note: an extended precision sum here has no effect on the final result.
        final double stop = sum.hi() * c.getEps();
        double tsum = 0;
        int i;
        for (i = 0; i < c.getM(); i++) {
            final double t = f * p * F[i];
            tsum += t;
            if (Math.abs(t) <= stop) {
                break;
            }
            // Update (a+n)^-(2k-1+s)
            p /= apn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
        }
        // Used to histogram convergence when testing
        lastM = i;
        M[i]++;
        return add.apply(sum, tsum).hi() * downscale;
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}
     * when {@code a} is negative. Uses {@link BigDecimal} arithmetic.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be negative.
     *
     * @param s Argument {@code s > 1} and integer
     * @param a Argument {@code a < 0}
     * @return zeta(s, a)
     */
    private static double zetaNegativeBD(int s, double a) {
        // a < 0 (non-integer) and s is a positive integer.
        // If s is odd then the negative series sum will be negative
        // and the addition of zeta(s, x > 0) has cancellation.
        // This is largest when a is close to half-integer.

        // Check case of total cancellation.
        final boolean odd = (s & 1) == 1;
        final double ca = Math.ceil(a);
        if (ca == a) {
            // The term 0^-s is infinity
            return Double.POSITIVE_INFINITY;
        }
        final double x = a - ca;
        // Intentional float comparison
        if (odd && x == -0.5) {
            // Use extended precision but evaluated with precision for a double result
            return zeta(s, BigDecimal.ONE.subtract(new BigDecimal(a)), null,
                Context.BD_DOUBLE).doubleValue();
        }

        // Handle cancellation as x -> 0.5
        // Note: 0.5^-1024 overflows.
        // Limit of [nextDown(0.5)^-s - nextUp(0.5)^-s] may have terms above 2^1024.
        // 0.5 +/- 2^-54 (requires extended precision as ulp(0.5) is 2^-53)
        // The largest odd s where the difference is finite:
        // var mc = MathContext.DECIMAL128
        // var a = new BigDecimal(Math.nextDown(0.5))
        // var b = BigDecimal.ONE.subtract(a)
        // var s = -1025
        // while (Double.isFinite(a.pow(s, mc).subtract(b.pow(s, mc), mc).doubleValue())) { s -= 2; }
        // s = -1067 : diff = 3.75e308
        if (s >= 1067) {
            // Use the dominant term using closest to zero
            return odd && x > -0.5 ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        }

        // Odd computation requires twice the precision of a double (17 digits).
        // Even requires some extra.
        final Context c = odd ? Context.BD_DOUBLE_DOUBLE : Context.BD_DOUBLE;

        // Compute the two terms either side of zero:
        // -1 < xn < 0 < xn + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        final BigDecimal xn = new BigDecimal(x);
        final BigDecimal xp = BigDecimal.ONE.add(xn);
        final MathContext mc = c.getMathContext();
        final BigDecimal pn = xn.pow(-s, mc);
        final BigDecimal pp = xp.pow(-s, mc);

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2.
        // Exit early if remaining terms cannot be added.
        final double d = pn.add(pp, mc).doubleValue();
        if (Math.abs(d) > 0x1p106) {
            // Adding to either side will not change a double result
            return d;
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) + (-1)^s * [ zeta(s, -x) - zeta(s, 1 - a) ]
        // z + za + zb

        BigDecimal z = zeta(s, xp, pp, c);
        BigDecimal za;
        BigDecimal zb;
        // A single call to zeta uses many pow operations;
        // use a direct sum when zeta will use more.
        if (x - a > 2 * (c.getN() + 2)) {
            // Note: The difference (za - zb) incurs cancellation.
            // This should not be an issue as x in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 1)  = 1.6449
            // zeta(2, 2)  = 0.6449
            // zeta(2, 30) = 0.0339
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            za = zeta(s, xn.negate(), pn.abs(), c);
            zb = zeta(s, BigDecimal.ONE.subtract(new BigDecimal(a)), null, c);
            // Both terms are positive. Correct the sign for final addition.
            if (odd) {
                za = za.negate();
            } else {
                zb = zb.negate();
            }
        } else {
            // Sum terms in ascending order of magnitude.
            // Use a double to track the iterations, and mirror with a BigDecimal.
            za = pn;
            zb = BigDecimal.ZERO;
            BigDecimal ba = new BigDecimal(a);
            for (double aa = a; aa < x; aa += 1.0) {
                zb = zb.add(ba.pow(-s, mc), mc);
                ba = ba.add(BigDecimal.ONE);
            }
        }

        // Sum in magnitude order. Use the scale for a fast comparison of magnitude
        // as all results have the same precision: base 10 exponent = precision - scale - 1
        // Smaller scale is a bigger value.
        if (za.scale() < z.scale()) {
            final BigDecimal tmp = za;
            za = z;
            z = tmp;
        }
        // za < z
        if (zb.scale() < z.scale()) {
            final BigDecimal tmp = zb;
            zb = z;
            z = tmp;
        }
        // za,zb < z
        return za.add(zb, mc).add(z, mc).doubleValue();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @param a0 {@code a^-s} (can be null)
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    private static BigDecimal zeta(int s, BigDecimal a, BigDecimal a0, Context c) {
        final int n = c.getN();
        final MathContext mc = c.getMathContext();
        BigDecimal apn = a.add(BigDecimal.valueOf(n));
        BigDecimal p = apn.pow(-s, mc);

        // Initialise sum with the first tail term
        BigDecimal sum = p.multiply(new BigDecimal(0.5));
        // S : k in [0, n-1]
        for (int k = n - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum = sum.add(a.add(BigDecimal.valueOf(k)).pow(-s, mc), mc);
        }
        // First term may be provided
        final BigDecimal t0 = a0 == null ? a.pow(-s, mc) : a0;

        // I : (a+p)^(1-s) / (s-1)
        final BigDecimal ti = apn.pow(1 - s, mc).divide(BigDecimal.valueOf(s - 1), mc);

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (t0.compareTo(ti) > 0) {
            sum = sum.add(ti, mc).add(t0, mc);
        } else {
            sum = sum.add(t0, mc).add(ti, mc);
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.

        // *** This must work for all s where (a+n)^-s is not an underflow ***
        // The function testTailConvergence examines the parameter space in DD arithmetic
        // and verifies convergence is possible for the provided n.

        // Rising factorial term : (s)_{2k-1}
        BigDecimal f = BigDecimal.valueOf(s);

        // Alternative implementations for (a+n)^-(2k-1+s).

        // Initialise (a+n)^-(2k-1+s) to (a+n)^-(1+s)
        p = apn.pow(-(1 + s), mc);
        // Compute using the power function
        final boolean usePow = c.isSet(Context.TAIL_POW);
        if (!usePow) {
            // Divide by (a+n)^2 using multiplication
            apn = apn.pow(-2, mc);
        }

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        BigDecimal tsum = BigDecimal.ZERO;
        final int stop = sum.scale() + mc.getPrecision();
        int i;
        for (i = 0; i < c.getM(); i++) {
            final BigDecimal t = f.multiply(p, mc).multiply(FBD[i], mc);
            tsum = tsum.add(t, mc);
            // Bigger scale is a smaller value
            if (t.scale() >= stop) {
                break;
            }
            // p = (a+n)^-(2k-1+s)
            if (usePow) {
                p = apn.pow(-(2 * i + 3 + s), mc);
            } else {
                p = p.multiply(apn, mc);
            }
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            f = f.multiply(BigDecimal.valueOf((s + (2L * i) + 1) * (s + (2L * i) + 2)), mc);
        }
        // Used to histogram convergence when testing
        M[i]++;
        return sum.add(tsum, mc);
    }

    /**
     * Calculates the sum of terms of the power series.
     *
     * <pre>
     *      b     1
     *   sum     ---
     *      k=a  k^s
     * </pre>
     *
     * <p>Assumes {@code a} and {@code b} are negative and separated by an integer
     * distance; and {@code exponent >= 2} and integer.
     *
     * <p>Large ranges may be evaluated using a difference of zeta functions.
     *
     * <p>This method is a reproduction of the logic in {@link #zetaNegativeBD(int, double)}
     * so the magnitude of terms that cancel can be computed in the search to find
     * the roots (see {@link #testDataZetaRoots(int, int, int, int)}).
     *
     * @param a First term in the series to calculate (negative non-integer).
     * @param x Last term in the series inclusive (negative non-integer); result provided.
     * @param s Exponent (positive integer).
     * @param c Evaluation context.
     * @return the sum
     */
    private static BigDecimal negativeSeriesSum(double a, double x, int s, Context c) {
        final MathContext mc = c.getMathContext();
        final BigDecimal pn = new BigDecimal(x).pow(-s, mc);
        // A single call to zeta uses many pow operations;
        // use a direct sum when zeta will use more.
        if (x - a > 2 * (c.getN() + 2)) {
            // Note: The difference incurs cancellation.
            // This should not be an issue as x in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 1)  = 1.6449
            // zeta(2, 2)  = 0.6449
            // zeta(2, 30) = 0.0339
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            final BigDecimal za = zeta(s, new BigDecimal(-x), pn.abs(), c);
            final BigDecimal zb = zeta(s, BigDecimal.ONE.subtract(new BigDecimal(a)), null, c);
            final BigDecimal r = za.subtract(zb, mc);
            return (s & 1) == 1 ? r.negate() : r;
        }

        // Sum terms in ascending order of magnitude.
        // Use a double to track the iterations, and mirror with a BigDecimal.
        BigDecimal zb = BigDecimal.ZERO;
        BigDecimal ba = new BigDecimal(a);
        for (double aa = a; aa < x; aa += 1.0) {
            zb = zb.add(ba.pow(-s, mc), mc);
            ba = ba.add(BigDecimal.ONE);
        }
        return zb.add(pn, mc);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}
     * when {@code a} is negative. Uses double-double ({@link DD}) arithmetic.
     *
     * <p>This method is very similar to the BigDecimal variant
     * {@link #zetaNegativeBD(int, double)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be negative.
     *
     * @param s Argument {@code s > 1} and integer
     * @param a Argument {@code a < 0}
     * @return zeta(s, a)
     */
    private static double zetaNegativeDD(int s, double a) {
        // a < 0 (non-integer) and s is a positive integer.
        // If s is odd then the negative series sum will be negative
        // and the addition of zeta(s, x > 0) has cancellation.
        // This is largest when a is close to half-integer.

        // Check case of total cancellation.
        final boolean odd = (s & 1) == 1;
        final double ca = Math.ceil(a);
        if (ca == a) {
            // The term 0^-s is infinity
            return Double.POSITIVE_INFINITY;
        }

        // Exponent for fractional representation
        final int[] exp = {0};

        final double x = a - ca;
        // Intentional float comparison
        if (odd && x == -0.5) {
            // Use extended precision but evaluated with precision for a double result.
            final DD r = zeta(s, DD.ONE.subtract(a), null, exp, Context.DD_DOUBLE);
            return Math.scalb(r.hi(), exp[0]);
        }

        // Handle cancellation as x -> 0.5
        // Note: 0.5^-1024 overflows.
        // Limit of [nextDown(0.5)^-s - nextUp(0.5)^-s] may have terms above 2^1024.
        // 0.5 +/- 2^-54 (requires extended precision as ulp(0.5) is 2^-53)
        // The largest odd s where the difference is finite:
        // var mc = MathContext.DECIMAL128
        // var a = new BigDecimal(Math.nextDown(0.5))
        // var b = BigDecimal.ONE.subtract(a)
        // var s = -1025
        // while (Double.isFinite(a.pow(s, mc).subtract(b.pow(s, mc), mc).doubleValue())) { s -= 2; }
        // s = -1067 : diff = 3.75e308
        if (s >= 1067) {
            // Use the dominant term using closest to zero
            return odd && x > -0.5 ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        }

        // Compute the two terms either side of zero:
        // -1 < x < 0 < x + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        final DD xn = DD.of(x);
        final DD xp = DD.ONE.add(x);
        final long[] expn = {0};
        final long[] expp = {0};
        final DD pn = DDMath.pow(xn, -s, expn);
        final DD pp = DDMath.pow(xp, -s, expp);
        final long diff = expp[0] - expn[0];
        if (Math.abs(diff) > 106) {
            // Cannot add these terms given the largest is [0.5, 1.0) with ulp 2^-106.
            // All other individual terms are smaller so further DD computation is not possible.
            return diff > 0 ?
                Math.scalb(pp.hi(), intExponent(expp[0])) :
                Math.scalb(pn.hi(), intExponent(expn[0]));
        }

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2; any individual term x^-s < 1.
        // Exit early if remaining terms cannot be added.
        // The result can be finite even if the terms are infinite. However
        // we do not support further computation from infinite terms so check
        // if anything can be added to these terms.
        if (Math.max(expn[0], expp[0]) > 106) {
            // Limit of double-double arithmetic
            // Add a scaled term avoiding overflow by re-scaling after the sum
            return Math.scalb(
                HurwitzZetaTest.accurateAdd(pn, pp.scalb((int) diff)).hi(),
                intExponent(expn[0])
            );
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) + (-1)^s * [ zeta(s, -x) - zeta(s, 1 - a) ]
        // z + za + zb

        // Worst case single term cancellation:
        // priority = exponent(max(a, b)) - exponent(a - b)
        // Cancellation is worst when s is small as the two terms are closer.
        // pow(0.5 + 2^-54, -3) - pow(0.5 - 2^-54, -3)  // 0.5+2^-54 requires extended precision
        // exponent(8.0) - exponent(5.33E-15) = 3 - -48 = 51 bits

        // We have the two largest power terms for each side.
        // The remaining terms are increasingly smaller. Computing with
        // double-double (DD) precision for the zeta evaluations should handle cancellation.

        // Evaluate zeta with extra precision
        final Context c = odd ? Context.DD_DOUBLE_DOUBLE : Context.DD_DOUBLE;

        // Use scaled terms in the DD implementation.
        // This allows correct convergence of the tail when sub-normal.

        exp[0] = intExponent(expp[0]);
        DD z = zeta(s, xp, pp, exp, c);
        int ez = exp[0];

        DD za;
        DD zb;
        int eza;
        int ezb;
        // A single call to zeta uses many pow operations;
        // use a direct sum when zeta will use more.
        // Prepare the exponent for the term we have.
        exp[0] = intExponent(expn[0]);
        if (x - a > 2 * (c.getN() + 2)) {
            // Note: The difference (za - zb) incurs cancellation.
            // This should not be an issue as x in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 1)  = 1.6449
            // zeta(2, 2)  = 0.6449
            // zeta(2, 30) = 0.0339
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            za = zeta(s, xn.negate(), pn.abs(), exp, c);
            eza = exp[0];
            zb = zeta(s, DD.ONE.subtract(a), null, exp, c);
            ezb = exp[0];
            // Both terms are positive. Correct the sign for final addition.
            if (odd) {
                za = za.negate();
            } else {
                zb = zb.negate();
            }
        } else {
            // Sum terms in ascending order of magnitude
            // Using a double to track the iterations is fine as (a+n) is exact until > x.
            za = pn;
            eza = exp[0];
            zb = DD.ZERO;
            ezb = 0;
            // Power function for the series (allows switching to a faster pow)
            final ScaledPowOperator powS = c.isSet(Context.EXT_POW_SERIES) ? DDMath::pow : DD::pow;
            final BiFunction<DD, DD, DD> addS = c.isSet(Context.AA_SERIES) ? HurwitzZetaTest::accurateAdd : DD::add;
            final long[] e = {0};
            for (double aa = a; aa < x; aa += 1.0) {
                final DD r = powS.apply(DD.of(aa), -s, e);
                zb = add(zb, ezb, r, intExponent(e[0]), exp, addS);
                ezb = exp[0];
            }
            // s >= 2. For any -1 < x < 0 : |x|^-s > zeta(s, |x|+1)
            // No sort required for the summation by magnitude
        }

        // Sum in magnitude order
        final BiFunction<DD, DD, DD> add = c.isSet(Context.AA_TERMS) ? HurwitzZetaTest::accurateAdd : DD::add;
        // zb < za : check z against za
        if (ez > eza) {
            // za,zb < z
            za = add(za, eza, zb, ezb, exp, add);
            eza = exp[0];
        } else {
            // z,zb < za
            z = add(z, ez, zb, ezb, exp, add);
            ez = exp[0];
        }
        // Note: DD result could be used in evaluating the polygamma function.
        // The support domain would have to be expanded as a sub-normal zeta can
        // be multiplied by a large factorial term:
        // psi(n, z) = (-1)^(n+1) n! zeta(n+1, z)
        z = add(za, eza, z, ez, exp, add);
        return Math.scalb(z.hi(), exp[0]);
    }

    /**
     * Helper function to add the fractional representation of two DD numbers.
     *
     * <p>Note a DD has a floating mantissa of 106+ bits. This method allows adding
     * a small term to an effectively 107-bit larger term.
     * See {@link #testFractionalAddAssumptions()}.
     *
     * @param a Addend a.
     * @param ae Exponent of a.
     * @param b Addend b.
     * @param be Exponent of b.
     * @param e Exponent of the result.
     * @param add Addition operator.
     * @return the fraction part of the result
     */
    private static DD add(DD a, int ae, DD b, int be, int[] e, BiFunction<DD, DD, DD> add) {
        // Add if the overlap is within 107-bits.
        final long diff = (long) ae - be;
        if (Math.abs(diff) > 107) {
            // Special case where the exponent of zero is 0.
            // This may be larger than the exponent of the other term.
            if (ae > be) {
                if (a.isZero()) {
                    e[0] = be;
                    return b;
                }
                e[0] = ae;
                return a;
            } else {
                if (b.isZero()) {
                    e[0] = ae;
                    return a;
                }
                e[0] = be;
                return b;
            }
        }
        final DD r = add.apply(b, a.scalb((int) diff)).frexp(e);
        e[0] = addExponent(e[0], be);
        return r;
    }

    /**
     * Helper function to add the two exponents and maintain the [min, max]
     * range of an integer.
     *
     * @param a Exponent.
     * @param b Exponent.
     * @return (a + b) within the range of an integer
     * @see Integer#MIN_VALUE
     * @see Integer#MAX_VALUE
     */
    private static int addExponent(int a, int b) {
        // Use the tested version
        return SpecialMath.addAsInt(a, b);
    }

    /**
     * Clip the long exponent to an integer range.
     *
     * <p>Note that use of a direct cast to {@code int} is a 32-bit truncation. If all
     * lower bits are zero but the high 32-bits contain non-zero bits then a cast is
     * incorrect. This method ensures {@code long} exponents generated from large
     * power computations are clipped to the range {@code [-(2^31), 2^31)}. This is
     * sufficient to detect infinite or zero values, and correctly scale results
     * using {@link Math#scalb(double, int)} or {@link DD#scalb(int)}.
     *
     * @param e Exponent.
     * @return the exponent as the closest integer
     */
    private static int intExponent(long e) {
        // Use the tested version
        return SpecialMath.toInt(e);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * <p>Support for sub-normal double values is provided by returning the
     * result using a fractional representation in {@code [0.5, 1)} with
     * a base-2 exponent.
     *
     * <p>The first term can be provided. This must be the fractional representation
     * with the exponent in {@code exp}.
     *
     * @param s Argument {@code s > 1}, expected {@code s < 1024}
     * @param a Argument {@code a > 0}
     * @param a0 {@code a^-s} (can be null; if non-null it is assumed to be finite)
     * @param exp base-2 exponent of the result.
     * @param c Evaluation context.
     * @return zeta(s, a) fractional part in [0.5, 1)
     * @see DD#frexp(int[])
     */
    private static DD zeta(int s, DD a, DD a0, int[] exp, Context c) {
        // Note:
        // The power function uses a long exponent to support the full range
        // of a double to an integer power: (2^1024)^(2^31) = 2^(2^10 * 2^31) = 2^(2^41)
        // Here we can clip the long to an integer. The final result is either
        // infinite or zero if the exponent is outside [-1075, 1024) so special
        // routines are used to clip the long exponents to integer and add them
        // without exceeding the range of an integer.
        // Only the first term is above 1; all other terms are < 1 and can be added
        // exactly even if sub-normal.
        // Since scaling of the power result would be required to add in standard representation
        // the cost of this is an additional re-scaling (scalb) after each addition. This
        // is negligible compared to the cost of a power function evaluation.

        final long[] e = {0};

        // Power function and add function (allows switching to a faster method for the series)
        final ScaledPowOperator pow = c.isSet(Context.EXT_POW_TERMS) ? DDMath::pow : DD::pow;
        final ScaledPowOperator powS = c.isSet(Context.EXT_POW_SERIES) ? DDMath::pow : DD::pow;
        final BiFunction<DD, DD, DD> add = c.isSet(Context.AA_TERMS) ? HurwitzZetaTest::accurateAdd : DD::add;
        final BiFunction<DD, DD, DD> addS = c.isSet(Context.AA_SERIES) ? HurwitzZetaTest::accurateAdd : DD::add;

        // First term may be provided
        final DD t0;
        final int et0;
        if (a0 == null) {
            t0 = pow.apply(a, -s, e);
            et0 = intExponent(e[0]);
        } else {
            t0 = a0;
            et0 = exp[0];
        }

        // Asymptotic Behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0, the I term and the first term of T.
        if (a.hi() > 0x1p105) {
            // a^(1-s) / (s - 1) + a^-s / 2
            // Note: If a^-s is sub-normal then a^(1-s)=(a^-s * a) may not be
            // sub-normal when s is small as a is large.
            DD r = pow.apply(a, 1 - s, e);
            long er = e[0];
            r = r.divide(s - 1).frexp(exp);
            er += exp[0];
            // Add 0.5 * t0 so adjust the exponent
            return add(r, intExponent(er), t0, addExponent(et0, -1), exp, add);
        }

        // Can overflow if 0 < a < 1.
        if (et0 > 106) {
            // Nothing can be added in double-double precision
            exp[0] = et0;
            return t0;
        }
        // Now any (a+n)^-s cannot overflow and the sum cannot overflow.
        // However the terms may be sub-normal and we can add them using a fractional
        // representation with a separate base 2 exponent.

        final int n = c.getN();
        DD apn = a.add(n);

        // Initialise sum with the first tail term: 0.5 * (a+n)^-s
        DD sum = powS.apply(apn, -s, e);
        int es = addExponent(intExponent(e[0]), -1);
        // S : k in [0, n-1]
        // Note k=0 already computed as t0
        for (int k = n - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            final DD r = powS.apply(a.add(k), -s, e);
            sum = add(sum, es, r, intExponent(e[0]), exp, addS);
            es = exp[0];
        }

        // I : (a+p)^(1-s) / (s-1)
        DD ti = pow.apply(apn, 1 - s, e);
        long eti = e[0];
        ti = ti.divide(s - 1).frexp(exp);
        eti += exp[0];

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (et0 > eti) {
            sum = add(sum, es, ti, intExponent(eti), exp, add);
            sum = add(sum, exp[0], t0, et0, exp, add);
        } else {
            sum = add(sum, es, t0, et0, exp, add);
            sum = add(sum, exp[0], ti, intExponent(eti), exp, add);
        }
        es = exp[0];

        // If the sum of terms is far below min double value then stop.
        // The smallest normalised exponent is -1023 - 52 = -1075.
        // Make sure the 106-bit double-double cannot be added to Double.MIN_VALUE.
        // This is disabled to allow computation on extremely small values for testing.
        //if (es < -1075 - 106) {
        //    exp[0] = 0;
        //    return DD.ZERO;
        //}

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        // Terms cannot overflow using the fractional representation. This also
        // allows the tail computation to continue further when an individual
        // part of each term computation would underflow.
        // Note:
        // This has a drawback when s > (a+n) and the pochammer function (s)_{2k-1} grows
        // faster than the term (a+n)^-(2k-1-s) reduces. For large s and/or a the
        // power term would normally underflow and the tail is ignored. Here the
        // convergence must be possible for all expected values of s and a.

        // *** This must work for all s where (a+n)^-s is not an underflow ***
        // The function testTailConvergence examines the parameter space and verifies
        // convergence is possible for the provided n and expected range of s.

        // Rising factorial term : (s)_{2k-1}
        DD f = DD.of(s).frexp(exp);
        long ef = exp[0];

        // Alternative implementations for (a+n)^-(2k-1+s).

        // Initialise (a+n)^-(2k-1+s) to (a+n)^-(1+s)
        DD p = pow.apply(apn, -1 - s, e);
        long ep = e[0];
        // Compute using the power function
        final boolean usePow = c.isSet(Context.TAIL_POW);
        long eapn = 0;
        if (!usePow) {
            // Divide by (a+n)^2 using multiplication
            apn = pow.apply(apn, -2, e);
            eapn = e[0];
        }

        // Set the scale for the tail series

        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Use the tail epsilon to obtain a minimum exponent for terms.
        final long stop = (long) es + Math.getExponent(c.getEps());
        // Initialise with the first term.
        // Avoids addition of fractional representation to zero (which has a zero exponent).
        DD tsum = f.multiply(p).multiply(FDD[0]).frexp(exp);
        long et = ef + ep + exp[0];
        int etsum = intExponent(et);
        int i;
        for (i = 1; i < c.getM(); i++) {
            if (et < stop) {
                break;
            }
            // p = (a+n)^-(2k-1+s)
            if (usePow) {
                p = powS.apply(apn, -(2 * i + 1 + s), e);
                ep = e[0];
            } else {
                p = p.multiply(apn).frexp(exp);
                ep += exp[0] + eapn;
            }
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            f = f.multiply((s + (2L * i) - 1) * (s + (2L * i))).frexp(exp);
            ef += exp[0];
            // FDD is finite. f*p in [0.25, 2).
            final DD t = f.multiply(p).multiply(FDD[i]).frexp(exp);
            et = ef + ep + exp[0];
            tsum = add(tsum, etsum, t, intExponent(et), exp, addS);
            etsum = exp[0];
        }
        // Used to histogram convergence when testing
        M[i]++;
        return add(sum, es, tsum, etsum, exp, add);
    }

    /**
     * Helper function to compute a BigDecimal result using the DD zeta implementation.
     *
     * @param s Argument {@code s > 1}, expected {@code s < 1024}
     * @param a Argument {@code a > 0}
     * @param c Evaluation context.
     */
    private static BigDecimal zeta(int s, DD a, Context c) {
        final int[] exp = {0};
        final DD r = HurwitzZetaTest.zeta(s, a, null, exp, c);
        return toBigDecimal(r, exp[0]);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}
     * when {@code a} is negative. Uses double arithmetic.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be negative.
     *
     * @param s Argument {@code s > 1} and integer
     * @param a Argument {@code a < 0}
     * @param firstTerm if true use double-double for the first terms
     * @return zeta(s, a)
     */
    private static double zetaNegative(int s, double a, boolean firstTerm) {
        // a < 0 (non-integer) and s is a positive integer.
        // If s is odd then the negative series sum will be negative
        // and the addition of zeta(s, x > 0) has cancellation.
        // This is largest when a is close to half-integer.

        // Check case of total cancellation.
        final boolean odd = (s & 1) == 1;
        final double ca = Math.ceil(a);
        if (ca == a) {
            // The term 0^-s is infinity
            return Double.POSITIVE_INFINITY;
        }

        // Evaluation context for zeta. Use the default optimised for double precision.
        final Context c = Context.DOUBLE;
        final PowOperator pow = Context::extPow;

        final double x = a - ca;
        // Intentional float comparison
        if (odd && x == -0.5) {
            return zeta(s, 1 - a, pow.apply(-a, 1, -s), c);
        }

        // Compute the dominant term using closest to zero: x or 1+x.
        // This optionally is performed using double-double precision.
        if (firstTerm) {
            // Copy from DD implementation
            if (s >= 1067) {
                // Use the dominant term using closest to zero
                return odd && x > -0.5 ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            }
            // Compute the two terms either side of zero:
            // -1 < xn < 0 < xn + 1 < 1
            // These are the largest terms and contain most of the error of the function.
            final DD xn = DD.of(x);
            final DD xp = DD.ONE.add(x);
            final long[] expn = {0};
            final long[] expp = {0};
            DD pn = DDMath.pow(xn, -s, expn);
            DD pp = DDMath.pow(xp, -s, expp);
            final long diff = expp[0] - expn[0];
            if (Math.abs(diff) > 106) {
                // Cannot add these terms given the largest is [0.5, 1.0) with ulp 2^-106.
                // All other individual terms are smaller so further DD computation is not possible.
                return diff > 0 ?
                    pp.scalb(intExponent(expp[0])).hi() :
                    pn.scalb(intExponent(expn[0])).hi();
            }
            // Add the smallest to the largest avoiding overflow by re-scaling after the sum
            final DD sum = pn.add(pp.scalb((int) diff)).scalb(intExponent(expn[0]));

            // Rescale terms
            pp = pp.scalb(intExponent(expp[0]));
            pn = pn.scalb(intExponent(expn[0]));

            // Here the remaining series above and below zero are effectively both zeta
            // evaluations with zeta(s >= 2, a > 1). This is always < 2; any individual term x^-s < 1.
            // Exit early if remaining terms cannot be added.
            if (Math.abs(sum.hi()) > 0x1p53) {
                // Limit of double arithmetic
                return sum.hi();
            }

            // The following is an adaption of the pure double implementation
            // to compute the zeta terms without the already known first terms in DD precision.

            // x = a - ceil(a) : x in -(1, 0)
            // (x+1)^-s + zeta(s, x + 2) + (-1)^s * [ -x^-s + zeta(s, 1-x) - zeta(s, 1 - a) ]
            // z + za + zb

            // We have the two largest power terms for each side.
            // The remaining terms are increasing smaller.

            DD z = pp.add(zeta(s, x + 2, pow.apply(x, 2, -s), c));

            DD za;
            double zb;
            // A single call to zeta uses many pow operations;
            // use a direct sum when zeta will use more.
            if (x - a > 2 * (c.getN() + 2)) {
                // Note: The difference (za - zb) incurs cancellation.
                // This should not be an issue as x in -(2, 1)
                // and the series is strongly converging, e.g.
                // zeta(2, 2)  = 0.6449
                // zeta(2, 30) = 0.0339
                // Significant cancellation (leading digits the same) is not possible.
                // Take care to change the sign of a provided result for the zeta method.
                za = pn.abs().add(zeta(s, 1 - x, pow.apply(-x, 1, -s), c));
                zb = zeta(s, 1 - a, pow.apply(-a, 1, -s), c);
                // Both terms are positive. Correct the sign for final addition.
                if (odd) {
                    za = za.negate();
                } else {
                    zb = -zb;
                }
            } else {
                // Sum terms in ascending order of magnitude
                // Using a double to track the iterations is fine as (a+n) is exact until > x.
                za = pn;
                zb = 0;
                // Instead of (a+k) use a+k = x+i with i in [a-x, 0) to use the extended precision power.
                for (int i = (int) (a - x); i < 0; i++) {
                    zb += pow.apply(x, i, -s);
                }
            }

            // Sum in magnitude order. Here z is positive.
            if (Math.abs(za.hi()) > z.hi()) {
                final DD tmp = za;
                za = z;
                z = tmp;
            }
            // za < z
            if (Math.abs(zb) > Math.abs(z.hi())) {
                return za.add(z).add(zb).hi();
            }
            // za,zb < z
            return za.add(zb).add(z).hi();
        }

        // Pure double implementation
        double d = pow.apply(x, x > -0.5 ? 0 : 1, -s);
        if (!Double.isFinite(d)) {
            return odd && x > -0.5 ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        }
        // Here the all other terms cannot overflow.
        // Add the other dominant term.
        double pn;
        double pp;
        if (x > -0.5) {
            pn = d;
            pp = pow.apply(x, 1, -s);
            d += pp;
        } else {
            pn = Math.pow(x, -s);
            pp = d;
            d += pn;
        }

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2; any individual term x^-s < 1.
        // Exit early if remaining terms cannot be added.
        if (Math.abs(d) > 0x1p53) {
            // Limit of double arithmetic
            return d;
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) + (-1)^s * [ zeta(s, -x) - zeta(s, 1 - a) ]
        // z + za + zb

        // We have the two largest power terms for each side.
        // The remaining terms are increasing smaller.

        double z = zeta(s, x + 1, pp, c);
        double za;
        double zb;
        // A single call to zeta uses many pow operations;
        // use a direct sum when zeta will use more.
        if (x - a > 2 * (c.getN() + 2)) {
            // Note: The difference (za - zb) incurs cancellation.
            // This should not be an issue as x in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 1)  = 1.6449
            // zeta(2, 2)  = 0.6449
            // zeta(2, 30) = 0.0339
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            za = zeta(s, -x, Math.abs(pn), c);
            zb = zeta(s, 1 - a, pow.apply(-a, 1, -s), c);
            // Both terms are positive. Correct the sign for final addition.
            if (odd) {
                za = -za;
            } else {
                zb = -zb;
            }
        } else {
            // Sum terms in ascending order of magnitude
            // Using a double to track the iterations is fine as (a+n) is exact until > x.
            za = pn;
            zb = 0;
            // Instead of (a+k) use a+k = x+i with i in [a-x, 0) to use the extended precision power.
            for (int i = (int) (a - x); i < 0; i++) {
                zb += pow.apply(x, i, -s);
            }
        }

        // Sum in magnitude order. Here z is positive.
        if (Math.abs(za) > z) {
            final double tmp = za;
            za = z;
            z = tmp;
        }
        // za < z
        if (Math.abs(zb) > Math.abs(z)) {
            final double tmp = zb;
            zb = z;
            z = tmp;
        }
        // za,zb < z
        return DD.ofSum(za, zb).add(z).hi();
    }

    /**
     * Test the factors required for the tail sum. These are computed from the numerator
     * and denominator of the Bernoulli numbers, and the factorial of 2k. The test asserts
     * that using 2k! / B_2k is more accurate than B_2k / 2k! when limited to double precision.
     */
    @Test
    void testFactors() {
        // Factorial of 2k. Initialise at k=0.
        BigInteger factorial = BigInteger.ONE;
        double sum1 = 0;
        double sum2 = 0;
        // Enough digits for quad-double precision.
        // Sufficient for testing the BigDecimal implementation with different precision.
        final MathContext mc = new MathContext(70);
        // Check precomputed factors required for N ~ M ~ P for P-bits of precision.
        final int maxK = 106;
        Assertions.assertTrue(BernoulliNumbers.isSupported(2 * maxK));
        for (int k = 1; k <= maxK; k++) {
            // 2k!
            factorial = factorial.multiply(BigInteger.valueOf((2L * k - 1) * (2L * k)));
            final BigFraction b2k = BernoulliNumbers.value(2 * k);

            // 2k! / B_2k
            final BigFraction factor1 = b2k.reciprocal().multiply(factorial);
            final double d1 = factor1.doubleValue();
            // Cross verify BigFraction vs BigDecimal
            BigDecimal v = factor1.bigDecimalValue(mc);
            Assertions.assertEquals(v.doubleValue(), d1);
            // Find ULP precision
            final double e1 = new BigDecimal(d1).subtract(v)
                .divide(new BigDecimal(Math.ulp(d1)), mc).doubleValue();
            sum1 += Math.abs(e1);

            Assertions.assertEquals(d1, RF[k - 1]);

            // Format to print the table:
            // "%s, // %s! / (%s / %s)%n", d1, 2 * k, num, denom

            // B_2k / 2k!
            final BigFraction factor2 = BigFraction.of(b2k.getNumerator(), factorial.multiply(b2k.getDenominator()));
            final double d2 = factor2.doubleValue();
            // Cross verify BigFraction vs BigDecimal
            v = factor2.bigDecimalValue(mc);
            Assertions.assertEquals(v.doubleValue(), d2);
            // Find ULP precision
            final double e2 = new BigDecimal(d2).subtract(v)
                .divide(new BigDecimal(Math.ulp(d2)), mc).doubleValue();
            sum2 += Math.abs(e2);

            Assertions.assertEquals(d2, F[k - 1]);
            Assertions.assertEquals(v, FBD[k - 1]);
            Assertions.assertEquals(DD.from(v), FDD[k - 1]);

            // Format to print the table:
            // "%s, // (%s / %s) / %s!%n", d2, num, denom, 2 * k

            // For any M, the cumulative error is lower using 2k! / B_2k
            // as the first 6/7 factors are exact and errors in the later factors
            // are comparable.
            Assertions.assertTrue(sum1 < sum2, "2k! / B_2k does not have lower combined error");
        }
    }

    @Test
    void testFractionalAddAssumptions() {
        // Adding double-double numbers.
        // With a fixed mantissa of 106 bits, numbers a+b cannot be added if they do not
        // overlap by 0.5 ulp. This is 2^-106 when the first term is in [1, 2) and b is smaller.
        // However a double-double has a floating mantissa of *at least* 106 bits.
        // When using normalised numbers such as a fractional representation the low part
        // exponent changes if the round-off from the high part has leading zeros. This expands
        // the effective mantissa:
        // x0 1-------
        // x1         10000000
        // x1         01-------
        // x1         001-------
        // A low part x1 with no leading zeros must have all trailing zeros
        // otherwise it will be larger than 0.5 ulp(x0).
        // So the practical working mantissa of a DD is 107+ bits in all but a single case.
        // With random low parts we expect the probability we can add the numbers
        // to reduce by 2 for each halving of the second term.
        // log2(scale)  p(add)
        // -107         1.000000
        // -108         0.500000
        // -109         0.250000
        // -110         0.125000
        // -111         0.062500
        // -112         0.031250
        // -107+b       2^-b

        // Special case of a 106-bit mantissa with a non-zero low part.
        // The low part is equal to 0.5 ulp(high part) and a number smaller
        // by a scale factor of 2^106 cannot be added
        final DD a = DD.ofSum(1.0, Math.scalb(1.0, -53));
        Assertions.assertEquals(a, a.add(DD.ONE.scalb(-106)), "Cannot add 1 to 1e-106");

        // Now test random mantissas.
        // Fixed seed for stability
        final SplittableRandom rng = new SplittableRandom(-2471564807353423733L);
        // A sign for each number: 1 - 0 or 1 - 2
        final IntSupplier sign = () -> 1 - (rng.nextInt() & 0x2);

        final int maxB = 10;
        final int n = 1 << (maxB + 4);
        for (int b = 0; b < maxB; b++) {
            int added = 0;
            for (int i = 0; i < n; i++) {
                // Create a number with all random bits.
                // The random numbers must have a leading 1 bit.
                // Create numbers +/- [1, 2) +/- (0, 1] scaled to be smaller that 0.5 ulp(1.0) = 2^-53.
                final double x0 = sign.getAsInt() * (1 + rng.nextDouble());
                final double x1 = sign.getAsInt() * Math.scalb(1 - rng.nextDouble(), -53);
                final double y0 = sign.getAsInt() * (1 + rng.nextDouble());
                final double y1 = sign.getAsInt() * Math.scalb(1 - rng.nextDouble(), -53);
                final DD x = DD.ofSum(x0, x1);
                final DD y = DD.ofSum(y0, y1);
                Assertions.assertEquals(x0, x.hi());
                Assertions.assertEquals(x1, x.lo());
                Assertions.assertEquals(y0, y.hi());
                Assertions.assertEquals(y1, y.lo());
                // Two double-double numbers in [1, 2).
                final DD s = x.add(y.scalb(-107 - b));
                if (!x.equals(s)) {
                    added++;
                }
            }
            // We could do a statistical test across all observations.
            // Here we just check the observed additions is approximately as expected.
            // This will fail if the additions are out by 2 which would occur if
            // the scaling base of 107 was incorrect.
            // p(added) = 2^-b
            final double expected = Math.scalb(1, -b);
            final double actual = (double) added / n;
            Assertions.assertEquals(expected, actual, expected * 0.25);
        }
    }

    @Test
    void testAddZero() {
        // Zero has an exponent of zero in the DD fractional representation.
        final int[] exp = {-42};
        final DD zero = DD.ZERO.frexp(exp);
        final int ezero = exp[0];
        Assertions.assertEquals(0, ezero);
        // Add any number to zero in fractional representation and it
        // should always be returned appropriately scaled
        final DD x = DD.ofSum(3.45, 5.67).frexp(exp);
        final int ex = exp[0];
        DD y;
        for (final int b : new int[] {0, 1, -13, -263784, 26378468, Integer.MIN_VALUE, Integer.MAX_VALUE}) {
            // Testing implementation (in this class)
            y = HurwitzZetaTest.add(zero, ezero, x, ex + b, exp, DD::add);
            Assertions.assertEquals(x, y, () -> Integer.toString(b));
            Assertions.assertEquals(ex + b, exp[0]);
            y = HurwitzZetaTest.add(x, ex + b, zero, ezero, exp, DD::add);
            Assertions.assertEquals(x, y, () -> Integer.toString(b));
            Assertions.assertEquals(ex + b, exp[0]);

            // Final version
            y = HurwitzZeta.add(zero, ezero, x, ex + b, exp);
            Assertions.assertEquals(x, y, () -> Integer.toString(b));
            Assertions.assertEquals(ex + b, exp[0]);
            y = HurwitzZeta.add(x, ex + b, zero, ezero, exp);
            Assertions.assertEquals(x, y, () -> Integer.toString(b));
            Assertions.assertEquals(ex + b, exp[0]);
        }

        // Zero added to zero.
        // Do not test scaled zero. DD.frexp will not return a scale for zero
        // so it will never be passed to the add utility function.
        y = HurwitzZetaTest.add(zero, ezero, zero, ezero, exp, DD::add);
        Assertions.assertEquals(zero, y);
        Assertions.assertEquals(ezero, exp[0]);
        y = HurwitzZeta.add(zero, ezero, zero, ezero, exp);
        Assertions.assertEquals(zero, y);
        Assertions.assertEquals(ezero, exp[0]);
    }

    /**
     * Test the precision of the double implementation of the zeta function.
     *
     * <p>Note convergence requires:
     * <pre>
     * M   N    eps
     * 8   10   2^-53
     * 10  8
     * 12  7
     * 15  6
     * 53  6
     * 15  8    2^-63
     * 26  7
     * 53  7
     * </pre>
     *
     * <p>Smaller M can prevent tail convergence before the last term and error
     * can increase due to a truncated tail. If N is too low for the size of the
     * tail it can diverge and error is extreme. This can be handled by
     * increasing N; or a emergency brake in the tail to stop when the terms
     * diverge. The second option will avoid extreme error. It will not provide
     * an accurate result for some part of the zeta domain.
     *
     * <p>The performance of the function is very dependent on N as it is
     * proportional to the number of power functjon calls. The ideal parameters
     * (N, M) should have low error across the entire domain, and ensure
     * convergence of the tail within M.
     *
     * @see #testTailConvergence2(double, double, boolean, int, int, int, int)
     */
    @ParameterizedTest
    @CsvSource({
        // Notes:
        // - N >= 7 has similar error. A lower N uses more M in the tail series.
        // - Using a high precision pow(n + x, y) lowers max error and rms.
        //   The power function computes (a+n)^-s using -s * round-off of (a+n).
        //   so the difference is greater when s is large and a is small
        //   (verified using test resources with s above 4).
        //   It does not impact a in [0, 1] where (a+n) is inexact as a^-s is magnitudes
        //   larger than all other terms.
        //   It could be applied when: a>1; s>2; and (a+n) is inexact. Here the round-off
        //   multiplied by s will be above 2^-53 and exp(-s * round-off) != 1.
        // - Using the power function to compute (a+n)^-(2k-1+s) in the tail makes no
        //   difference and will be expensive.
        // - Extended precision sum makes a small difference. Max error is similar and RMS drops.
        //   Terms are added in ascending magnitude which reduces cumulative error in double
        //   precision unless the largest terms have error. Can be helped using high precision pow.
        // - Extended precision tail sum makes no difference. The major error is in the other part
        //   of the computation.
        // - Using a smaller epsilon requires more terms in the tail but does not impact
        //   error unless epsilon is too high and the tail series can diverge.
        // - RMS drops with increasing N then plateaus. The plateau is at larger N when
        //   using the higher precision power **and** sum, e.g. N=8 vs N=10. Using one or
        //   the other N=8 is OK.
        // - The tail power term can be be very small and underflow. This can be mitigated by
        //   scaling up the power term, and scaling down another factor by an equivalent amount
        //   to extend the number of terms before underflow, or convergence.

        // TODO:
        // Add option to record errors and M to a file and plot them.
        // Require domain data for a in [0, 1].
        // Require a sub-normal curve at 2^1070 for reference.
        //
        // Test the domain data for all a > 2^b with the asymptotic formula.
        // Compare to the full evaluation then the threshold is OK.

        // The full domain data requires tail scaling: +256.
        // Without this option the tail cannot be computed when
        // the power terms -> 0.

//        // Extended precision power function (big improvement if used for all terms)
//        "7, 15, 15, -53, 256, false",
//        "7, 15, 15, -53, 257, false",
//        "7, 15, 15, -53, 259, false",
//        // With extended precision sum: +4
//        // Difference is not visible due to other errors.
//        "7, 15, 15, -53, 260, false",
//        "7, 15, 15, -53, 262, false",
//        "7, 15, 15, -53, 263, false",
//
//        // Sub-normal term: +512. Large drop in max error.
//        "7, 12, 15, -53, 771, false",
//        "7, 12, 26, -53, 771, false",
//        "7, 12, 53, -53, 771, false",
//        // With extended precision sum: +4
//        // The extended precision sum makes a small difference only when
//        // the major errors have been removed.
//        // RMS drops as N increases. Max error is variable.
//        "7, 12, 15, -53, 775, false",
//        "7, 12, 26, -53, 775, false",
//        "7, 12, 53, -53, 775, false",
//
//        // M convergence
//        // Negligible error change unless too high (> -53) where the
//        // tail may have remaining terms that can be added to the sum.
//        // If set much lower than 53-bits then N and M must be increased
//        // or tail divergence is observed.
//        "7, 10, 15, -50, 775, false", // Worse due to truncated tail
//        "8, 10, 15, -63, 775, false",
//        "7, 10, 26, -63, 775, false",
//
//        // Sub-normal zeta result and/or terms.
//        // Requires scaling of the tail series terms (option 256).
//        "7, 12, 15, -53, 263, true", // 7 + 256
//        // Using summation of sub-normal terms to prevent underflow.
//        "7, 12, 15, -53, 775, true", // 7 + 256 + 512
//        // Using summation of sub-normal terms to prevent underflow.
//        "7, 12, 18, -63, 775, true", // 7 + 256 + 512

        // Match the final implementation here
        "7, 12, 15, -53, 775, false",
    })
    void testPrecisionDouble(int ln, int un, int maxM, int bt, int options, boolean subNormal)
        throws IOException {
        // Default to not use variable convergence on series sum S.
        // Early convergence can be very bad and requires adjusting the tolerance.
        // The zeta domain data is very challenging.
        // Lower tolerance when using extended precision.
        final int extFlags = Context.EXT_SUM | Context.EXT_POW_TERMS | Context.EXT_POW_SERIES;
        final boolean extPrecision = (options & extFlags) == extFlags;
        final double maxTolerance = extPrecision ? 10 : 100;
        final double rmsTolerance = extPrecision ? 1 : 10;
        final double eps = Math.scalb(1.0, bt);
        for (int n = ln; n <= un; n++) {
            // Reset
            Arrays.fill(N, 0);
            Arrays.fill(M, 0);
            String resultFilename = null;
            if (doReporting()) {
                resultFilename = String.format("zeta_n%d_m%d_b%d_options%d.csv",
                    n, maxM, -bt, options);
            }

            final Context c = Context.of(n, maxM)
                .withEpsilon(eps)
                .withOptions(options);
            final DoublePrecisionTestCase test = new DoublePrecisionTestCase() {
                @Override
                public double getTolerance() {
                    return maxTolerance;
                }

                @Override
                public double getRmsTolerance() {
                    return rmsTolerance;
                }

                @Override
                public DoubleBinaryOperator getFunction() {
                    return (s, a) -> HurwitzZetaTest.zeta(s, a, c);
                }

                @Override
                public String[] getFilenames() {
                    // Use combined data from multiple resources in order to
                    // find N and M values. Sub-normal data is for testing tail scaling
                    return subNormal ? SUB_NORMAL_RESOURCES : TEST_RESOURCES;
                }

                @Override
                public String toString() {
                    return String.format("ZETA %s  2^%d  %3d",
                        zetaName(c.getN()), bt, options);
                }
            };
            assertFunction(test, true, resultFilename);
            // When not reporting cut the test short
            if (!doReporting()) {
                break;
            }
        }
    }

    /**
     * Test the precision of the BigDecimal implementation of the zeta function.
     *
     * <p>Note convergence requires:
     * <pre>
     * M   N
     * 106 24
     * 53  25
     * </pre>
     * <p>Smaller M can prevent tail convergence before the last term and error
     * can increase.
     *
     * @see #testTailConvergenceDD(int, boolean, int, int, int, int)
     */
    @ParameterizedTest
    @CsvSource({
        // Full double-double precision (~34 digits)
        // Requires N=13. Larger N uses smaller M.
        "13, 30, 34, 0, 105",
        // Use pow in the tail series (no difference)
        "13, 30, 34, 16, 105",
        // If the computation is limited to fewer bits it cannot achieve double-double precision.
        // This has implications for the DD version because DD arithmetic is typically
        // performed to a few eps of 2^-106.
        "13, 30, 33, 0, 105",
        // Not enough
        "13, 30, 32, 0, 105",
        // Full double precision (~17 digits)
        // Requires N=7. Larger N uses smaller M.
        "7, 15, 18, 0, 52",
        "7, 15, 17, 0, 52",
        // Not enough
        "7, 15, 16, 0, 52",
        // Quad double precision. Requires test resources to have more than 68 digits of precision
        "30, 40, 70, 0, 211",
        // Not enough
        "35, 40, 65, 0, 211",
    })
    // TODO - update with a maxM
    void testPrecisionBigDecimal(int ln, int un, int digits, int options, int scale)
        throws IOException {
        for (int n = ln; n <= un; n++) {
            // Reset M
            Arrays.fill(M, 0);

            final Context c = Context.of(n, FBD.length)
                .withMathContext(new MathContext(digits))
                .withOptions(options);
            final ExtendedPrecisionTestCase test = new ExtendedPrecisionTestCase() {
                @Override
                public double getTolerance() {
                    // Do not fail individual cases
                    return 1e5;
                }

                @Override
                public double getRmsTolerance() {
                    // Within 10 bits of the target precision
                    return Math.scalb(1.0, 10);
                }

                @Override
                public int scale() {
                    return scale;
                }

                @Override
                public BiFunction<Integer, Double, BigDecimal> getFunction() {
                    return (s, a) -> HurwitzZetaTest.zeta(s, new BigDecimal(a), null, c);
                }

                @Override
                public String[] getFilenames() {
                    // Use resources to push zeta(s, a) beyond the limit where it will
                    // realistically be used
                    return EXT_INT_TEST_RESOURCES;
                }

                @Override
                public String toString() {
                    return String.format("ZETA %s  dps=%3d  %2d", zetaName(c.getN()), digits, options);
                }
            };
            assertFunction(test, 0, true);
            // When not reporting cut the test short
            if (!doReporting()) {
                break;
            }
        }
    }

    /**
     * Test the precision of the DD implementation of the zeta function.
     *
     * <p>Note convergence requires:
     * <pre>
     * M   N
     * 106 12
     * 26  13
     * </pre>
     * <p>Smaller M can prevent tail convergence before the last term and error
     * can increase.
     *
     * <p>Differences between some options may only be observed when using a sub-set of
     * the {@link #EXT_INT_TEST_RESOURCES}. The current resources include data for zeta(s, a)
     * below the sub-normal range with a ~ s and much larger than 1. The current use case
     * for the negative a evaluation will use {@code 0 < a < 1} for the major terms on
     * both sides. Extended precision for positive a may be required for a polygamma
     * implementation but is not required for the zeta implementation.
     *
     * <p>These data contain cases where the result is not accurate to 104+ bits using DD
     * pow even when the first term is computed using DDMath pow. The next term in the
     * series (a+1)^-s can have multiple bits of error when s is large and this can
     * propagate up to the exact first term a^-s when a is far above 1. This can be fixed
     * by using DDMath for all terms, or optionally choosing to use DDMath when terms
     * (a+k)^-s are expected to have bit errors that overlap with a^-s. The error bits can
     * be approximated using ~ 16 * s, although this is a conservative upper bound and it
     * may be possible to use a few bits lower for the bound. E.g. (a+1)^-100 is expected to
     * have 1600 ulp error, or log2(1600) ~ 10.6 bits. If this overlaps with a^-100 by
     * (106 - 11) = 95-bits or more then compute using DDMath. Continue with DDMath
     * until (a+k)^-s has no overlap with a^-s. When using the ascending sum starting at
     * (a+N)^-s the running sum must be compared with the known term a^-s. When it approaches
     * overlap then switch to DDMath pow for remaining terms.
     *
     * @see #testTailConvergenceDD(int, boolean, int, int, int, int)
     */
    @ParameterizedTest
    @CsvSource({
        // Extended precision (ulp ~ 2^-105)
        // DD.pow
        "13, 30, -106, 0, 105, 0",
        // DD.pow and DDMath pow are similar accuracy when s = 2 and a < 1.
        // When s is larger and a > 1 the DDMath pow gains a few bits in the result
        // but the max is ~105 bits
        "13, 30, -106, 3, 105, 0",
        // DDMath.pow with DD.pow in the sum of the series.
        // This is worse than DDMath when a is above 1, i.e. DD.pow cannot
        // be selectively used for the *same* precision. However zeta evaluations
        // with a above 1 are only used: (1) as a term added to a much larger
        // zeta evaluation and precision does not require all the bits; (2)
        // in the total cancellation path which requires double precision
        "13, 30, -106, 1, 105, 0",
        // Ignore below MIN_VALUE = 2^-1074
        // Has precision to ~100 bits.
        "15, 20, -106, 1, 105, -1074",

        // Maximum precision. Accurate add and DDMath.pow.
        // Max ulp is ~1.
        // Note: Cancellation in the negative zeta implementation has not been observed
        // above 55-bits. This level of precision may only gain ~ 1 ulp around the roots.
        "15, 20, -106, 195, 105, 0", // 3 vs 3 + 64 + 128

        // Use accurate add for all addition. Does not change accuracy when
        // the power term may have several bits of error.
        "15, 20, -106, 193, 105, 0", // 1 vs 1 + 64 + 128

        // Use pow in the tail series (no difference)
        "13, 30, -106, 19, 105, 0", // 3 vs 3 + 16

        // Higher tail convergence.
        // No difference to max error. RMS is lower. Higher M.
        // "13, 20, -108, 3, 105, 0",
        // "13, 20, -110, 3, 105, 0",

        // Double precision (ulp ~ 2^-52)
        // Ignore below MIN_VALUE = 2^-1074
        "7, 15, -53, 0, 52, -1074",
        // Extra tail terms lower error in double precision
        "7, 15, -58, 0, 52, -1074",
        // Not enough
        "7, 15, -48, 0, 52, -1074",
    })
    // TODO - update with a maxM
    void testPrecisonDD(int ln, int un, int b, int options, int scale, int ignoreB)
        throws IOException {
        final double eps = Math.scalb(1.0, b);
        for (int n = ln; n <= un; n++) {
            // Reset M
            Arrays.fill(M, 0);

            final Context c = Context.of(n, FDD.length)
                .withEpsilon(eps)
                .withOptions(options);
            final ExtendedPrecisionTestCase test = new ExtendedPrecisionTestCase() {
                @Override
                public double getTolerance() {
                    // Do not fail individual cases
                    return 1e5;
                }

                @Override
                public double getRmsTolerance() {
                    // Within 10 bits of the target precision
                    return Math.scalb(1.0, 10);
                }

                @Override
                public int scale() {
                    return scale;
                }

                @Override
                public BiFunction<Integer, Double, BigDecimal> getFunction() {
                    return (s, a) -> HurwitzZetaTest.zeta(s, DD.of(a), c);
                }

                @Override
                public String[] getFilenames() {
                    // Use resources to push zeta(s, a) beyond the limit where it will
                    // realistically be used
                    return EXT_INT_TEST_RESOURCES;
                }

                @Override
                public String toString() {
                    return String.format("ZETA %s  2^%-4d  %3d", zetaName(c.getN()), b, options);
                }
            };
            assertFunction(test, ignoreB, true);
            // When not reporting cut the test short
            if (!doReporting()) {
                break;
            }
        }
    }

    /**
     * Create a formated name for the test zeta method using the N and M used during evaluation.
     *
     * @param fixedN Fixed N.
     * @return the name
     */
    private static String zetaName(int fixedN) {
        // Get mean/max
        int maxN = 0;
        double meanN = 0;
        int n = 0;
        int maxM = 0;
        double meanM = 0;
        int m = 0;
        if (fixedN == 0) {
            for (int i = 0; i < N.length; i++) {
                if (N[i] != 0) {
                    n += N[i];
                    maxN = i;
                    meanN = i * N[i];
                }
            }
            if (n != 0) {
                maxN++;
                meanN = 1 + meanN / n;
            }
        }
        for (int i = 0; i < M.length; i++) {
            if (M[i] != 0) {
                m += M[i];
                maxM = i;
                meanM = i * M[i];
            }
        }
        // Summarize. Note the values are indices from 0 so add 1.
        if (m != 0) {
            maxM++;
            meanM = 1 + meanM / m;
        }
        if (fixedN != 0) {
            return String.format("%2d %2d [%5.2f]", fixedN, maxM, meanM);
        }
        return String.format("%2d %2d [%5.2f %5.2f]", maxN, maxM, meanN, meanM);
    }

    /**
     * Convert a fractional representation of a DD to BigDecimal.
     *
     * @param f the fractional representation in [0.5, 1)
     * @param exp the exponent
     * @return the BigDecimal
     */
    private static BigDecimal toBigDecimal(DD f, int exp) {
        // Convert to BigDecimal without losing precision for sub-normals.
        // r in [0.5, 1.0). Since the low part of the DD is floating
        // the effective mantissa is 106+ bits.
        // Try scaling the low part. If a normal number then scaling is allowed.
        // Ignore scaling of a zero low part.
        if (Math.scalb(Math.abs(f.lo()), exp) >= Double.MIN_NORMAL ||
            f.lo() == 0 && Math.scalb(Math.abs(f.hi()), exp) >= Double.MIN_NORMAL) {
            return f.scalb(exp).bigDecimalValue().round(MathContext.DECIMAL128);
        }
        // scale by 2^exp = 1 << exp
        final BigDecimal v = f.bigDecimalValue();
        if (exp < 0) {
            return v.divide(new BigDecimal(BigInteger.ONE.shiftLeft(-exp)), MathContext.DECIMAL128);
        }
        return v.multiply(new BigDecimal(BigInteger.ONE.shiftLeft(exp)), MathContext.DECIMAL128);
    }

    @ParameterizedTest
    @CsvSource({
        // s = 1 : infinity
        "1.0, 1.0, Infinity",
        "1.0, 345.6, Infinity",
        "1.0, -345.6, Infinity",
        // s < 1 : not convergent so unsupported (nan)
        "0.23, 1.0, NaN",
        "-1.23, 1.0, NaN",
        "-Infinity, 1.0, NaN",
        // s = infinity; a < 1 : infinity
        "Infinity, 0.5, Infinity",
        // s = infinity; a = 1 : 1 (matches Riemann zeta)
        "Infinity, 1.0, 1.0",
        // s = infinity; a > 1 : 0
        "Infinity, 1.5, 0.0",
        // a = negative integer or 0 : infinity
        "1.23, 0.0, Infinity",
        "1.23, -1.0, Infinity",
        "1.23, -1e300, Infinity",
        // a <= 0 not integer; s != integer : complex result (nan)
        "1.23, -0.5, NaN",
        "1.23, -1.5, NaN",
        // a = -infinity : nan
        "2.0, -Infinity, NaN",
        // nan arguments : nan
        "2.0, NaN, NaN",
        "NaN, 1.0, NaN",
        "NaN, NaN, NaN",
    })
    void testZetaSpecial(double s, double a, double z) {
        Assertions.assertEquals(z, HurwitzZeta.value(s, a));
    }

    /**
     * Spot tests for the zeta function to check various points in the domain and extreme values.
     */
    @ParameterizedTest
    @MethodSource(value = "testZetaSpot")
    void testZetaSpot(double s, double a, double z, int ulp) {
        assertClose(HurwitzZeta::value, s, a, z, ulp);
    }

    /**
     * Spot tests for the zeta function using {@code a = 1}. Uses data from the Riemman zeta test.
     */
    @ParameterizedTest
    @MethodSource(value = "org.apache.commons.numbers.gamma.RiemannZetaTest#testZetaSpot")
    void testZeta1Spot(double s, double z, int ulp) {
        // Cannot be called with s < 1
        Assumptions.assumeTrue(s >= 1);
        assertClose(HurwitzZeta::value, s, 1.0, z, ulp);
    }

    /**
     * Spot tests for the zeta function using {@code s = 2}. Uses data from the trigamma test
     * for a polygamma (psi) implementation:
     * <pre>
     * psi(n, x) = (-1)^(n+1) n! zeta(n+1, x).
     * psi(1, x) = -1^2 * 1! * zeta(2, x)
     *           = zeta(2, x)
     * </pre>
     */
    @ParameterizedTest
    @MethodSource(value = "org.apache.commons.numbers.gamma.TrigammaTest#testTrigamma")
    void testZetaTrigamma(double x, double expected) {
        // Note: The Trigamma function passes with relative/absolute error 1e-9.
        // The zeta function is within 1 ULP.
        TestUtils.assertEquals(expected, HurwitzZeta.value(2, x), 1, null, () -> Double.toString(x));
    }

    /**
     * Spot tests for the double-precision test implementation.
     */
    @ParameterizedTest
    @MethodSource(value = "testZetaSpot")
    void testZetaDoubleSpot(double s, double a, double z, int ulp) {
        if (a > 0) {
            // The final version calls RiemmanZeta and has different ULP tolerances
            final int u = a == 1 ? ulp + 1 : ulp;
            assertClose((x, y) -> HurwitzZetaTest.zeta(x, y, Context.DOUBLE), s, a, z, u);
        }
    }

    /**
     * Spot tests for the negative zeta test implementations.
     */
    @ParameterizedTest
    @MethodSource(value = "testZetaSpot")
    void testZetaNegativeSpot(double s, double a, double z, int ulp) {
        if (a < 0 && s < Integer.MAX_VALUE) {
            assertClose((x, y) -> HurwitzZetaTest.zetaNegativeBD((int) x, y), s, a, z, ulp);
            assertClose((x, y) -> HurwitzZetaTest.zetaNegativeDD((int) x, y), s, a, z, ulp);
        }
    }

    static Stream<Arguments> testZetaSpot() {
        final double inf = Double.POSITIVE_INFINITY;
        return Stream.of(
            // Reference values from mpmath version 1.4.1.
            // from mpmath import mp, zeta
            // mp.dps = 30; mp.pretty = True
            // def f(s, a):
            //   print(f'Arguments.of({s}, {a}, {zeta(s, a)}, 0),')
            // f(1.001, 1) etc.
            Arguments.of(1.001, 0.999, 1000.57893532226893974131793946, 1),
            Arguments.of(1.001, 1, 1000.57728847601162684806668989, 0),
            Arguments.of(1.001, 1.001, 1000.57564403708096426589811664, 0),
            Arguments.of(1.001, 3, 999.077634929516400548962369402, 0),
            Arguments.of(1.001, 156.78, 994.961087152463264230812967825, 1),
            Arguments.of(1.001, 345600, 987.327939500118316838908716646, 1),
            Arguments.of(1.001, 100000000.0, 981.747943025003254769065166246, 1),
            Arguments.of(1.001, 1000000000.0, 979.48998540929872849188230117, 1),
            Arguments.of(1.001, 10000000000.0, 977.237220955969649930501114761, 0),
            Arguments.of(1.001, 8.374e+19, 955.1620677642774549247109353, 0),
            Arguments.of(1.001, 7.2834e+238, 576.949320020341985622075858479, 0),
            Arguments.of(1.1678, 1, 6.54877176355186372373480299218, 1),
            Arguments.of(1.1678, 2.5, 5.29477174574315816738362779184, 1),
            Arguments.of(1.1678, 5.765, 4.50848224904053531500311915417, 1),
            Arguments.of(1.1678, 87698, 0.882622082916184877125466190512, 1),
            Arguments.of(1.1678, 1098765, 0.577489707637839655591145476126, 1),
            Arguments.of(1.1678, 100000000.0, 0.27089940045667952429754209707, 0),
            Arguments.of(1.1678, 1000000000.0, 0.184080609461973448176589687532, 1),
            Arguments.of(1.1678, 10000000000.0, 0.125085809513774573322512635067, 1),
            Arguments.of(1.1678, 6.786e+70, 7.75645430994177324455600419939e-12, 0),
            Arguments.of(1.3567, 1, 3.40603490577277492861753138946, 1),
            Arguments.of(1.3567, 12, 1.17295098623816176759602458034, 1),
            Arguments.of(1.3567, 267, 0.382340735582659556038655214163, 0),
            Arguments.of(1.3567, 100000000.0, 0.00392732544562110257556960028487, 0),
            Arguments.of(1.3567, 1000000000.0, 0.0017274158124759523689829065186, 0),
            Arguments.of(1.3567, 6780000000.0, 0.000872764873809386121835245990354, 1),
            Arguments.of(1.3567, 10000000000.0, 0.000759795803739606906900648768122, 1),
            Arguments.of(1.3567, 2.394279e+25, 2.48282011395930739748191376465e-09, 0),
            Arguments.of(1.3567, 1.37e+201, 5.03765781298315620677652110278e-72, 0),
            Arguments.of(1.9183, 1.234, 1.31039727875809754102393633608, 0),
            Arguments.of(1.9183, 2.234, 0.642314727500847848717986207886, 1),
            Arguments.of(1.9183, 32, 0.0458231265360237920443268300232, 1),
            Arguments.of(1.9183, 189, 0.00886331331419119681212176796288, 0),
            Arguments.of(1.9183, 26378, 9.48410976427270589834458588242e-05, 0),
            Arguments.of(1.9183, 1484793, 2.34193449523082289840265727697e-06, 1),
            Arguments.of(1.9183, 100000000.0, 4.90473353191884545392736492566e-08, 1),
            Arguments.of(1.9183, 1000000000.0, 5.91991424826323516752379877396e-09, 1),
            Arguments.of(1.9183, 10000000000.0, 7.14521688264220831940065076935e-10, 1),
            Arguments.of(1.9183, 7.18923e+79, 5.06551945929778382002806813887e-74, 1),
            Arguments.of(1.9183, 2.3423e+159, 4.87385597076733868223619309795e-147, 0),
            Arguments.of(1.9183, 3.4535e+303, 1.98532564832533048061338098412e-279, 1),

            // Large s
            Arguments.of(1001, 1, 1.0, 0),
            Arguments.of(1001, 1.3, 8.76403979411431591130552997055e-115, 0),
            Arguments.of(1001, 2, 4.66631809251609439495044772362e-302, 0),
            Arguments.of(1001, 3, 0, 0),

            Arguments.of(26783400000000.0, 1.0000000000000002, 0.994070539579705196633883729806, 0),
            Arguments.of(26783400000000.0, 1.0000000000002, 0.00470868956401873280652247209542, 0),
            Arguments.of(26783400000000.0, 1.0000000002, 0, 0),

            Arguments.of(1e+18, 1, 1.0, 0),
            Arguments.of(1e+18, 1.0000000000000002, 3.69192903582901397705695784205e-97, 0),
            Arguments.of(1e+18, 1.0000000000000004, 1.36303400055980247979692825583e-193, 1),

            Arguments.of(1e+19, 1, 1.0, 0),
            Arguments.of(1e+19, 1.0000000000000002, 0, 0),

            // Large s with large a
            Arguments.of(50, 20, 9.7412887436036448537718048483189669464194616553866e-66, 1),
            // XXX: mpmath disagrees with matlab
            Arguments.of(50, 200, 4.0877905547003724197061862269315573442276635335729e-115, 26000),
            // On this case using mp.dps=50,80,100 has differences in the first 12 digits
            Arguments.of(50, 2000, 3.6698119961414421907531687439376998787130928517259e-164, 600000),
            Arguments.of(50, 20000, 3.6296607820625266676514381487393861416851065066219e-213, 265),
            // --- matlab agrees with our implementation
            Arguments.of(50, 200, 4.0877905546774699327612047878897e-115, 0),
            Arguments.of(50, 2000, 3.6698119957034991027055454908981e-164, 0),
            Arguments.of(50, 20000, 3.6296607820623517558146318095834e-213, 0),
            // --- end matlab
            Arguments.of(50, 200000, 3.6256621473059150264825344395019993470971219502611e-262, 1),
            Arguments.of(50, 2000000, 3.625262448698370064130045968154440065060042043605e-311, 2),
            Arguments.of(50, 3000000, 8.5283690506160333026898323949253691719144047850959e-320, 3),
            Arguments.of(50, 3500000, 4.4717013098147129077153665995697691837355569545233e-323, 0),
            Arguments.of(50, 4000000, 0, 0), // 6.43e-326
            Arguments.of(5, 10000000000.0, 2.5000000005000000000416666666666666666666666666529e-41, 1),
            Arguments.of(5, 100000000000.0, 2.5000000000500000000004166666666666666666666666667e-45, 0),
            Arguments.of(5, 1000000000000.0, 2.5000000000050000000000041666666666666666666666667e-49, 1),
            Arguments.of(5, 10000000000000.0, 2.5000000000005000000000000416666666666666666666667e-53, 1),
            Arguments.of(5, 100000000000000.0, 2.5000000000000500000000000004166666666666666666667e-57, 1),
            Arguments.of(5, 1000000000000000.0, 2.5000000000000050000000000000041666666666666666667e-61, 0),

            // a^-s underflows
            // a = 2^-53 - 100 : Below asymptotic threshold
            Arguments.of(21, 9007199254740892.0, 4.04738577073239488356976513074e-321, 0),
            Arguments.of(21.1, 9007199254740892.0, 1.02223221121124030233647572377e-322, 0),
            Arguments.of(21.15, 9007199254740892.0, 1.62457991690996194160345547236e-323, 0),
            // Zeta ~ Double.MIN_VALUE
            Arguments.of(108, 1000, 9.85479262642312000966259642922e-324, 0),
            Arguments.of(107, 1060, 2.0597259802867861753987491314e-323, 0),

            // Requires M=12 when N=9.
            // This is largest M noted during development when the RMS error is close to optimal.
            Arguments.of(31.76, 8.23, 8.6811830191090714059294777379e-30, 1),
            Arguments.of(31.76, 11.23, 4.68180272588528321856530594913e-34, 0),
            Arguments.of(31.76, 12.23, 3.17258072050987981477830566397e-35, 0),
            Arguments.of(31.76, 13.23, 2.66312289027731475813712555339e-36, 1),
            Arguments.of(31.76, 14.23, 2.68377485329827353841985109724e-37, 2),
            Arguments.of(31.76, 15.23, 3.1669792006042583105605159352e-38, 1),
            Arguments.of(31.76, 17.23, 6.55526711653597399706064463689e-40, 0),
            Arguments.of(31.76, 19.23, 2.08909657211945256234757218604e-41, 1),

            Arguments.of(61.76, 30.23, 4.27875492887441652081612955782e-92, 2),

            // s -> large, a == 1
            // Asymptote of Riemann zeta function: 1 + 2^-s
            Arguments.of(25.67, 1, 1.00000001873152459253171691257, 0),
            Arguments.of(29.67, 1, 1.00000000117069191296622655835, 0),
            Arguments.of(39.67, 1, 1.00000000000114324712238181431, 0),
            Arguments.of(59.67, 1, 1.00000000000000000109028530523, 0),
            Arguments.of(69.67, 1, 1.00000000000000000000106473174, 0),
            Arguments.of(89.67, 1, 1.00000000000000000000000000102, 0),

            // a in [0, 1]
            Arguments.of(1.5, 0.5, 4.77653794755483324857662766936, 1),
            Arguments.of(1.5, 0.1, 34.0529755150756003469433380579, 1),
            Arguments.of(1.5, 0.9, 2.83731486390441065293824846471, 0),
            Arguments.of(1.234, 0.9, 5.06661932437970945786497350413, 0),
            Arguments.of(1.234, 0.567, 6.18461364443178731103712531596, 1),
            Arguments.of(1.234, 0.1567, 14.4639283884787232559726455653, 0),

            // a < 0 (non-integer) and s is a positive integer
            Arguments.of(2, -0.1567, 42.8461979972498360068058012169, 1),
            Arguments.of(3, -0.1567, -257.977656635792432293119990256, 0),
            Arguments.of(4, -0.1567, 1660.62037088089365709488623999, 0),
            Arguments.of(2, -5.1567, 44.004434644084193441570355425, 1),
            Arguments.of(3, -5.1567, -258.776505481263060632140556707, 0),
            Arguments.of(4, -5.1567, 1661.24004757051316284917148769, 1),

            Arguments.of(2, -1.00000000001567, 4072555449182211754851.99822909, 0),
            Arguments.of(3, -1.00000000001567, -2.59896546788095560752643092304e+32, 0),
            Arguments.of(2, -1.0000000000000002, 2.0282409603651670423947251286e+31, 0),

            // a close to half integer
            Arguments.of(2, -11.499, 9.7864096532064037170689523814012513283574378908625, 0),
            Arguments.of(2, -21.499, 9.8242530207688014634659878248899978896965690833573, 1),

            // -0.5 < a < 0 ; 1 + a may be inexact ; s is large : a -> -0.5
            Arguments.of(26, -0.49999999999999994, 134217728.00002640146373806618122596018906243147027, 0),
            // This has very large cancellation of terms either side of 0.
            // Computed with Math.pow(x, -s) + Math.pow(1+x, -s) the relative error is 0.023.
            Arguments.of(27, -0.49999999999999994, 1.6796301108582691430718404414343395058370129217631e-05, 5),

            // a is odd has cancellation.
            // ULP tolerance is very dependent on the JDK pow implementation.
            Arguments.of(3, -21.499, -0.09637775418460447271221850760411324846873791133448, 1),
            Arguments.of(3, -21.501, 0.098442803974649377073746600966748291488497220978961, 1),
            Arguments.of(5, -21.499, -0.64094298652061283507484770384318604647832428321319, 1),
            Arguments.of(5, -21.501, 0.64094511727200934057870930923285174211166495014722, 1),
            Arguments.of(3, -7.5000000001, 0.007782265638668063994730919817396152634504258389332, 1),
            Arguments.of(3, -7.4999999999, 0.0077822461572358593294768713462220939956242365981849, 1),
            Arguments.of(3, -7.499999999999999, 0.0077822558978654467309133842121074379966549149936996, 2),
            Arguments.of(3, -7.500000000000001, 0.0077822558980384765932871763254914638339941060968054, 3),

            // s is odd and a is half-integer -> total cancellation
            Arguments.of(3, -12.5, 0.0029542182928941203954486528445780501312428003881805, 0),
            Arguments.of(5, -23.5, 7.5243259573223049042413305168073933156343070954263e-07, 1),
            // s is even and a is half-integer -> no cancellation (sum below 0 is mirrored above 0)
            Arguments.of(2, -12.5, 9.792719176487780248390899866936566735096179245964, 0),
            Arguments.of(4, -23.5, 32.469672919579233053403733014746340086149182248433, 0),

            // Overflow the sum of the series with a < 0
            Arguments.of(18, -1.0000000000000002, 5.8086597987413400890549316334e+281, 0),
            Arguments.of(19, -1.0000000000000002, -2.61598781051334795153424084243e+297, 0),
            Arguments.of(20, -1.0000000000000002, inf, 0), // 1.18e313
            Arguments.of(21, -1.0000000000000002, -inf, 0), // -5.31e328

            // Overflow the sum of the series with a > 0
            Arguments.of(18, 2e-16, 3.81469726562500143524108492804e+282, 0),
            Arguments.of(19, 2e-16, 1.90734863281250075748835037869e+298, 0),
            Arguments.of(20, 2e-16, inf, 0), // 9.54e313

            Arguments.of(18, -0.9999999999999998, 5.8086597987413400890549316334e+281, 0),
            Arguments.of(19, -0.9999999999999998, 2.61598781051334795153424084243e+297, 0),
            Arguments.of(20, -0.9999999999999998, inf, 0), // 1.178e313
            Arguments.of(21, -0.9999999999999998, inf, 0), // 5.31e328
            // s > 1024 0.5^-s overflows with odd or even s
            Arguments.of(1025, -1.25, -inf, 0), // -1.29e617
            Arguments.of(1024, -1.75, inf, 0), // 3.23e616
            // s > 1065 (0.5+/-u)^-s with u=2^-54 overflows with odd or even s
            Arguments.of(1095, -1.25, -inf, 0), // -4.61e+661
            Arguments.of(1098, -1.75, inf, 0), // 1.15e+661

            // Eventually [0.5 +/- ulp(a)]^-s - [0.5 -/+ ulp(a)]^-s will overflow.
            // Even s is limited by 0.5^-1024, or 2^-1023
            Arguments.of(1020, -0.49999999999999994, 2.24711641857789488466163150291e+307, 0),
            Arguments.of(1022, -0.49999999999999994, 8.98846567431157953864652601186e+307, 0),
            Arguments.of(1024, -0.49999999999999994, inf, 0), // 3.56e+308
            Arguments.of(1020, -0.5000000000000001, 2.24711641857789488466163154618e+307, 0),
            Arguments.of(1022, -0.5000000000000001, 8.98846567431157953864652618561e+307, 0),
            Arguments.of(1024, -0.5000000000000001, inf, 0), // 3.56e+308
            // For code coverage
            Arguments.of(1068, -0.49999999999999994, inf, 0), // 6.32e+321
            Arguments.of(1068, -0.5000000000000001, inf, 0), // 6.32e+321
            // Odd s has cancellation and s > 1024
            Arguments.of(1061, -0.49999999999999994, -5.82077763699010369293454836251e+306, 0),
            Arguments.of(1063, -0.49999999999999994, -2.332699954051076428120424094e+307, 0),
            Arguments.of(1065, -0.49999999999999994, -9.34835541322444551626811537202e+307, 0),
            Arguments.of(1067, -0.49999999999999994, -inf, 0), // -3.75e+308
            Arguments.of(1061, -0.5000000000000001, 1.1641555273980207385869096725e+307, 0),
            Arguments.of(1063, -0.5000000000000001, 4.66539990810215285624084818801e+307, 0),
            Arguments.of(1065, -0.5000000000000001, inf, 0), // 1.87e+308
            Arguments.of(1067, -0.5000000000000001, inf, 0), // 7.49e+308
            // very large s: will be odd if cast to a long which will incorrectly trigger the
            // cancellation path with s half-integer; s should be detected as even.
            Arguments.of(1e+19, -1.25, inf, 0), // 1.88e+6020599913279623904
            Arguments.of(1e+19, -1.5, inf, 0), // 2.74e+3010299956639811952
            Arguments.of(1e+19, -1.75, inf, 0), // 1.88e+6020599913279623904

            // When a is < -1 and the number of terms is small a summation is used.
            // Hit the edge case where the first term has an exponent much smaller than 0.
            // Adding the fractional representation to zero should be the first term
            Arguments.of(182, -1.5000000000000002, 1.22599643269271108668667762575e+55, 0),

            // Large a and/or s have (a+n)^-s as sub-normal
            Arguments.of(143, 124.5, 3.61590594362399955089661789179e-300, 0),

            // Large a and/or s have a^-s as sub-normal
            Arguments.of(147, 124.5, 1.48298686245279271827660650487e-308, 2),
            Arguments.of(148, 124.5, 1.18689831491506998542539609447e-310, 4),
            Arguments.of(149, 124.5, 9.499650589522059112318304392e-313, 1),
            Arguments.of(152, 124.5, 4.87185978234820120308892496184e-319, 1),

            // Note: large negative a can be evaluated as complex using mpmath.
            // To obtain a real result requires the digits of precision (dps)
            // to be set above |a|. Larger a are tested using Matlab.
            // mp.dps = 70
            Arguments.of(2, -66.26738, 17.78437226907544744370996767080669047664124490204306137519070727420819, 0),
            Arguments.of(3, -66.26738, -50.12250472489562851240391661869509127911964240532672215389564430564016, 0),

            // s is odd and a is half-integer -> total cancellation
            // Test with large s: zeta(s, a) = zeta(s, 1 - a)
            // Note: mpmath or MATLAB do not recognise this case and compute
            // with the configured decimal precision. This can suffer cancellation
            // so cases are checked using f(s, 1 - a).
            // Cross verified using https://www.wolframalpha.com/input?i=HurwitzZeta
            Arguments.of(1199, -0.5, 7.35496130541141006844938414845e-212, 0),
            Arguments.of(1599, -0.5, 2.69201062691590466247243270165e-282, 0),
            Arguments.of(1799, -0.5, 1.62863948177150213207306826673e-317, 0),
            Arguments.of(1821, -0.5, 2.17679363860391259822822464064e-321, 0),
            Arguments.of(1831, -0.5, 3.77489320044438771289217773716e-323, 0),
            Arguments.of(1835, -0.5, 7.45657916137163005015738812279e-324, 0),
            Arguments.of(1837, -0.5, 3.31403518283183557784772805457e-324, 0),
            Arguments.of(1839, -0.5, 0, 0), // 1.47e-324
            Arguments.of(1841, -0.5, 0, 0), // 6.55e-325
            Arguments.of(25, -1.5, 1.12615058406163209167415852865e-10, 0),
            Arguments.of(13, -11.5, 8.94160424473182258196042781849e-15, 0),

            // Large s computed using zeta(s, 1 - a)
            // a = -1.5 : zeta(s, 2.5).
            // Verified using WolframAlpha HurwitzZeta without the identity (1 - a)
            // https://www.wolframalpha.com/input?i=HurwitzZeta%28327%2C-1.5%29
            Arguments.of(327, -1.5, 7.47510270791220464622169555878e-131, 0),
            Arguments.of(679, -1.5, 6.29120897946407243273186899998e-271, 0),
            Arguments.of(779, -1.5, 1.01095830534845946704119374011e-310, 0),
            Arguments.of(791, -1.5, 1.69610658558250597457949882757e-315, 0),
            Arguments.of(801, -1.5, 1.77849665908376178480067256262e-319, 0),
            Arguments.of(805, -1.5, 4.55295144725443016908972176031e-321, 0),
            Arguments.of(1099, -1.5, 0, 0), // 4.61e-438
            // a = -123.5 : zeta(s, 124.5)
            // *** Note some error here from mpmath ***
            // *** Mismatch bits: ceil(log2(ulp)) ***
            // *** Repeated below using WolframAlpha ***
            Arguments.of(33, -123.5, 3.19237825509754559298777006642e-69, 3400000), // 22-bits
            Arguments.of(133, -123.5, 3.36756171921967259297351128082e-279, 0),
            Arguments.of(143, -123.5, 3.61590594362399955089661789179e-300, 1750), // 11-bits
            Arguments.of(153, -123.5, 3.89993978254564232599612355239e-321, 0),
            Arguments.of(155, -123.5, 0, 0), // 2.50e-325

            // Computed using https://www.wolframalpha.com/input?i=HurwitzZeta
            Arguments.of(33, -123.5, 3.192378253351764442262703039055841594302566e-69, 0),
            // Note: WolframAlpha is the same using HurwitzZeta(143, 124.5)
            // This case requires using a scaled power in the tail (whole summation)
            // due to (a+n)^-(1+s) being sub-normal.
            Arguments.of(143, -123.5, 3.615905943623999550732130915723115628045590e-300, 0),

            // Copied from a = 124.5 cases
            Arguments.of(147, -123.5, 1.48298686245279271827660650487e-308, 0),
            Arguments.of(149, -123.5, 9.499650589522059112318304392e-313, 0),

            // -------

            // Reference values using Matlab R2026a Symbolic Math Toolbox
            //   vpa(hurwitzZeta(sym(s, 'f'), sym(a, 'f')))
            // Note: The use of 'f' uses the floating-point conversion as N * 2^e
            // where N is the mantissa and e is the exponent.

            // large negative a
            Arguments.of(2, -126.26738, 17.791460939023473438160367544592, 0),
            Arguments.of(3, -126.26738, -50.122585766014119719214807501752, 0),
            Arguments.of(2, -12326.26738, 17.79926823859865890222951879278, 1),
            Arguments.of(3, -12326.26738, -50.122616876585709689690102522698, 0),
            Arguments.of(2, -12326.06738, 223.58067369905866158686397531508, 0),
            Arguments.of(3, -12326.06738, -3268.4963860622797455849029641482, 0),

            // very large negative a (not feasible using a sum of the negative terms)
            Arguments.of(2, -12326165757157.066, 230.08687018030868985783440062902, 0),
            Arguments.of(3, -12326165757157.066, -3414.424543876135548542972795997, 0),

            // Worst case of extreme cancellation. When a -> half-integer
            // and the low series is the full range possible given 0.5+/-2^-b
            Arguments.of(3, -0.5 + 0x1p-54, 0.41439832211715462961751618626938, 0),
            Arguments.of(3, -0.5 + 0x1p-53, 0.41439832211714926143686524195861, 0),
            Arguments.of(3, -0.5 - 0x1p-53, 0.41439832211717073415946901920171, 0),
            Arguments.of(3, -1.5 + 0x1p-52, 0.11810202582084209719727895331851, 0),
            Arguments.of(3, -1.5 - 0x1p-52, 0.11810202582088530580646271524921, 0),
            Arguments.of(3, -3.5 + 0x1p-51, 0.0307784106604705956811455131349, 0),
            Arguments.of(3, -3.5 - 0x1p-51, 0.030778410660557098867785659805987, 0),
            Arguments.of(3, -7.5 + 0x1p-50, 0.0077822558978654467309133842121074, 0),
            Arguments.of(3, -7.5 - 0x1p-50, 0.0077822558980384765932871763254915, 0),
            Arguments.of(3, -15.5 + 0x1p-49, 0.0019512219787038490075394821332913, 0),
            Arguments.of(3, -15.5 - 0x1p-49, 0.0019512219790499147520220765569762, 0),
            Arguments.of(3, -31.5 + 0x1p-48, 0.00048816210820002952659777173497611, 0),
            Arguments.of(3, -31.5 - 0x1p-48, 0.00048816210889216253017517339304625, 0),
            Arguments.of(3, -63.5 + 0x1p-47, 0.00012206286228806035641514115942213, 0),
            Arguments.of(3, -63.5 - 0x1p-47, 0.00012206286367232674283574332594035, 0),
            Arguments.of(3, -127.5 + 0x1p-46, 0.000030517111096024468697569156501944, 0),
            Arguments.of(3, -127.5 - 0x1p-46, 0.000030517113864557336393645086698446, 0),
            Arguments.of(3, -255.5 + 0x1p-45, 0.0000076293626591457113750252501943586, 0),
            Arguments.of(3, -255.5 - 0x1p-45, 0.0000076293681962114704832983463164781, 0),
            Arguments.of(3, -511.5 + 0x1p-44, 0.0000019073412767613820522958206941628, 0),
            Arguments.of(3, -511.5 - 0x1p-44, 0.0000019073523508929061980225612021874, 0),
            Arguments.of(3, -1023.5 + 0x1p-43, 0.00000047682597038482563656491683673066, 0),
            Arguments.of(3, -1023.5 - 0x1p-43, 0.00000047684811864787541032292535846052, 0),
            Arguments.of(3, -2047.5 + 0x1p-42, 0.00000011918713418230492155747450649894, 0),
            Arguments.of(3, -2047.5 - 0x1p-42, 0.00000011923143070840483965021033635372, 0),
            Arguments.of(3, -4095.5 + 0x1p-41, 0.000000029758025417506153675797113283396, 0),
            Arguments.of(3, -4095.5 - 0x1p-41, 0.000000029846618469706082505485151864752, 0),
            // Limit of a double-double summation of opposing terms from 0.5 +/- 2^-40
            Arguments.of(3, -8191.5 + 0x1p-40, 0.0000000073619875169683123404158617742165, 0),
            Arguments.of(3, -8191.5 - 0x1p-40, 0.0000000075391736213681931608483274726682, 0),
            Arguments.of(3, -16383.5 + 0x1p-39, 0.0000000016854590430963498434783046783686, 0),
            Arguments.of(3, -16383.5 - 0x1p-39, 0.0000000020398312518961172746074881289298, 0),
            // 50 ulp = Loss of ~ 6-bits of precision
            Arguments.of(3, -65535.5 + 0x1p-37, -0.00000000059232909577937793989464549320949, 0),
            Arguments.of(3, -65535.5 - 0x1p-37, 0.00000000082515973941969504164666837084463, 0),
            Arguments.of(3, -1048575.5 + 0x1p-33, -0.000000011339455934241697907112956454121, 0),
            Arguments.of(3, -1048575.5 - 0x1p-33, 0.000000011340365428943470628555741917531, 0),

            // Check some larger powers
            Arguments.of(7, -31.5 + 0x1p-48, 0.00000000014222080058151244961535769914742, 0),
            Arguments.of(9, -127.5 - 0x1p-46, 0.00000000026193893956547650755982737479515, 0),
            Arguments.of(5, -1023.5 - 0x1p-43, 0.00000000007309223831961703808738193264568, 0),

            Arguments.of(1.5, 4789, 0.0289021565574206125831622859402, 1),
            Arguments.of(2.345, 12.789, 0.0254390524135780630410689989495, 1),
            Arguments.of(1.345, 12.789, 1.2196695365743183226009917322, 1),
            Arguments.of(1.345, 1278.9, 0.245686258677206272154248922088, 1),
            Arguments.of(4.345, 28697.9, 0.000000000000000366539298049937530342298075842, 1),
            Arguments.of(1.345, 1278562927.9, 0.00209103729002248575358360674936, 0),
            // large s before underflow
            Arguments.of(23.45, 12.789, 0.0000000000000000000000000134520611728481431909082787144, 0),
            Arguments.of(23.45, 1278.9, 8.01879331040682555147782226582e-72, 1),
            Arguments.of(23.45, 1278562927.9, 1.59541010013538002585127908428e-206, 1),
            Arguments.of(234.5, 12.789, 2.7978232305220904641253847499e-260, 0),
            Arguments.of(234.5, 17.89, 1.83179298527375942363255194085e-294, 0),
            // small s -> 1
            Arguments.of(1.00001, 1.789, 99999.7231466483928005870213366, 1),
            Arguments.of(1.00001, 17.89, 99997.1440070978817356297446315, 1),
            Arguments.of(1.00001, 1278562927.89, 99979.0331951153772621341875866, 1),
            Arguments.of(1.0000000000000002, 1.789, 4503599627370495.72314713725233, 0),
            Arguments.of(1.0000000000000002, 1278562927.89, 4503599627370475.03099742881301, 0)
        );
    }

    @ParameterizedTest
    @EnumSource(value = ZetaTestCase.class)
    void testZeta(ZetaTestCase tc) {
        // Cases starting with T_ are for testing code.
        // This can be configured to skip most cases to reduce runtime.
        assertFunction(tc, tc.name().startsWith("T_"), null);
    }

    @ParameterizedTest
    @EnumSource(value = ExtZetaTestCase.class)
    void testZeta(ExtZetaTestCase tc) {
        // Cases starting with T_ are for testing code.
        // This can be configured to skip most cases to reduce runtime.
        assertFunction(tc, 0, tc.name().startsWith("T_"));
    }

    /**
     * Assert the function is close to the expected value.
     *
     * @param fun Function
     * @param x Input value
     * @param y Input value
     * @param expected Expected value
     * @param tolerance the tolerance
     */
    private static void assertClose(DoubleBinaryOperator fun, double x, double y, double expected, int tolerance) {
        final double actual = fun.applyAsDouble(x, y);
        TestUtils.assertEquals(expected, actual, tolerance, null, () -> x + ", " + y);
    }

    /**
     * Assert the function using extended precision.
     *
     * @param tc Test case
     * @param allowSkips Allow skipping cases when not reporting
     * @param resultFilename The name of a result file to record the error
     */
    private static void assertFunction(DoublePrecisionTestCase tc, boolean allowSkips,
        String resultFilename) {
        final TestUtils.ErrorStatistics stats = new TestUtils.ErrorStatistics();
        final boolean noSkips = !allowSkips || doReporting() || resultFilename != null;
        final BooleanSupplier skip = createSkipSupplier(noSkips, 7);
        // Allow optional save of result to file
        try (PrintStream out = resultFilename == null ? null : getPrintStream(resultFilename)) {
            final StringBuilder sb = new StringBuilder(512);
            if (out != null) {
                out.println("s,a,z,ulp_error,m");
            }
            long time = System.nanoTime();
            for (final String filename : tc.getFilenames()) {
                try (DataReader in = new DataReader(filename)) {
                    while (in.next()) {
                        if (skip.getAsBoolean()) {
                            continue;
                        }
                        try {
                            final double s = in.getDouble(0);
                            final double a = in.getDouble(1);
                            final double actual = tc.getFunction().applyAsDouble(s, a);
                            final BigDecimal expected = in.getBigDecimal(2);
                            double e = TestUtils.assertEquals(expected, actual, tc.getTolerance(), stats::add,
                                () -> tc + " s=" + s + ", a=" + a);
                            if (out != null) {
                                // No special formatting required.
                                // Print the raw fields and the error.
                                final String[] fields = in.getFields();
                                sb.setLength(0);
                                for (int i = 0; i < 3; i++) {
                                    sb.append(fields[i]).append(',');
                                }
                                sb.append(e).append(',').append(lastM);
                                out.println(sb.toString());
                            }
                        } catch (final NumberFormatException ex) {
                            Assertions.fail("Failed to load data: " + Arrays.toString(in.getFields()), ex);
                        }
                    }
                } catch (final IOException ex) {
                    Assertions.fail("Failed to load data: " + filename, ex);
                }
            }
            time = System.nanoTime() - time;
            assertRms(tc, stats, time);
        } catch (IOException e) {
            Assertions.fail("Failed to write to result file: " + resultFilename, e);
        }
    }

    /**
     * Assert the function using extended precision.
     *
     * <p>Can ignore results below {@code 2^ignoreB} when {@code ignoreB < 0}.
     * Use to ignore extremely sub-normal results that cannot affect the zeta result
     * in double precision.
     *
     * @param tc Test case
     * @param ignoreB Base-2 exponent of results to ignore
     * @param allowSkips Allow skipping cases when not reporting
     */
    private static void assertFunction(ExtendedPrecisionTestCase tc, int ignoreB, boolean allowSkips) {
        final TestUtils.ErrorStatistics stats = new TestUtils.ErrorStatistics();
        BigDecimal ignore = BigDecimal.ZERO;
        if (ignoreB < 0) {
            ignore = BigDecimal.ONE.divide(new BigDecimal(BigInteger.ONE.shiftLeft(-ignoreB)), MathContext.DECIMAL64);
        }
        final boolean noSkips = !allowSkips || doReporting();
        final BooleanSupplier skip = createSkipSupplier(noSkips, 7);
        final long start = System.nanoTime();
        for (final String filename : tc.getFilenames()) {
            try (DataReader in = new DataReader(filename)) {
                while (in.next()) {
                    if (skip.getAsBoolean()) {
                        continue;
                    }
                    try {
                        final double x = in.getDouble(0);
                        final int s = (int) x;
                        Assertions.assertEquals(x, s, "Expecting integer s");
                        final double y = in.getDouble(1);
                        final BigDecimal expected = in.getBigDecimal(2);
                        if (ignoreB < 0 && expected.compareTo(ignore) <= 0) {
                            continue;
                        }
                        final BigDecimal actual = tc.getFunction().apply(s, y);
                        TestUtils.assertEquals(expected, actual, tc.getTolerance(), tc.scale(),
                            stats::add,
                            () -> tc + " x=" + x + ", y=" + y);
                    } catch (final NumberFormatException ex) {
                        Assertions.fail("Failed to load data: " + Arrays.toString(in.getFields()), ex);
                    }
                }
            } catch (final IOException ex) {
                Assertions.fail("Failed to load data: " + filename, ex);
            }
        }
        assertRms(tc, stats, System.nanoTime() - start);
    }

    /**
     * Assert the Root Mean Square (RMS) error of the function is below the allowed
     * maximum for the specified TestError.
     *
     * @param te Test error
     * @param stats Error statistics
     * @param nanos Duration in nanoseconds
     */
    private static void assertRms(TestError te, TestUtils.ErrorStatistics stats, long nanos) {
        final double rms = stats.getRMS();
        debugRms(te.toString(), stats.getMaxAbs(), rms, stats.getMean(), stats.size(), nanos);
        Assertions.assertTrue(rms <= te.getRmsTolerance(),
            () -> String.format("%s RMS %s < %s", te, rms, te.getRmsTolerance()));
    }

    /**
     * Output the maximum and RMS ulp for the named test. Used for reporting the
     * errors and setting appropriate test tolerances. This is relevant across
     * different JDK implementations where the java.util.Math functions used in
     * BoostGamma may compute to different accuracy.
     *
     * @param name Test name
     * @param maxAbsUlp Maximum |ulp|
     * @param rmsUlp RMS ulp
     * @param meanUlp Mean ulp
     * @param size Number of measurements
     * @param nanos Duration in nanoseconds
     */
    private static void debugRms(String name, double maxAbsUlp, double rmsUlp, double meanUlp,
        int size, long nanos) {
        if (doReporting()) {
            report("%-38s   max %14.6g   RMS %14.6g   mean %14.6g  n %5d  (%.3gms)%n",
                name, maxAbsUlp, rmsUlp, meanUlp, size, nanos * 1e-6);
        }
    }

    /**
     * Check if reporting to stdout. Prints the JDK version on first return of true.
     *
     * <p>Some test fixtures also use the reporting setting to configure random
     * skipping of tests to reduce build times.
     *
     * @return true if reporting
     * @see #reporting
     */
    private static boolean doReporting() {
        if (reporting < 0) {
            return false;
        }
        if (reporting == 0) {
            reporting = 1;
            report("JDK %s %s%n",
                System.getProperty("java.vm.vendor"),
                System.getProperty("java.vm.version")
            );
        }
        return true;
    }

    /**
     * Print the formatted message to stdout.
     *
     * <p><strong>WARNING: Use only in testing code</strong>.
     * No message will be printed when the reporting flag is not configured.
     *
     * @param format Format string.
     * @param args Arguments.
     * @see java.io.PrintStream#printf(String, Object...)
     */
    private static void report(String format, Object... args) {
        if (reporting >= 0) {
            // CHECKSTYLE: stop regexp
            System.out.printf(format, args);
            // CHECKSTYLE: stop regexp
        }
    }

    @Test
    void testSkipSupplier() {
        for (int b = 0; b < 5; b++) {
            // No skips when reporting
            final BooleanSupplier skip = createSkipSupplier(true, b);
            final int interval = 1 << b;
            int count = 0;
            final int total = 10 * interval;
            for (int i = total; i-- >= 0;) {
                if (skip.getAsBoolean()) {
                    count++;
                }
            }
            Assertions.assertEquals(0, count);
        }
        for (int b = 0; b < 5; b++) {
            // Skips
            final BooleanSupplier skip = createSkipSupplier(false, b);
            final int interval = 1 << b;
            int count = 0;
            final int total = 10 * interval;
            for (int i = total; i-- >= 0;) {
                if (skip.getAsBoolean()) {
                    count++;
                }
            }
            // p(not skip) = 1/interval
            final double expected = Math.scalb(1, -b);
            // observed (not skip)
            final double actual = (total - count) / (double) total;
            final int c = count;
            // This may not be exact as the skips are pseudo-random
            Assertions.assertEquals(expected, actual, expected * 0.1,
                () -> String.format("%d : %d/%d", interval, c, total));
        }
    }

    /**
     * Create a supplier that will return {@code false} at regular intervals
     * to signal that a processing a series should not skip the item.
     *
     * <p>If {@code noSkips} is {@code true} this always returns {@code false}.
     *
     * <p>If {@code noSkips} is {@code false} this returns {@code false} at intervals
     * of {@code 2^b}.
     *
     * <p>This allows the class to exercise a subset of test cases for tests
     * that target testing code. The testing code will be run but the time cost
     * is reduced to maintain reasonable build times.
     *
     * @param noSkips the no skips flag
     * @param b the exponent of the interval
     * @return the supplier
     */
    private static BooleanSupplier createSkipSupplier(boolean noSkips, int b) {
        if (noSkips) {
            // Do not skip
            return () -> false;
        }
        // A test that is not reporting. Run at intervals of 2^b
        final int[] counter = {ThreadLocalRandom.current().nextInt()};
        final int mask = (1 << b) - 1;
        return () -> (++counter[0] & mask) != 0;
    }

    /**
     * Get the scaled tail terms used in the zeta evaluation.
     * The terms are scaled to be 2^106 larger.
     *
     * @param m the maximum number of terms
     * @param apn the value of (a+n)
     * @param s the value of s
     * @return the tail terms
     */
    private static double[] getScaledTailTerms(int m, double apn, double s) {
        // (s)_{2k-1} * (a+n)^-(2k-1+s) * B_2k / 2k!
        final double[] t = new double[m];
        final double apn2 = apn * apn;
        double f = 0x1p-916 * s;
        double p = 0x1p1022 * Math.pow(apn, -(1 + s));
        if (p < Double.MIN_NORMAL) {
            p = Math.pow(apn, -(1 + s) * 0.5);
            p = p * 0x1p1022 * p;
            if (p == 0) {
                return t;
            }
        }
        double k2 = 1;
        for (int i = 0; i < m; i++) {
            t[i] = f * p * F[i];
            // If f overflows and p underflow this will create inf*0 = NaN
            if (Double.isNaN(t[i])) {
                t[i] = 0;
                break;
            }
            // Update (a+n)^-(2k-1+s)
            p /= apn2;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
        }
        return t;
    }

    /**
     * Get the tail terms used in the zeta evaluation.
     * Uses a DD fractional representation to prevent underflow.
     *
     * @param m the maximum number of terms
     * @param apn the value of (a+n)
     * @param s the value of s
     * @param et exponent of t (output)
     * @return the tail terms
     */
    private static DD[] getTailTerms(int m, DD apn, int s, int[] et) {
        // (s)_{2k-1} * (a+n)^-(2k-1+s) * B_2k / 2k!
        final DD[] t = new DD[m];
        final long[] e = {0};
        final int[] exp = {0};
        // Rising factorial term : (s)_{2k-1}
        DD f = DD.of(s).frexp(exp);
        long ef = exp[0];
        // Initialise (a+n)^-(2k-1+s) to (a+n)^-(1+s)
        DD p = DDMath.pow(apn, -(1 + s), e);
        long ep = e[0];
        // Divide by (a+n)^2 using multiplication
        final DD apn2 = DDMath.pow(apn, -2, e);
        long eapn = e[0];
        for (int k = 0; k < m; k++) {
            // FDD is finite. f*p in [0.25, 2).
            t[k] = f.multiply(p).multiply(FDD[k]).frexp(exp);
            et[k] = intExponent(ef + ep + exp[0]);
            // p = (a+n)^-(2k-1+s)
            p = p.multiply(apn2).frexp(exp);
            ep += exp[0] + eapn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            // This is cast to double is is fine when s < 2^26.
            f = f.multiply((s + (2L * k) + 1) * (s + (2L * k) + 2)).frexp(exp);
            ef += exp[0];
        }
        return t;
    }

    // TODO - Remove this
//    /**
//     * For a given (s, a), test if convergence is possible in the tail series T within M
//     * terms to double precision using the configured N for the series sum S.
//     *
//     * <p>This uses real {@code s > 1.0}. Note that the tail series eventually diverges as
//     * the Pochammer term (s)_2k-1 (eventually) grows faster than the power and Bernoulli
//     * factor terms reduce. If the tail does not converge within a set M then the result
//     * is error limited by the omitted residual term R.
//     *
//     * <p>As N increases the magnitude of the tail terms reduce compared to the sum S + I.
//     * At the limit when N is infinite then this is simply the infinite zeta series and no
//     * other terms are used. In practice choosing too small an N can give very inaccurate
//     * results if the tail diverges significantly.
//     *
//     * <p>The terms S using N and T using M evaluate the same number of power terms. So N
//     * ~ M is reasonable. However the tail power term can be computed using iterative
//     * updates and is faster. A performant implementation can use {@code N < M}.
//     *
//     * <p>This test computes the point (a+n) where an update to find the magnitude of the
//     * M-th term of the tail series is 1. Any (a+n) after this threshold will have a
//     * converging tail series in the first M terms. Any (a+n) below this will diverge. The
//     * choice of N then requires that terms in the tail series are below the epsilon of
//     * zeta(s, a) before it diverges. Use of the tail can then be stopped before
//     * it diverges.
//     *
//     * <p>This test will assert that the tail for (a+N) is not used before it diverges.
//     * It does not provide errors for an implementation using N and M and the optimal
//     * choice for N may be higher. It does provide the minimum N that should be used
//     * when testing errors of different implementations for a given table size M. Any lower
//     * N can have very high error.
//     *
//     * <p>Note: The test M is used to find where the M-th and (M-1)-th terms are the same.
//     * If a test passes at M=9, all sampled points converged within (M-1) tail terms. So
//     * the test M is 1 higher than the required table size of pre-computed tail factors
//     * for an implementation to avoid tail divergence. If the table size is larger
//     * the implementation may use extra table values for higher (a+n) which can impact
//     * precision.
//     *
//     * <p>Note as s increases the series is strongly convergent for small a, and will
//     * underflow a double for large a. As (a+n) increases the magnitude of the tail
//     * terms reduces compared to the other terms S+I. Thus choosing N requires a value
//     * large enough to avoid excess use of the tail when a is close to 0, and moderate
//     * use when a is small to medium.
//     *
//     * <p>This test is slow when enumerating all supported s. Tests that should pass are
//     * configured to skip samples from a when not in reporting mode. Skipping is periodic
//     * with a random start. When run many times this should visit all points since the
//     * samples are deterministic. Tests that fail as N is too low run without skipping.
//     */
//    @ParameterizedTest
//    @CsvSource({
//        // double precision
//        // Deterministic failure; then enumerate and random sampling when expected to pass.
//        // Use of a smaller tail M requires higher N (this test finds the smallest N
//        // to avoid using diverging tail terms within size M; if the tail is smaller
//        // then N must be higher). Accuracy of the function before/beyond the
//        // point (a+n) where the tail converges within M is not tested.
//
//        // Enumerate and random sampling when expected to pass.
//        // Timings when no reporting is configured are < 1s per test.
//
//        // Note: N ~ M ~ P for P-bits of precision
//        // M=53 requires N >= 6 to converge
//        "1.0000000000000002, 1075, false, -250, 53, 53, -5",
//        "1.0000000000000002, 1075, false, -250, 53, 53, 6",
//        "1.0000000000000002, 1075, false, 250, 53, 53, 6",
//
//        // Smaller M will truncate the tail before it diverges.
//        // May require higher N to converge within the M-th term.
//        "1.0000000000000002, 1075, false, -250, 26, 53, -5",
//        "1.0000000000000002, 1075, false, -250, 26, 53, 6",
//        "1.0000000000000002, 1075, false, 250, 26, 53, 6",
//
//        "1.0000000000000002, 1075, false, -250, 13, 53, -6",
//        "1.0000000000000002, 1075, false, -250, 13, 53, 7",
//        "1.0000000000000002, 1075, false, 250, 13, 53, 7",
//
//        // double-double precision
//        // Note the tail may underflow in this implementation where the DD
//        // implementation would not. The limit is the same as the issue
//        // occurs for s -> 1 and small a.
//        // M=106 requires N >= 18 to converge
//        "1.0000000000000002, 1075, false, -250, 106, 106, -11",
//        "1.0000000000000002, 1075, false, -250, 106, 106, 12",
//        "1.0000000000000002, 1075, false, 250, 106, 106, 12",
//
//        // Smaller M
//        "1.0000000000000002, 1075, false, -250, 53, 106, -11",
//        "1.0000000000000002, 1075, false, -250, 53, 106, 12",
//        "1.0000000000000002, 1075, false, 250, 53, 106, 12",
//        "1.0000000000000002, 1075, false, -250, 26, 106, -12",
//        "1.0000000000000002, 1075, false, -250, 26, 106, 13",
//        "1.0000000000000002, 1075, false, 250, 26, 106, 13",
//    })
//    void testTailConvergence(double minS, double maxS, boolean uniform, int samples, int m,
//            int b, int testN) {
//        // Negative N for tests that should fail
//        final boolean expectedToFail = testN < 0;
//        final boolean noSkips = doReporting() || expectedToFail;
//        final int expectedN = Math.abs(testN);
//
//        Assertions.assertTrue(minS > 1);
//        Assertions.assertTrue(maxS > minS);
//        Assertions.assertTrue(m > 0);
//        Assertions.assertNotEquals(0, samples);
//        Assertions.assertTrue(b > 40);
//        Assertions.assertTrue(expectedN > 3 && expectedN <= 53);
//
//        // Run 1 in 1024 cases
//        final BooleanSupplier skip = createSkipSupplier(noSkips, 10);
//
//        // Evaluate zeta with larger N and support for sub-normal terms.
//        // The purpose of this test is to find the smallest n where the tail can converge before m.
//        // Accuracy when using small N is evaluated elsewhere.
//        final double eps = Math.scalb(1.0, -b);
//        // Limit N ~ M ~ P bits of precision
//        final Context context = Context.of(53, 53).withTailEpsilon(eps);
//
//        // Note negative samples will enumerate and this RNG may not be used.
//        final SplittableRandom rng = new SplittableRandom();
//        DoubleSupplier nextS;
//        if (samples < 0) {
//            nextS = createRange(samples, minS, maxS, uniform);
//        } else {
//            nextS = createSampler(rng, minS, maxS, uniform);
//        }
//
//        // Start below the expected N
//        int n = Math.min(3, expectedN - 1);
//
//        final int inc = samples < 0 ? 1 : -1;
//        for (int i = samples; i != 0; i += inc) {
//            final double s = nextS.getAsDouble();
//            // Tail term: (a+n)^-s * [ (s)_{2k-1} * (a+n)^-(2k-1) * B_2k/2k! ]
//            // For k<=M find (a+n) where the Pochammer term can cause tail divergence as it
//            // is a larger magnitude than the reduction. To compute magnitudes as the
//            // reciprocal this inverts (a+n)^-(2k-1) as (a+n)^(2k-1) and B_2k/2k! as 2k!/B_2k.
//            //
//            // First term: Check if larger
//            // (s)_1 >= (a+n)^(1+s) * 2k!/B_2k
//            // a+n = exp(log(s * B_2k/2k!) / (1+s))
//            //
//            // Successive terms: Check if the pochammer increase is larger than the decrease
//            // (s+2k-1)*(s+2k) >= (a+n)^2 * 2k!/B_2k / 2(k-1)!/B_2(k-1)
//            // a+n = sqrt((s+2k-1)*(s+2k) * 2(k-1)!/B_2(k-1) * B_2k/2k!)
//            //
//            // This only computes the (a+n) for the last term k=M.
//            // An assertion checks terms are descending magnitude up to this point.
//
//            // Cannot overflow for small m
//            final double p = (s + 2 * m - 3) * (s + 2 * m - 2);
//            final double apn = Math.sqrt(p * -RF[m - 2] * F[m - 1]);
//            // This is the value of (a+n) that prevents tail divergence at t[k],
//            // (i.e. where t[k] == t[k-1])
//            double[] t = getScaledTailTerms(m, apn, s);
//            // Sanity check for double precision
//            if (m <= 53) {
//                // Test the final terms are equal
//                if (Math.abs(t[m - 2]) > Double.MIN_NORMAL) {
//                    Assertions.assertEquals(Math.abs(t[m - 2]), Math.abs(t[m - 1]), Math.abs(t[m - 2]) * 1e-15);
//                }
//                // Test the terms reduce in magnitude
//                for (int k = 1; k < m - 1; k++) {
//                    Assertions.assertTrue(Math.abs(t[k - 1]) >= Math.abs(t[k]), "Tail not convergent");
//                }
//            }
//
//            // The aim is to choose large enough n so that t[k] is not required when
//            // the tail will diverge after t[k].
//            // We know (a+n) where the tail converges at M terms for s.
//            // Any a > (a+n) has a converging tail.
//            // Below this the tail will not converge. Ensure the M-th term of the tail
//            // is not required. Compute zeta(s, a) for a -> 0 using a reference
//            // in steps. For each a find the M-th tail term for a+n < ulp(zeta(s, a)).
//
//            // Too many representable values in [0, 0.5): use a spot check
//            for (final double a : new double[] {0x1p-53, 0x1p-43, 0x1p-33, 0.01, 0.25}) {
//                final double z = zeta(s, a, context);
//                if (z > 0x1p106) {
//                    // Zeta cannot be zero with a < 1
//                    continue;
//                }
//                // For comparison to the scaled tail term at the configured epsilon
//                final double stop = z * 0x1p106 * eps;
//                t = getScaledTailTerms(m, a + n, s);
//                while (!canConverge(t, stop)) {
//                    n++;
//                    if (n > expectedN) {
//                        if (expectedToFail) {
//                            return;
//                        }
//                        Assertions.fail(
//                            String.format("Exceeded expected N=%d, s=%s, a=%s, z=%s: t=%s",
//                                expectedN, s, a, z, Arrays.toString(t)));
//                    }
//                    t = getScaledTailTerms(m, a + n, s);
//                }
//            }
//            // Sample a in [0.5, a+n)
//            // Tail may be zero.
//            // First lower a until a tail exists.
//            double ua = apn;
//            while (t[0] == 0) {
//                ua -= 0.5;
//                t = getScaledTailTerms(1, ua, s);
//            }
//            DoubleSupplier nextA;
//            if (samples < 0) {
//                nextA = createRange(samples, 0.5, ua, uniform);
//            } else {
//                nextA = createSampler(rng, 0.5, ua, uniform);
//            }
//            for (int j = samples; j != 0; j += inc) {
//                if (skip.getAsBoolean()) {
//                    continue;
//                }
//                final double a = nextA.getAsDouble();
//                final double z = zeta(s, a, context);
//                if (z > 0x1p106 || z == 0) {
//                    // Scaled tail was non-zero but after scaling the zeta result is zero
//                    continue;
//                }
//                // For comparison to the scaled tail term at the configured epsilon
//                final double stop = z * 0x1p106 * eps;
//                t = getScaledTailTerms(m, a + n, s);
//                while (!canConverge(t, stop)) {
//                    n++;
//                    if (n > expectedN) {
//                        if (expectedToFail) {
//                            return;
//                        }
//                        Assertions.fail(
//                            String.format("Exceeded expected N=%d, s=%s, a=%s, z=%s: t=%s",
//                                expectedN, s, a, z, Arrays.toString(t)));
//                    }
//                    t = getScaledTailTerms(m, a + n, s);
//                }
//            }
//        }
//        Assertions.assertFalse(expectedToFail,
//            () -> String.format("Should fail M=%d, maxS=%s, N=%d", m, maxS, testN));
//    }

    /**
     * Test to identify the maximum (a+n) value for s where the tail can diverge within M terms.
     *
     * <P>When the tail is longer it can diverge (terms get larger). This occurs when s is
     * close to or above the (a+n) value. Large divergence creates error in the zeta evaluation.
     *
     * <p>This is a simple test that computes the value (a+n) where a tail of length M
     * is equal in the last two terms for chosen s. Any (a+n) above this will not
     * create large errors. The tail is asserted to converge up to the M-th term.
     *
     * <p>This test shows the N required to avoid divergence for the entire zeta domain
     * from a=0 when limited to a non-zero tail of size M. It does not evaluate whether
     * the N is optimum for speed or accuracy. A lower N and/or M will be faster but
     * divergence may occur in the tail for some part of the domain.
     *
     * <p>Note that the highest (a+n) values are for large s when the zeta result is non-zero.
     * This can use s in (1, MAX_VALUE] but a more realistic range before the tail terms
     * are zero is 150 to 250. When the tail is zero the zeta result may not be zero, but
     * tail divergence is not an issue so the test just evaluates the tail which is very fast.
     */
    @ParameterizedTest
    @CsvSource({
        // A longer tail requires a higher (a+n) to ensure it never diverges.

        // These parameters can be printed by this test when the reporting is on.
        // The report add 0.1 to the maximum N.

        // For the minimum S this is the N required so that a tail of length M does not diverge.
        // It does not evaluate whether a short tail will be accurate.
        "1.0000000000000002, 1.0000000000000002, false, 1, 53, 16.73", // s=1.0000000000000002
        "1.0000000000000002, 1.0000000000000002, false, 1, 26, 8.14", // s=1.0000000000000002
        "1.0000000000000002, 1.0000000000000002, false, 1, 18, 5.59", // s=1.0000000000000002
        "1.0000000000000002, 1.0000000000000002, false, 1, 15, 4.64", // s=1.0000000000000002
        "1.0000000000000002, 1.0000000000000002, false, 1, 12, 3.68", // s=1.0000000000000002
        // For the maximum S (found by search of large s where tail is supported)
        // this is the N required so that a tail of length M does not diverge.
        // Note that for M=53 the (a+n) < 49. Thus for the entire zeta domain from a=0
        // a tail of length 53 will never diverge using N=49.
        // Use of N=49, M=53 in the Context will not have large errors.
        "150.0, 250.0, false, 2500, 53, 49.87", // s=209.18367346935452
        "150.0, 250.0, false, 2500, 39, 46.10", // s=213.50540216082868
        "150.0, 250.0, false, 2500, 26, 42.67", // s=217.94717887151046
        "150.0, 250.0, false, 2500, 18, 40.59", // s=220.9083633452983
        "150.0, 250.0, false, 2500, 15, 39.82", // s=222.06882753097193
        "150.0, 250.0, false, 2500, 12, 39.06", // s=223.2693077230481
    })
    void testTailConvergenceLimit(double minS, double maxS, boolean uniform, int samples, int m,
            double testN) {
        Assertions.assertTrue(minS > 1);
        Assertions.assertTrue(maxS >= minS);
        Assertions.assertTrue(m > 0);
        Assertions.assertNotEquals(0, samples);

        final DoubleSupplier nextS = minS == maxS ? () -> minS : createRange(samples, minS, maxS, uniform);

        double maxapn = 0;
        double maxs = 1;

        for (int i = 0; i < samples; i++) {
            final double s = nextS.getAsDouble();
            // Tail term: (a+n)^-s * [ (s)_{2k-1} * (a+n)^-(2k-1) * B_2k/2k! ]
            // For k<=M find (a+n) where the Pochammer term can cause tail divergence as it
            // is a larger magnitude than the reduction. To compute magnitudes as the
            // reciprocal this inverts (a+n)^-(2k-1) as (a+n)^(2k-1) and B_2k/2k! as 2k!/B_2k.
            //
            // First term: Check if larger
            // (s)_1 >= (a+n)^(1+s) * 2k!/B_2k
            // a+n = exp(log(s * B_2k/2k!) / (1+s))
            //
            // Successive terms: Check if the pochammer increase is larger than the decrease
            // (s+2k-1)*(s+2k) >= (a+n)^2 * 2k!/B_2k / 2(k-1)!/B_2(k-1)
            // a+n = sqrt((s+2k-1)*(s+2k) * 2(k-1)!/B_2(k-1) * B_2k/2k!)
            //
            // This only computes the (a+n) for the last term k=M.
            // An assertion checks terms are descending magnitude up to this point.

            // Cannot overflow for small m
            final double p = (s + 2 * m - 3) * (s + 2 * m - 2);
            final double apn = Math.sqrt(p * -RF[m - 2] * F[m - 1]);
            // This is the value of (a+n) that prevents tail divergence at t[k],
            // (i.e. where t[k] == t[k-1])
            double[] t = getScaledTailTerms(m, apn, s);
            // Sanity check for double precision
            if (m <= 53) {
                // Test the final terms are equal
                if (Math.abs(t[m - 2]) > Double.MIN_NORMAL) {
                    Assertions.assertEquals(Math.abs(t[m - 2]), Math.abs(t[m - 1]), Math.abs(t[m - 2]) * 1e-15);
                }
                // Test the terms reduce in magnitude
                for (int k = 1; k < m - 1; k++) {
                    Assertions.assertTrue(Math.abs(t[k - 1]) >= Math.abs(t[k]), "Tail not convergent");
                }
            }

            // We know (a+n) after which any zeta evaluation using a tail of M will not diverge.
            // Check if the tail can be used. We could evaluate zeta here with sub-normal support.
            // Instead we check if the tail terms can be used.
            if (t[0] == 0) {
                // With ascending s no more tail computation is possible
                // (all tail terms for larger s and (a+n) will also be zero).
                // If the tail * 2^106 is zero then the zeta evaluation is
                // zero.
                break;
            }

            Assertions.assertTrue(testN > apn,
                () -> String.format("M=%d, s=%s, (a+n)=%s, N=%s", m, s, apn, testN));

            if (maxapn < apn) {
                maxapn = apn;
                maxs = s;
            }
        }
        // Print out the format for the arguments. Add 0.1 to the a+n so the test passes.
        report("\"%s, %s, %b, %d, %d, %.2f\", // s=%s%n",
            minS, maxS, uniform, samples, m, maxapn + 0.1, maxs);
    }

    // TODO:
    // Create the function and plot the score over the domain of zeta.
    // Can use the same plotting as for the planned error plot.
    // Should be easy as getScaledTailTerms accepts real (a+n).
    // Effectively evaluate zeta. Then bracket getScaledTail terms
    // with variable n until it can converge. It will show the hotspot
    // in the domain for convergence issues.

    // Function goals
    //
    // [1] We want the tail to always converge within M when used.
    // Tested in testTailConvergence2.
    //
    // [2] We want to never rely on requiring an accurate tail T.
    // TODO: Untested
    //
    // The first tail term should always be below the magnitude of zeta
    // for all s and a:
    // (a+n)^-(1+s) * s / 12 < zeta(s, a)
    //
    // Possible approximation that could be solved:
    // (a+n)^-(1+s) * s / 12 < a^-s + (a+1)^(1-s)/(s-1)
    // Convergence is slowest at large a=2^53
    // Maximum s to support at a=2^53
    // (2^53)^-s = 2^-1075 (effectively zero)
    // s = 1075 / 53 = 20.28
    // As s reduces the tail is more important.
    // find n where (a+n)^(1-s)/(s-1) = 0.5 * zeta(s, a) and the tail converges.
    //
    // slowest convergence s = 1 + 2^-52
    // fastest convergence s = 20.28
    // critical convergence I ~ 0.5 * zeta(s, a)
    // find n where (a+n)^(1-s)/(s-1) = 0.5 * zeta(s, a) and the tail converges in m terms at a=0


    // TODO: Write an equivalent for the DD version.
    // To maximise (a+n) requires fixed a.
    // But s is not univariate for an optimiser. We can enumerate s.

    /**
     * A test to search for the maximum n required to converge the tail.
     *
     * <ul>
     * <li>When a is much larger than s then the tail cannot realistically diverge
     * with a fixed length of M terms.
     * <li>When s is much larger than a then the series converges within N terms
     * and there is no tail.
     * <li>Between these extremes is a point where the tail can diverge unless the
     * series is evaluated to a higher number of terms N. This increases the current
     * sum's magnitude in relation to the magnitude of the tail.
     *
     * <p>To search the entire parameter space (s, a) requires use of an optimiser for
     * a function of (s, a) to maximise n. Here we sample random a within the domain
     * and optimise a univariate function of s to n as a real value. This function
     * is itself a binary search within a supported n=[1, maxN] for the smallest
     * n where the tail can converge.
     *
     * <p>The highest n occurs with larger s where {@code s > a}. This is limited by
     * the tail becoming zero as zeta approaches zero. So the bound on a only
     * requires a small range.
     *
     * <p>Note this test provides a minimum bound on N for a zeta evaluation using M
     * tail terms. It does not assess the accuracy.
     */
    @ParameterizedTest
    @CsvSource({
        // Note:
        // Worst case for a=0 is to use N where convergence is possible for a=1
        // and use (N+1). However when a -> 0 the usage of the tail reduces as
        // a^-s increasingly dominates the result when both (a+1)^-s can be added and
        // I = (a+1)^-s * (a+1) / (s-1) is not significant.
        // The largest tail usage is for large a.
        // No requirement to test down to a=0.

        // These parameters can be printed by this test when the reporting is on.
        // The report add 0.1 to the maximum N.

        // TODO - Add some tests that fail.
        // Increase speed when not reporting with lower epsilon.
        // Repeat with much higher samples.
        // Cross-verify with real data from the domain.

        "25.0, 130.0, false, 250, 8, 53, 9.78", // (159.0228355327598, 95.58232931726678) = 1.192e-283 * 2^-106 = 1.470e-315 : t_0=2.656e-291 (z * 2^-25.42)
        "25.0, 130.0, false, 250, 10, 53, 7.56", // (172.10361659747323, 66.5783132530103) = 1.386e-282 * 2^-106 = 1.709e-314 : t_0=2.874e-291 (z * 2^-28.85)
        "25.0, 130.0, false, 250, 12, 53, 6.60", // (182.3319254903897, 51.732931726907005) = 2.807e-281 * 2^-106 = 3.460e-313 : t_0=6.090e-291 (z * 2^-32.10)
        "25.0, 130.0, false, 250, 15, 53, 6.02", // (194.34071488630175, 39.43775100401567) = 5.779e-279 * 2^-106 = 7.123e-311 : t_0=3.246e-291 (z * 2^-40.70)
        "25.0, 130.0, false, 250, 26, 53, 5.86", // (203.45493627683237, 32.502008032128245) = 2.016e-276 * 2^-106 = 2.485e-308 : t_0=3.398e-291 (z * 2^-49.08)
        "25.0, 130.0, false, 250, 53, 53, 5.86", // (203.45493627683217, 32.502008032128245) = 2.016e-276 * 2^-106 = 2.485e-308 : t_0=3.398e-291 (z * 2^-49.08)
        "25.0, 130.0, false, 250, 15, 54, 6.16", // (193.47877526730463, 40.06827309236907) = 6.382e-279 * 2^-106 = 7.867e-311 : t_0=3.231e-291 (z * 2^-40.85)
        "25.0, 130.0, false, 250, 26, 54, 5.88", // (204.29789206451426, 31.62048192771072) = 2.948e-275 * 2^-106 = 3.634e-307 : t_0=1.365e-290 (z * 2^-50.94)
        "25.0, 130.0, false, 250, 15, 55, 6.29", // (189.4338618296996, 43.85140562248948) = 7.379e-280 * 2^-106 = 9.096e-312 : t_0=3.164e-291 (z * 2^-37.76)
        "25.0, 130.0, false, 250, 26, 55, 6.08", // (203.96998643693334, 31.93574297188742) = 1.210e-275 * 2^-106 = 1.491e-307 : t_0=3.407e-291 (z * 2^-51.66)
        // Critical size M=15 to converge tail at 53-55 bits precision is ~ n=6
        // Check higher M
        "25.0, 130.0, false, 250, 17, 54, 5.89", // (198.60538165878268, 35.02409638554185) = 1.553e-275 * 2^-106 = 1.914e-307 : t_0=3.980e-289 (z * 2^-45.15)
        "25.0, 130.0, false, 250, 18, 54, 5.87", // (201.4685639548269, 33.13253012048165) = 4.234e-275 * 2^-106 = 5.219e-307 : t_0=1.615e-289 (z * 2^-47.90)
        "25.0, 130.0, false, 250, 17, 55, 6.10", // (193.89624860931625, 39.75301204819237) = 6.322e-279 * 2^-106 = 7.793e-311 : t_0=3.238e-291 (z * 2^-40.83)
        "25.0, 130.0, false, 250, 18, 55, 6.10", // (201.32469559214005, 33.76305220883505) = 1.575e-276 * 2^-106 = 1.942e-308 : t_0=3.362e-291 (z * 2^-48.73)
        // Converge tail beyond 10-bits of the result ~ n=7
        "25.0, 130.0, false, 250, 18, 63, 7.12", // (192.77329062798356, 39.75301204819237) = 3.954e-277 * 2^-106 = 4.873e-309 : t_0=6.439e-291 (z * 2^-45.80)
        "25.0, 130.0, false, 250, 26, 63, 6.96", // (202.91004799242464, 31.77811244979907) = 1.297e-273 * 2^-106 = 1.599e-305 : t_0=3.389e-291 (z * 2^-58.41)

        // double-double precision
        // This is optimistic as the tail terms may become zero due to underflow.
        // The DD implementation does not underflow.
        "25.0, 130.0, false, 250, 26, 106, 12.47", // (182.8974681390024, 44.79718875501958) = 8.079e-271 * 2^-106 = 9.958e-303 : t_0=6.109e-291 (z * 2^-66.84)
        "25.0, 130.0, false, 250, 53, 106, 11.70", // (197.23117481581144, 31.30522088353402) = 8.469e-264 * 2^-106 = 1.044e-295 : t_0=6.588e-291 (z * 2^-90.05)
        "25.0, 130.0, false, 250, 26, 108, 12.86", // (179.34483214935423, 49.5261044176701) = 9.155e-273 * 2^-106 = 1.128e-304 : t_0=2.995e-291 (z * 2^-61.41)
        "25.0, 130.0, false, 250, 53, 108, 11.90", // (195.92410211834715, 32.186746987951544) = 3.311e-264 * 2^-106 = 4.081e-296 : t_0=3.272e-291 (z * 2^-89.71)
    })
    void testTailConvergence(double minA, double maxA, boolean uniform, int samples, int m,
            int b, double testN) {
        // Negative N for tests that should fail
        final boolean expectedToFail = testN < 0;
        final double expectedN = Math.abs(testN);

        Assertions.assertTrue(minA > 0x1p-53);
        Assertions.assertTrue(maxA > minA);
        Assertions.assertTrue(m > 0);
        Assertions.assertNotEquals(0, samples);
        Assertions.assertTrue(b > 40);
        Assertions.assertTrue(expectedN > 3 && expectedN <= 53);

        // Evaluate zeta with larger N and support for sub-normal terms.
        // The purpose of this test is to find the smallest n where the tail can converge before m.
        // Accuracy when using small N is evaluated elsewhere.
        final double eps = Math.scalb(1.0, -b);
        final Context context = Context.DOUBLE_REF
            .withEpsilon(eps)
            .addOptions(Context.SCALED_RESULT);

        double n = 0;
        double[] max = {};

        // For selected a in the supported range: [2^-1023, 2^53]
        // Search for s where (a+n) to converge the tail is the maximum.
        // Caution:
        // Lowering the convergence tolerance increases speed.
        // The result is only required to find an integer so a few digits
        // is sufficient. When the result is close to an integer then the
        // tolerance can be raised to check to higher accuracy.

        // when not reporting change epsilon to 1e-3
        final double relativeTolerance = 1e-15;
        final BrentOptimizer optimizer = new BrentOptimizer(relativeTolerance, Math.ulp(1.0));

        // The maximum n is always at the upper limit of a.
        // So run in ascending order. When there is nothing to optimise then stop.
        final DoubleSupplier nextA = createRange(samples, minA, maxA, uniform);
        for (int i = 0; i < samples; i++) {
            final double a = nextA.getAsDouble();

            // Find the largest n for (a+n) used in tail convergence at point a.
            // Function to maximise is s -> n.
            final DoubleUnaryOperator fun = s -> {
                final double z = zeta(s, a, context);
                // Note the result is scaled 2^106
                if (z > LARGE_T0 * UPSCALE_106 || z == 0) {
                    // Only occurs with a < 0. s is too large. No tail usage.
                    return 0;
                }
                // Find a+n required for tail convergence.
                // Search n in [1, N] where N is used by the zeta reference.
                final double stop = z * eps;
                double l = 1;
                double u = context.getN();
                // Converge to higher accuracy than the solver is using
                final double tol = relativeTolerance * 0.5;
                while ((u - l) > u * tol) {
                    double middle = (l + u) * 0.5;
                    double[] t = getScaledTailTerms(m, a + middle, s);
                    if (canConverge(t, stop)) {
                        u = middle;
                    } else {
                        l = middle;
                    }
                    // Safety check to avoid infinite loop when the relative tolerance
                    // is too low. Stop when 1 ulp apart.
                    if (Math.nextUp(l) == u) {
                        break;
                    }
                }
                // u is the smallest n that could converge for this s
                return u;
            };

            // Optimization start point:
            //
            // Find where a tail can be used with n=1.
            // We want the tail to overlap with the sum (S+I) at the configured epsilon.
            //
            // This function knows the highest use of the tail is at larger a
            // where s > a and zeta(s, a) -> 0. So the bounds on s can be fixed
            // where zeta(s, a) is definitely zero and s is much larger.
            // The start s can be set to ~a and increased until the tail is used.
            // This can be input to the optimiser.
            // Note:
            // The I term represents all terms to infinity that are not summed in n-1 terms.
            // I ~ zeta(s, a+n)
            // I = (a+n)^-s * (a+n) / (s-1)
            // First term of zeta(s, a+n) = (a+n)^-s.
            // When the I term equals this value then the series converges in n terms.
            // (a+n)^-s ~ (a+n)^-s * (a+n) / (s-1)
            // (a+n) = (s-1)
            // For n=1: (a+1) = (s-1) : s = a+2
            final double minS = Math.nextUp(1.0);
            final double maxS = 400;
            double startS = a + 2;

            // Increase s until the function is higher
            double n0 = fun.applyAsDouble(startS);
            while (n0 < 1.1) {
                startS += 1;
                if (startS > maxS) {
                    Assertions.assertNotEquals(0, n,
                        () -> String.format("Nothing to optimize at min a: M=%d, a=%s%n", m, minA));
                    break;
                }
                n0 = fun.applyAsDouble(startS);
            }
            if (startS > maxS) {
                // Reached the upper limit
                break;
            }

            // TODO - Investigate this
            // Highest tail usage is where the result ~ MIN_NORMAL.
            // Need better resolution at high a.
            // This is where the scaled result in double-double precision does not
            // overlap MIN_NORMAL. Is this coincidence.
            // A plot of the used M across the domain is required to see hot spots.

            final UnivariatePointValuePair point = optimizer.optimize(GoalType.MAXIMIZE, new SearchInterval(minS, maxS, startS),
                 new UnivariateObjectiveFunction(fun::applyAsDouble), new MaxEval(10000));
            final double s = point.getPoint();
            n0 = point.getValue();
            // Check that use of n above this threshold matches the expected result
            final double z = zeta(s, a, context);
            // Add some tolerance
            final int nn = (int) Math.ceil(n0 + 0.05);
            final double zn = zeta(s, a, context.withN(nn));
            // When n is very low the result can differ a few ULP
            final int ulp = nn <= 5 ?
                10 :
                // When the scaled result is < 2^53 * sub-normal the result is zero.
                // Check the computation well below this level.
                (z < 0x1p10 * Double.MIN_NORMAL) ? 5 : 1;
            TestUtils.assertEquals(z, zn, ulp,
                String.format("zeta(%s, %s) (n=%d, m=%d) : %s : %s != %s (n=%d; %s)",
                    s, a, context.getN(), m, z * DOWNSCALE_106, z, zn, nn, n0));
            if (n < n0) {
                n = n0;
                max = new double[] {s, a, z, n0};
                if (n > expectedN) {
                    if (expectedToFail) {
                        return;
                    }
                    Assertions.fail(
                        String.format("A=[%s, %s) M=%d : Exceeded expected N=%s < %s at zeta(%s, %s)",
                            minA, maxA, m, expectedN, n, s, a));
                }
            }
            // TODO - remove this
//            report("%25s  %25s, %25s + %25s = %s * 2^106 = %s%n", startS, s, a, n0, z * DOWNSCALE_106, z);
        }

        if (doReporting()) {
            final double s = max[0];
            final double a = max[1];
            final double z = max[2];
            // Here the M-th scaled tail term is < z * eps
            // Report magnitude difference in bits for the first term:
            // log2(z / t_0)
            // log2(z) - log2(t_0)
            final double[] t = getScaledTailTerms(1, a + n, s);
            final double bits = (Math.log(z) - Math.log(t[0])) / Math.log(2);
            // Print out the format for the arguments. Add 0.1 to the N so the test passes.
            report("\"%s, %s, %b, %d, %d, %d, %.2f\", // (%s, %s) = %.4g * 2^-106 = %.4g : t_0=%.4g (z * 2^-%.2f)%n",
                minA, maxA, uniform, samples, m, b, n + 0.1,
                s, a, z, z * DOWNSCALE_106, t[0], bits);
        }
        if (expectedToFail) {
            Assertions.fail(String.format("A=[%s, %s) M=%d : Expected to fail N=%s > %s at zeta(%s, %s)",
                    minA, maxA, m, expectedN, n, max[0], max[1]));
        }
    }

    @Test
    void test() {
        // TODO - add these to spot checks
//        double s = 1121.1037200772555, a = 9.007199254740992E15, expected = 0;
//        double s = 1121.1037200772555, a = Math.nextUp(1.0), expected = 0.999999999999751064967395493924;
//        double s=1501.7021967661672, a=1.535965979848382, expected =1.2892791106364060286402216376328398051685475849534E-280;
//        double s=798.4098080701187, a=2.44681658669993, expected =5.4552937300588262652119314762621199432693920384755E-311;
//        double  s=1.0046876033700256, a=2.774626370847173E13, expected=184.51544775024293707146458999993276990921062542111;
//        double s=20.23754877788927, a=9.134550490612894E15, expected=4.6919362125812846992028258016880507980787512899492E-309;
//          double s=180.26951483738634, a=48.74884443673853, expected=5.2823717105788599902950123870564151804047370926344E-305;
//          double  s=1.0427289765022538, a=3.94996781319737E-16, expected=11523336262829300.1663588738791442776;

          double s=27316.932399926933, a=1.000594477194544, expected=8.9014344678386187888608455646373495579064616407205E-8;


//        double s = Math.nextUp(1.0), a = 9.007199254740992E15, expected = 4503599627370459.5;
//        TestUtils.assertEquals(expected, zeta(s, a, Context.DOUBLE), 0);
        TestUtils.assertEquals(expected, zeta(s, a, Context.DOUBLE), 1);
//        TestUtils.assertEquals(expected * UPSCALE_106, zeta(s, a,
//            Context.DOUBLE.addOptions(Context.SCALED_RESULT)), 200);
    }

    /**
     * Test if the tail can converge on the provided stop magnitude.
     *
     * @param t Tail terms.
     * @param stop Stop magnitude
     * @return true if the tail converges
     */
    private static boolean canConverge(double[] t, double stop) {
        for (int k = 0; k < t.length; k++) {
            // Tail terms alternate sign
            if (Math.abs(t[k]) <= stop) {
                return true;
            }
        }
        return false;
    }

    /**
     * For a given (s, a), test if convergence is possible in the tail series T within M
     * terms to double-double precision using the configured N for the series sum S.
     *
     * <p>This id the DD equivalent of
     * {@link #testTailConvergence(double, double, boolean, int, int, int, int)}.
     *
     * <p>The test purpose is to assert that an implementation is possible using N and M.
     * It does not assess the error of such an implementation. In practice too low an N
     * can have significant error if the tail diverges.
     *
     * <p>This enumerates integer {@code s > 1}.
     *
     * <p>This test is slow when enumerating all supported s. Tests that should pass are
     * configured to skip samples from a when not in reporting mode. Skipping is periodic
     * with a random start. When run many times this should visit all points since the
     * samples are deterministic. Tests that fail as N is too low run without skipping.
     */
    @ParameterizedTest
    @CsvSource({
        // Limit of s.
        // Set using 2^-s where the full precision of the result of any term (a+k)^-s
        // is below the min double value:
        // double 1075
        // double-double 1075+53 = 1128
        // Maximum s used in the total cancellation case (only requires double precision)
        // 1837

        // Enumerate and random sampling when expected to pass.
        // Timings when no reporting is configured are < 1s per test.

        // Note: N ~ M ~ P for P-bits of precision
        // For double precision
        "1837, false, -250, 53, 53, -5",
        "1837, false, -250, 53, 53, 6",
        "1837, false, 250, 53, 53, 6",
        "1837, false, -250, 13, 53, -6",
        "1837, false, -250, 13, 53, 7",
        "1837, false, 250, 13, 53, 7",
        // double-double precision
        "1128, false, -250, 106, 106, -11",
        "1128, false, -250, 106, 106, 12",
        "1128, false, 250, 106, 106, 12",
        "1128, false, -250, 26, 106, -12",
        "1128, false, -250, 26, 106, 13",
        "1128, false, 250, 26, 106, 13",
        // quad-double precision (for the BigDecimal implementation testing)
        "1128, false, -250, 106, 212, -23",
        "1128, false, -250, 106, 212, 24",
        "1128, false, 250, 106, 212, 24",
        "1128, false, -250, 53, 212, -24",
        "1128, false, -250, 53, 212, 25",
        "1128, false, 250, 53, 212, 25",
    })
    void testTailConvergenceDD(int maxS, boolean uniform, int samples, int m,
            int b, int testN) {
        // Negative N for tests that should fail (run all cases)
        final boolean expectedToFail = testN < 0;
        final boolean noSkips = doReporting() || expectedToFail;
        final int expectedN = Math.abs(testN);

        Assertions.assertTrue(maxS > 1);
        Assertions.assertTrue(m > 0);
        Assertions.assertNotEquals(0, samples);
        Assertions.assertTrue(b > 40);
        Assertions.assertTrue(expectedN > 3 && expectedN <= 106);

        // Run 1 in 1024 cases
        final BooleanSupplier skip = createSkipSupplier(noSkips, 10);

        // Evaluate zeta with larger N and support for sub-normal terms.
        // The purpose of this test is to find the smallest n where the tail can converge before m.
        // Accuracy when using small N is evaluated elsewhere.
        final double eps = Math.scalb(1.0, -b);
        // Limit N ~ M ~ P bits of precision
        final Context context = Context.of(Math.min(106, expectedN * 3), 106).withEpsilon(eps);

        // Note negative samples will enumerate and this RNG may not be used.
        final SplittableRandom rng = new SplittableRandom();

        // Start below the expected N
        int n = Math.min(3, expectedN - 1);
        // For fractional representations
        final int[] et = new int[m];
        final int[] exp = {0};

        final int inc = samples < 0 ? 1 : -1;
        for (int s = 2; s <= maxS; s++) {
            // This follows the logic of the double precision test using DD arithmetic.
            // The computation of (a+n) that creates equal magnitude updates is limited
            // to the M-th term.

            // Cannot overflow for small m
            final long p = (s + 2L * m - 3) * (s + 2L * m - 2);
            // Note F factors alternate signs so 1 negation creates a positive.
            final DD apn = FDD[m - 1].divide(FDD[m - 2].negate()).multiply(DD.of(p)).sqrt();
            // Check against the simpler double version
            TestUtils.assertEquals(Math.sqrt(p * -RF[m - 2] * F[m - 1]), apn.hi(), 10);

            // This is the value of (a+n) that prevents tail divergence at t[k]
            // (the value of (a+n) where t[k] == t[k-1])
            DD[] t = getTailTerms(m, apn, s, et);
            // Sanity check against the simpler double version ...
            if (s < 100) {
                final int ss = s;
                final double[] tt = getScaledTailTerms(Math.min(m, 3), apn.hi(), s);
                for (int k = 0; k < tt.length; k++) {
                    if (Math.abs(tt[k]) < Double.MIN_NORMAL) {
                        break;
                    }
                    // The double tail terms are scaled by 2^106
                    TestUtils.assertEquals(tt[k], Math.scalb(t[k].hi(), 106 + et[k]), 100, null,
                        () -> String.format("Mismatch with double version: s=%d, a=%s", ss, apn.hi()));
                }
            }

            // Test the final terms are equal
            Assertions.assertEquals(Math.abs(t[m - 2].hi()), Math.abs(t[m - 1].hi()));
            Assertions.assertEquals(et[m - 2], et[m - 1]);
            // Test the terms reduce in magnitude
            for (int k = 1; k < m - 1; k++) {
                Assertions.assertTrue(et[k - 1] > et[k] ||
                    Math.abs(t[k - 1].hi()) >= Math.abs(t[k].hi()), "Tail not convergent");
            }

            // Too many representable values in [0, 0.5): use a spot check
            for (final double a : new double[] {0x1p-53, 0x1p-43, 0x1p-33, 0.01, 0.25}) {
                final DD z = zeta(s, DD.of(a), null, exp, context);
                // For comparison to the scaled tail term at the configured epsilon
                final int stop = exp[0] - b;
                t = getTailTerms(m, DD.ofSum(a, n), s, et);
                while (!canConverge(et, stop)) {
                    n++;
                    if (n > expectedN) {
                        if (expectedToFail) {
                            return;
                        }
                        Assertions.fail(
                            String.format("Exceeded expected N=%d, s=%s, a=%s, z=%s * 2^%d: t=%s",
                                expectedN, s, a, z, exp[0], Arrays.toString(et)));
                    }
                    t = getTailTerms(m, DD.ofSum(a, n), s, et);
                }
            }
            // Sample a in [0.5, a+n)
            // A tail will always be non-zero due to use of the fractional representation
            // Lower a until a 106-bit double-double cannot be added to Double.MIN_VALUE.
            double ua = apn.hi();
            while (et[0] < -1075 - 106) {
                ua -= 0.5;
                t = getTailTerms(m, DD.of(ua), s, et);
            }
            DoubleSupplier nextA;
            if (samples < 0) {
                nextA = createRange(samples, 0.5, ua, uniform);
            } else {
                nextA = createSampler(rng, 0.5, ua, uniform);
            }
            for (int j = samples; j != 0; j += inc) {
                if (skip.getAsBoolean()) {
                    continue;
                }
                final double a = nextA.getAsDouble();
                final DD z = zeta(s, DD.of(a), null, exp, context);
                // For comparison to the scaled tail term at the configured epsilon
                final int stop = exp[0] - b;
                t = getTailTerms(m, DD.ofSum(a, n), s, et);
                while (!canConverge(et, stop)) {
                    n++;
                    if (n > expectedN) {
                        if (expectedToFail) {
                            return;
                        }
                        Assertions.fail(
                            String.format("Exceeded expected N=%d, s=%s, a=%s, z=%s * 2^%d: t=%s",
                                expectedN, s, a, z, exp[0], Arrays.toString(et)));
                    }
                    t = getTailTerms(m, DD.ofSum(a, n), s, et);
                }
            }
        }
        Assertions.assertFalse(expectedToFail,
            () -> String.format("Should fail M=%d, maxS=%s, N=%d", m, maxS, testN));
    }

    /**
     * Test if the tail can converge on the provided stop magnitude.
     *
     * @param t Exponent of the tail terms.
     * @param stop Exponent of the stop magnitude
     * @return true if the tail converges
     */
    private static boolean canConverge(int[] t, int stop) {
        for (int k = 0; k < t.length; k++) {
            if (t[k] < stop) {
                return true;
            }
        }
        return false;
    }

    /**
     * Assert that the tail fails to converge. This uses failure cases accrued during testing.
     * It ensures that cases that failed the test implementation can be reviewed if
     * the implementation changes.
     */
    @ParameterizedTest
    @MethodSource
    void testFailedTailConvergence(double s, double a, int n, int m) {
        // Evaluate zeta with large N and support for sub-normal terms.
        // The purpose of this test is to examine if the tail can converge with
        // a much smaller N. Accuracy when using small N is evaluated elsewhere.
        final Context context = Context.DOUBLE_REF;
        final double z = zeta(s, a, context);
        final double[] t = getScaledTailTerms(m, a + n, s);
        // Convergence epsilon is 2^-53, tail terms are scaled 2^106.
        final double stop = z * 0x1p53;
        Assertions.assertFalse(canConverge(t, stop), () -> String.format(
            "Unexpected tail convergence: zeta(%s, %s) %s n=%d, m=%d. Tail = 2^-106 * %s",
                s, a, z, n, m, Arrays.toString(t)));
    }

    static Stream<Arguments> testFailedTailConvergence() {
        return Stream.of(
            Arguments.of(30.784011044176562, 7.205411373936762, 6, 13),
            Arguments.of(3.326653306613046, 1.0, 6, 12),
            Arguments.of(4.0, 1.0, 6, 12),
            Arguments.of(16.0, 2.2244488977955688, 7, 7),
            Arguments.of(16.0, 4.777555110220352, 7, 10),
            Arguments.of(64.0, 18.43687374749453, 7, 10),
            Arguments.of(5.977311652199914, 3.3806723392907463, 7, 9),
            Arguments.of(2.910160291044077, 1.1639528045623206, 7, 9),
            Arguments.of(10.784409764595118, 3.442693930852058, 8, 8),
            Arguments.of(1.4490282054459427, 1.048502805426968, 8, 8),
            Arguments.of(2.929848286313245, 1.1163961491886014, 9, 7),
            Arguments.of(2.1123819632831826, 1.2045804557520403, 9, 7)
        );
    }

    /**
     * Test the cancellation for negative half-integer a using {@code a = 0.5 - 2^b}.
     */
    @ParameterizedTest
    @CsvSource({
        "3, 9, 30",
        "21, 21, 30",
    })
    @Disabled("Used to show cancellation at half-integer a")
    void testCancellation(int ls, int us, int maxB) {
        Assertions.assertTrue(ls >= 3);
        Assertions.assertTrue(us >= ls);
        Assertions.assertTrue(maxB >= 0);
        int maxCancellation = 0;
        for (int s = ls; s <= us; s += 2) {
            for (int b = 0; b <= maxB; b++) {
                final double a = 0.5 - Math.scalb(1.0, b);
                // Cancellation: max(exponent(a), exponent(b)) - exponent(a - b)
                // We know the largest term is the zeta evaluation at 0.5
                final double z = HurwitzZetaTest.zeta(s, 0.5, Context.DOUBLE);
                final double r = HurwitzZetaTest.zeta(s, 1 - a, Context.DOUBLE);
                final int cancellation = Math.getExponent(z) - Math.getExponent(r);
                maxCancellation = Math.max(maxCancellation, cancellation);
                report("|%d|%s|%.4g|%.4g|%d|%n", s, a, z, r, cancellation);
            }
        }
        Assertions.assertTrue(maxCancellation > 53);
    }

    /**
     * Create test data by finding roots of the the zeta function.
     *
     * <p>The approximate cancellation of the positive and negative terms is computed in
     * bits. This does not exceed 55 bits, and trends down as s and |a| increase due to
     * the increasing magnitude of ulp(a) affecting the largest terms.
     *
     * <p>This sets the limit on double-double computation as 2 ulp in the double result.
     * This would require exact 106-bit double-double arguments, which is not possible.
     * The DD zeta function can be optimised to ~105 bits of precision.
     *
     * <p>This uses positive a parameters for convenience. Cases with cancellation below a
     * threshold are ignored.
     */
    @ParameterizedTest
    @CsvSource({
        // Notes on cancellation in the computation of sum_(n=0)^infty (a+n)^-s
        // when a < 0 and s is odd.
        //
        // Sides of the computation (i.e. positive and negative sums)
        //   x = a - ceil(a) : x in -(1, 0)
        //   zeta(s, x + 1) - [ zeta(s, -x) - zeta(s, 1 - a) ]
        //
        // The most significant terms on each side are: -1 < x < 0 < 1 + x < 1.
        // The positive side is larger than the negative side when x = -0.5
        // as the series is infinite on one side and not the other.
        // The root must be x >= -0.5 to bias the most significant term to be larger
        // on the negative side.
        //
        // As s increases the magnitude of each side increases as the terms with (0, 1)^-s
        // increase in size. For the two sides to cancel x -> -0.5.
        // When s is large the two terms (0.5 +/- ulp(a))^-s dominate the two sides.
        // If s >= 1067 the difference between terms overflows for any a.
        // The true root is never at a half-integer as total cancellation will result
        // in the positive value zeta(s, 1 - a). However in double precision
        // the value of x is limited by ulp(a). So the smallest value may be computed
        // when x is half-integer. Thus zeta(s, 1 - a) is the upper bound on the value at
        // the root.
        //
        // If |a| is large the cancellation is less around (0.5 +/- ulp(a))^-s as the ulp
        // of |a| pushes the dominant terms ~0.5^-s away from each other. However the
        // sides have more terms and cancellation can still be large.
        //
        // When a single negative term (-0.5+ulp(a))^-s is larger than the entire positive
        // side (zeta(s, 0.5+ulp(a)) then the roots must be half-integer.
        // Detection of half-integer roots:
        // s  2^b
        // 31 0
        // 29 1
        // 27 4
        // 25 7
        // 23 11
        // 21 14   <-- Impractical to search after this point
        // 19 17
        // 17 20
        // 15 24
        // 13 27
        // 11 31
        // 9 34
        // 7 38
        // 5 41
        // 3 46
        //
        // The following cases should find roots where cancellation is high and record
        // them to be verified using e.g. mpmath or MATLAB zeta functions.
        // Note: Inspection of the output from an independent zeta evaluation should
        // see a sign change for each case as they are roots in double precision.
        // The zeta values computed by this implementation are recorded in the output
        // so visual verification is possible.
        //
        // When s is large, search small |a|.
        // When s is small, search larger |a|.
        // When |a| is very large, all roots are at half-integer

        // Max cancellation 55-bits (~1.9s)
        "3, 21, 0, 6",
        // Max cancellation 52-bits (~0.5s)
        "23, 1067, 0, 6",
        // Max cancellation 52-bits (~36s)
        "3, 25, 7, 11",
        // Max cancellation 51-bits (~250s)
        "3, 21, 12, 14",
        // Higher b is impractical to search.
        // This requires optimisation of the BigDecimal implementation to detect
        // when power terms are sub-normal and stop.
        // "3, 19, 15, 17",
    })
    @Disabled("Used to generate test data")
    void testDataZetaRoots(int ls, int us, int minb, int maxb) throws IOException {
        // Validate arguments
        Assertions.assertTrue(ls > 2, "s must be >= 3");
        Assertions.assertTrue(minb >= 0, "2^b must be >= 1");
        Assertions.assertTrue(maxb <= 51, "max a is too large as ulp(2^b - 0.5) >= 0.5");
        Assertions.assertEquals(1, ls & 1, "s must be odd");
        // Context for final evaluation around the root
        // quad-double precision should be able to evaluate to a double-double result
        // Changing this to DOUBLE_DOUBLE speeds up the test by a factor of 0.6.
        Context context;
        context = Context.BD_QUAD_DOUBLE;
        // context = Context.BD_DOUBLE_DOUBLE;
        // Maximum cancellation
        double maxc = 0;
        boolean nonHalfIntegerRoot = false;
        // Cases to record
        final ArrayList<String> cases = new ArrayList<>();
        // Threshold to include the case in the results
        final double threshold = 45;
        // Used to store information about the maximum a for all s
        final double[] roots = new double[us + 2];
        final double[] cancellation = new double[us + 2];
        // Lowest tolerance allowed
        final BrentSolver solver = new BrentSolver(0, 0, 0);
        // Search for roots in (2^-b, 0) for all b.
        // Each increment in b increases the ulp which eventually forces all roots to be half-integer.
        int upperS = us;
        long time = System.nanoTime();
        for (int b = minb; b <= maxb; b++) {
            if (upperS < ls) {
                // Done
                break;
            }
            final double ua = Math.scalb(1.0, b);
            final double la = ua > 1 ? ua * 0.5 : 0;
            final double ulp = Math.ulp(ua - 0.5);
            // Iterate s
            for (int s = ls; s <= upperS; s += 2) {
                boolean halfIntegerRoots = roots[s] > 0;
                if (!halfIntegerRoots) {
                    // If -0.5+ulp is so small that s creates a term larger than the entire
                    // positive side (zeta(s, 0.5+ulp)) then eliminate larger s and a.
                    final BigDecimal t1 = zeta(s, new BigDecimal(0.5 + ulp), null, context);
                    final BigDecimal t2 = new BigDecimal(0.5 - ulp).pow(-s, context.getMathContext());
                    if (t2.compareTo(t1) >= 0) {
                        // Any a at this ulp or larger is eliminated for this s.
                        // Set a flag for half integer roots to skip the root finder.
                        // The root will still be tested using the cancellation and included
                        // if above threshold.
                        roots[s] = b + 1;
                        halfIntegerRoots = true;
                    }
                }
                final int ss = s;
                // Assume the function is optimised for accuracy
                final DoubleUnaryOperator f = x -> HurwitzZetaTest.zetaNegativeBD(ss, x);
                // Check if any cases were recorded
                final int casesCount = cases.size();
                for (double ta = la; ta < ua; ta += 1) {
                    // test a is integer: bracket -(a, a+1)
                    // The root should always have a - Math.ceil(a) >= -0.5.
                    // The half-integer point is a good first approximation.
                    // Setting mid as the next value will short circuit the root finding
                    // when a is large and the root is at half-integer.
                    double xx;
                    final double min = -ta - 0.5;
                    if (halfIntegerRoots) {
                        // The assertion below will verify this is the root
                        xx = min;
                    } else {
                        final double mid = Math.nextUp(min);
                        final double max = Math.nextDown(-ta);
                        xx = solver.findRoot(f, min, mid, max);
                    }
                    // Check the solver found a bracket
                    final double x0 = Math.nextDown(xx);
                    final double x1 = Math.nextUp(xx);
                    final double f0 = f.applyAsDouble(x0);
                    final double fx = f.applyAsDouble(xx);
                    final double f1 = f.applyAsDouble(x1);
                    // Root should be bracketed by a sign change, and be the minimum
                    Assertions.assertTrue(f0 > 0 && f1 < 0,
                        () -> String.format("No sign change: %d %s %s %s %s%n", ss, xx, f0, fx, f1));
                    Assertions.assertTrue(Math.abs(fx) < Math.min(f0, -f1),
                        () -> String.format("Not the minimum: %d %s %s %s %s%n", ss, xx, f0, fx, f1));
                    nonHalfIntegerRoot |= xx - Math.ceil(xx) != -0.5;
                    // Compute the cancellation using sides of the computation:
                    // x = a - ceil(a) : x in -(1, 0)
                    // zeta(s, x + 1) +/- [ zeta(s, -x) - zeta(s, 1 - a) ]
                    // Cancellation is the power of 2 magnitude difference.
                    final double[] args = {x0, xx, x1};
                    final double[] results = {f0, fx, f1};
                    // Record the root (if not half-integer) and both sides.
                    // Only do this when cancellation of is above the threshold for 1 of the results.
                    final ArrayList<String> record = new ArrayList<>();
                    double casec = 0;
                    for (int i = 0; i < 3; i++) {
                        final double a = args[i];
                        final double x = a - Math.ceil(a);
                        if (x == -0.5 || !Double.isFinite(results[i])) {
                            // Skip computing cancellation.
                            // Include the case so it can be evaluated for test resource data.
                            record.add(String.format("%s, %s%n", s, a));
                            continue;
                        }
                        // Get the terms that cancel
                        final BigDecimal z1 = zeta(s, BigDecimal.ONE.add(new BigDecimal(x)), null, context);
                        final BigDecimal z2 = negativeSeriesSum(a, x, s, context);
                        // Verify the terms are correct
                        TestUtils.assertEquals(results[i],
                            z1.add(z2, context.getMathContext()).doubleValue(), 0, null,
                            () -> String.format("%d %s %s", ss, a, z1.doubleValue()));
                        // Cancellation is the number of matching leading bits:
                        // r = x - y
                        // max(exponent(x), exponent(y)) - exponent(r)
                        final BigDecimal zz = z1.compareTo(z2.abs()) > 0 ? z1 : z2.abs();
                        final double z = zz.doubleValue();
                        double lz;
                        if (Double.isFinite(z)) {
                            lz = Math.getExponent(z);
                        } else {
                            // floor(log2(max(|x|, |y|))) - floor(log2(r)) ~ log2(max(|x|, |y|) / r)
                            // log2(z) == log10(z) / log10(2)
                            // precision - scale = floor(log10(z))
                            // The floor operation is before conversion to base 2 so is approximate
                            lz = (zz.precision() - zz.scale()) / Math.log10(2);
                        }
                        // If result is 0 the exponent is -1023. The cancellation is total and
                        // computed as the number of binary digits in z with trailing zeros.
                        final double lr = Math.getExponent(results[i]);
                        final double cx = lz - lr;
                        casec = Math.max(casec, cx);
                        maxc = Math.max(maxc, cx);
                        // Record the case
                        record.add(String.format("# %s : %s  %s%n%s, %s%n",
                            zz.round(new MathContext(4)).toEngineeringString(), shortFormat(cx),
                            results[i], s, a));
                    }
                    // Only include if at least one is above threshold
                    if (casec >= threshold) {
                        cases.addAll(record);
                        // Store largest a for this s
                        roots[s + 1] = xx;
                        // Record max cancellation for s
                        if (casec >= cancellation[s]) {
                            cancellation[s] = casec;
                            cancellation[s + 1] = xx;
                        }
                    }
                }
                if (casesCount == cases.size() && halfIntegerRoots) {
                    // No cases recorded for this s.
                    // If half-integer roots have been detected, no larger a
                    // will have non-half-integer roots since cancellation of the 0.5+/-ulp
                    // sides will reduce at higher ulp.
                    upperS = s - 2;
                    // For this ulp no higher s will have more cancellation
                    break;
                }
            }
        }
        time = System.nanoTime() - time;
        final String msg = String.format("max cancellation %s; non-half-integer root=%s", shortFormat(maxc), nonHalfIntegerRoot);
        Assertions.assertFalse(cases.isEmpty(), () -> "No test cases recorded: " + msg);
        Assertions.assertTrue(nonHalfIntegerRoot);
        Assertions.assertTrue(maxc <= 55, "Exceeded 55 bits: " + msg);

        try (PrintStream out = getPrintStream(
            String.format("hzeta_root_s%d_%d_na2b%d_2b%d.txt", ls, us, minb, maxb))) {
            final String la = minb > 0 ? "2^" + (minb - 1) : "0";
            out.printf("# s in [%d, %d] : a in [%s, 2^%d]%n", ls, us, la, maxb);
            out.printf("# Runtime: %.3fs%n", time * 1e-9);
            out.printf("# Cancellation of terms (x - y) computed using:%n");
            out.printf("# max(exponent(x), exponent(y)) - exponent(x - y)%n");
            out.printf("# Cancellation threshold = %s%n", shortFormat(threshold));
            out.printf("# Maximum cancellation (a - ceil(a) != -0.5) = %s%n", shortFormat(maxc));
            out.printf("# N = %d%n", cases.size());
            // Record detection of half-integer-roots
            boolean header = false;
            for (int s = 3; s < roots.length; s += 2) {
                if (roots[s] > 0) {
                    if (!header) {
                        header = true;
                        out.printf("# Half-integer roots:%n");
                        out.printf("# s max(|a|)%n");
                    }
                    out.printf("# %d 2^%d%n", s, (int) roots[s] - 1);
                    if (roots[s] == 1) {
                        out.printf("# [%d, %d] 2^0%n", s + 2, us);
                        break;
                    }
                }
            }
            // Summarise cases min a for all s
            out.printf("# Cases above threshold:%n");
            out.printf("# s min(a) : max_cancellation a%n");
            for (int s = 3; s < roots.length; s += 2) {
                if (roots[s + 1] < 0) {
                    out.printf("# %d %s : %s %s%n", s, roots[s + 1],
                        shortFormat(cancellation[s]), cancellation[s + 1]);
                }
            }
            out.printf("# Comment shows max(|x|, |y|), cancellation bits and computed value of zeta%n");
            cases.forEach(out::print);
        }
    }

    /**
     * Create test data for integer {@code s} and {@code a in [la, ua)}.
     * Samples can follow a log-uniform limiting distribution with full randomisation
     * of the 52-bit mantissa.
     * Used to generate data for the high-precision zeta functions.
     */
    @ParameterizedTest
    @CsvSource({
        // zeta called with 0 < a < 1. Use a close to 1:
        // 0.9999999999988898 =  1.0 - 10000 * 0x1p-53
        "2, 11, 0.9999999999988898, 1, false",
        // zeta called with a > 1.
        // Occurs when a > 2N with N the number of power terms in a zeta evaluation.
        "2, 11, 40, 41, true",
    })
    @Disabled("Used to generate test data")
    void testDataSampleIntegerS(int ls, int us,
                                double la, double ua, boolean uniforma) throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Validate arguments
        Assertions.assertTrue(ls > 1);
        Assertions.assertTrue(us >= ls);
        Assertions.assertTrue(la >= 0);
        Assertions.assertTrue(ua > la);
        // Create samplers
        final int range = us - ls + 1;
        final DoubleSupplier s = () -> ls + rng.nextInt(range);
        final DoubleSupplier a = createSampler(rng, la, ua, uniforma);

        final int size = 3000;
        try (PrintStream out = getPrintStream(
            String.format("hzeta_is%s_%s_a%s_%s.txt",
                shortFormat(ls), shortFormat(us), shortFormat(la), shortFormat(ua)))) {
            for (int i = 0; i < size; i++) {
                out.printf("%s, %s%n", s.getAsDouble(), a.getAsDouble());
            }
        }
    }

    /**
     * Create test data for {@code s in [ls, us)} and {@code a in [la, ua)}.
     * Samples can follow a log-uniform limiting distribution with full randomisation
     * of the 52-bit mantissa.
     */
    @ParameterizedTest
    @CsvSource({
        "1.0000000000000002, 4, false, 1, 8, false",
        "1.0000000000000002, 4, false, 8, 32, false",
        "1.0000000000000002, 4, false, 32, 2.147483648E9, true",
        // Up to the asymptotic approximation.
        // (2^-53)^-19 is ~ 7.3e-304 so cannot generate sub-normal terms.
        // A separate method creates data with sub-normal terms and results.
        "1.0000000000000002, 19, true, 2.147483648E9, 1e16, true",
        "1.0000000000000002, 4, false, 0, 1, true",
        "1.0000000000000002, 4, false, 1e-16, 1e-14, false",
        "4, 32, true, 1, 8, false",
    })
    @Disabled("Used to generate test data")
    void testDataSample(double ls, double us, boolean uniforms,
                        double la, double ua, boolean uniforma) throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Validate arguments
        Assertions.assertTrue(ls > 1);
        Assertions.assertTrue(us > ls);
        Assertions.assertTrue(la >= 0);
        Assertions.assertTrue(ua > la);
        // Create samplers
        final DoubleSupplier s = createSampler(rng, ls, us, uniforms);
        final DoubleSupplier a = createSampler(rng, la, ua, uniforma);

        final int size = 3000;
        try (PrintStream out = getPrintStream(
            String.format("hzeta_s%s_%s_a%s_%s.txt",
                shortFormat(ls), shortFormat(us), shortFormat(la), shortFormat(ua)))) {
            for (int i = 0; i < size; i++) {
                out.printf("%s, %s%n", s.getAsDouble(), a.getAsDouble());
            }
        }
    }

    /**
     * Create test data for {@code s > a + N} where the tail series may eventually diverge
     * as the Pochammer term (s)_2k-1 (eventually) grows faster than the power term reduces.
     *
     * @see #testTailConvergenceDD(int, boolean, int, int, int, int)
     */
    @ParameterizedTest
    @CsvSource({
        // Set using 2^-s where the full precision of the result of any term (a+k)^-s is below
        // the min double value:
        // double 1075
        // double-double 1075+53 = 1128
        "1075, 9",
        "1128, 15",
        // Maximum s used in the total cancellation case (only requires double precision)
        "1837, 9",
    })
    @Disabled("Used to generate test data")
    void testDataTailConvergence(int maxS, int n) throws IOException {
        // Validate arguments
        Assertions.assertTrue(maxS > n);
        Assertions.assertTrue(n > 1);
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Create samplers
        final IntSupplier s = () -> n + rng.nextInt(maxS - n);

        final int size = 3000;
        // Limit a to the range where a^-s will not be significantly larger than (a+n)^-s
        final double minA = 1.5;
        try (PrintStream out = getPrintStream(
            String.format("hzeta_s%d_%d_tail.txt", n, maxS))) {
            for (int i = 0; i < size; i++) {
                final double ss = s.getAsInt();
                final double a = rng.nextDouble(Math.max(minA, ss - 2 * n), ss);
                out.printf("%s, %s%n", ss, a);
            }
        }
    }

    /**
     * Create data within a line for the zeta domain for (s, a) == x.
     * The line is loaded from a CSV file using records: s, a, z.
     * The file should contain reference data from an independent zeta implementation.
     * The expected value must be provided as a power of 2.
     *
     * <p>Any points not close to the line are discard. This can occur
     * when 2^b approaches zero. Trapezoids are created under each successive
     * pair of points (a_i, a_(i+1)). The size of the trapezoid is computed
     * and then sampled are generated across the entire domain under the line
     * by sampling a trapezoid proportional to its area, then sampling randomly
     * with the trapezoid by rejection sampling.
     *
     * <p>Note that data can only be generated within the bounds defined
     * by the input file. Points on the curve can be discarded using the parameters
     * [minS, maxS] and [minA, maxA] to bracket the range. Ranges are ignored
     * if the lower and upper bounds match.
     *
     * @see #testDataSubNormalDomain(int)
     */
    @ParameterizedTest
    @CsvSource({
        // "hzeta_2p-1022.csv, -1022, 1e-6, 50000, 0, 0, 0, 0",
        "hzeta_2p-1075.csv, -1075, 1e-6, 50000, 0, 0, 0, 0",
    })
    @Disabled("Used to generate test data")
    void testDataWithinDomain(String filename, int b, double eps, int samples,
        double minS, double maxS, double minA, double maxA) throws IOException {
        // Compute the filtering of values using normal numbers
        int offset = Double.MIN_EXPONENT + 53 - b;
        final boolean scaled = offset > 0;
        double t;
        BigDecimal scale = null;
        if (scaled) {
            // Scale back to within a normal number
            t = Math.scalb(1.0, b + offset);
            scale = new BigDecimal(BigInteger.ONE.shiftLeft(offset));
        } else {
            t = Math.scalb(1.0, b);
        }
        Assertions.assertTrue(t > 0, "Outside the scaling range of a double");
        Assertions.assertTrue(eps >= 0);
        final DoublePredicate filterS = minS < maxS ? s -> s < minS || s > maxS : s -> false;
        final DoublePredicate filterA = minA < maxA ? a -> a < minA || a > maxA : a -> false;
        final ArrayList<double[]> points = new ArrayList<>();
        final double delta = t * eps;
        int removed = 0;
        try (DataReader in = new DataReader(filename)) {
            while (in.next()) {
                try {
                    double z;
                    if (scaled) {
                        z = in.getBigDecimal(2).multiply(scale, MathContext.DECIMAL64).doubleValue();
                    } else {
                        z = in.getDouble(2);
                    }
                    // Ignore z not close to the target
                    if (Math.abs(z - t) > delta) {
                        removed++;
                        continue;
                    }
                    final double s = in.getDouble(0);
                    final double a = in.getDouble(1);
                    // Ignore outside the configured ranges
                    if (filterS.test(s) || filterA.test(a)) {
                        removed++;
                        continue;
                    }
                    points.add(new double[] {s, a, z});
                } catch (final NumberFormatException ex) {
                    Assertions.fail("Failed to load data: " + Arrays.toString(in.getFields()), ex);
                }
            }
        } catch (final IOException ex) {
            Assertions.fail("Failed to load data: " + filename, ex);
        }

        // If any a or s are the same then use the record closest to t
        final Comparator<double[]> cs = Comparator.comparingDouble(k -> k[0]);
        final Comparator<double[]> ca = Comparator.comparingDouble(k -> k[1]);
        // Sort by s then a
        points.sort(cs.thenComparing(ca));
        int current = 0;
        for (int i = 1; i < points.size(); i++) {
            if (points.get(i)[1] == points.get(current)[1]) {
                // Keep the result closest to t
                if (Math.abs(points.get(i)[2] - t) < Math.abs(points.get(current)[2] - t)) {
                    points.set(current, null);
                    current = i;
                } else {
                    points.set(i, null);
                }
                removed++;
            } else {
                current = i;
            }
        }
        points.removeIf(x -> x == null);

        // Sort by a then s
        points.sort(ca.thenComparing(cs));
        current = 0;
        for (int i = 1; i < points.size(); i++) {
            if (points.get(i)[1] == points.get(current)[1]) {
                // Keep the result closest to t
                if (Math.abs(points.get(i)[2] - t) < Math.abs(points.get(current)[2] - t)) {
                    points.set(current, null);
                    current = i;
                } else {
                    points.set(i, null);
                }
                removed++;
            } else {
                current = i;
            }
        }
        points.removeIf(x -> x == null);

        // Data is now sorted by ascending a. It should be descending for s.
        // If not then the curve is not always at 2^b and data generation will not
        // represent the space under the curve a 2^b. We could continue
        // and sample anyway. But we fail and require a better input curve.
        double maxD = 0;
        for (int i = 0; i < points.size(); i++) {
            final double d = Math.abs(points.get(i)[2] - t);
            maxD = Math.max(maxD, d);
            if (i != 0) {
                Assertions.assertTrue(points.get(i)[0] < points.get(i - 1)[0]);
            }
        }
        Assertions.assertTrue(points.size() > 1, "No points to create samples under the curve");

        // Create right trapezoids under the curve for all a-axis points.
        // The two right angles are between the base along the a-axis and the heights in s.
        // We have n-1 regions. First convert the double to raw bits
        // and subtract from the minimum to sample. This is the number of double
        // values between the min and the curve point.
        // Then compute the the width and area.
        final long[] w = new long[points.size()];
        final long[] h = new long[points.size()];
        final long bitsS = Double.doubleToRawLongBits(Math.nextUp(1.0));
        for (int i = 0; i < points.size(); i++) {
            double[] p = points.get(i);
            // The height is sampled from the minimum to the highest point of the trapezoid
            h[i] = Double.doubleToLongBits(p[0]) - bitsS;
            // The width is sampled using the difference between parallel side
            w[i] = Double.doubleToLongBits(p[1]);
            Assertions.assertTrue(h[i] > 0);
            Assertions.assertTrue(w[i] > 0);
        }
        final double[] probability = new double[w.length - 1];
        double total = 0;
        for (int i = 0; i < probability.length; i++) {
            // Area = width * height at the middle
            probability[i] = (double) (w[i + 1] - w[i]) * ((h[i + 1] + h[i]) >>> 1);
            total += probability[i];
        }
        // Sample probability = area / total
        for (int i = 0; i < probability.length; i++) {
            probability[i] /= total;
        }
        // Sample trapezoids uniformly. Then sample within
        // the chosen trapezoid uniformly. Rejection should be < 50%.
        // For consistency we keep the same random source as other tests.
        final SplittableRandom rng = new SplittableRandom(SEED);
        final SharedStateDiscreteSampler sampler =
            AliasMethodDiscreteSampler.of(rng::nextLong, probability);
        int rejects = 0;
        // This could write directly to the output file.
        // For convenience the data are sorted by a to allow slicing the
        // data to domains supported by reference implementations.
        final ArrayList<double[]> data = new ArrayList<>(samples);
        final int ilogb0 = ilogb(0);
        final int[] values = doReporting() ? new int[ilogb(Double.POSITIVE_INFINITY) - ilogb0 + 1] : null;
        for (int n = 0; n < samples; n++) {
            final int i = sampler.sample();
            final long h1 = h[i];
            final long h2 = h[i + 1];
            final long wa = w[i + 1] - w[i];
            // Rejection sampling is worse case 50% rejection if h2 ~ 0 and this is
            // a right angled triangle. h2 < h1 due to the prior check in descending s.
            //
            // h1
            // |      h2
            // |      |
            // |      |
            // +--wa--+
            long s;
            long a;
            for (;;) {
                s = rng.nextLong(h1);
                a = rng.nextLong(wa);
                if (s > h2) {
                    // Check if the height is below the line joining h1 to h1.
                    // Interpolate
                    final double x = (double) a / wa;
                    final long y = (long) ((1 - x) * h1 + x * h2);
                    if (s >= y) {
                        rejects++;
                        continue;
                    }
                }
                break;
            }
            // Map the raw bits back to double values
            final double[] x = {Double.longBitsToDouble(bitsS + s),
                                Double.longBitsToDouble(w[i] + a)};
            data.add(x);
            if (values != null) {
                // Get the value
                final double z = zeta(x[0], x[1], Context.DOUBLE_REF);
                values[ilogb(z) - ilogb0]++;
            }
        }
        if (values != null) {
            for (int i = 0; i < values.length; i++) {
                if (values[i] != 0) {
                    report("2^%d = %d%n", i + ilogb0, values[i]);
                }
            }
        }
        data.sort(ca.thenComparing(cs));
        // Write the data directly to file. Name the file using range filters and the target power.
        String rangeS = "";
        String rangeA = "";
        if (minS < maxS) {
            rangeS = String.format("_s%s_%s", shortFormat(minS), shortFormat(maxS));
        }
        if (minA < maxA) {
            rangeA = String.format("_s%s_%s", shortFormat(minA), shortFormat(maxA));
        }
        try (PrintStream out = getPrintStream(
            String.format("hzeta_domain_%s%s2p%d.txt", rangeS, rangeA, b))) {
            // Report statistics on the data
            out.printf("# Input: %s%n", filename);
            out.printf("# Target: 2^%d%n", b);
            if (minS < maxS) {
                out.printf("# Filtered: s in [%s, %s]%n", minS, maxS);
            }
            if (minA < maxA) {
                out.printf("# Filtered: a in [%s, %s]%n", minA, maxA);
            }
            out.printf("# Points: %d / %d%n", points.size(), points.size() + removed);
            out.printf("# Max relative deviation = %.5g%n", maxD / t);
            final double[] copy = probability.clone();
            Arrays.sort(copy);
            out.printf("# Region sample probability: Mean %s; Median %s%n",
                1.0 / probability.length, copy[copy.length >>> 1]);
            out.printf("# Samples: %d%n", samples);
            out.printf("# Rejection rate: %.5g%n", (double) rejects / samples);
            data.forEach(x -> out.printf("%s, %s%n", x[0], x[1]));
        }
    }

    /**
     * Get integer exponent of a floating-point value.
     *
     * <p>In contrast to {@link Math#getExponent(double)} this handles
     * sub-normal numbers by computing the number of leading zeros in the mantissa
     * and shifting the unbiased exponent. The result is that for all finite, non-zero,
     * numbers {@code x}, the magnitude of {@code scalb(x, -ilogb(x))} is
     * always in the range {@code [1, 2)}.
     *
     * <p>This method is a functional equivalent of the c function ilogb(double)
     * modified to not raise exceptions for special cases of numbers that are not
     * finite non-zero.
     *
     * <p>Special cases:
     *
     * <ul>
     * <li>If the argument is NaN or infinite, then the result is
     * {@link Double#MAX_EXPONENT} + 1 = 1024.</li>
     * <li>If the argument is zero, then the result is
     * {@link Double#MIN_EXPONENT} - 53 = -1075.</li>
     * </ul>
     *
     * @param x the value
     * @return The unbiased exponent of the value to be used for scaling
     * @see Math#getExponent(double)
     * @see Math#scalb(double, int)
     * @see <a href="https://www.cplusplus.com/reference/cmath/ilogb/">ilogb</a>
     */
    private static int ilogb(double x) {
        // Only interested in the exponent and mantissa so remove the sign bit
        final long bits = Double.doubleToRawLongBits(x) & Long.MAX_VALUE;
        // Get the unbiased exponent
        int exp = (int) (bits >>> 52) - 1023;

        // No case to distinguish nan/inf
        // Handle sub-normal numbers
        if (exp == Double.MIN_EXPONENT - 1) {
            // Special case for zero
            if (bits == 0) {
                return Double.MIN_EXPONENT - 53;
            }
            // A sub-normal number has an exponent below -1022. The amount below
            // is defined by the number of shifts of the most significant bit in
            // the mantissa that is required to get a 1 at position 53 (i.e. as
            // if it were a normal number with assumed leading bit)
            final long mantissa = bits & 0x000f_ffff_ffff_ffffL;
            exp -= Long.numberOfLeadingZeros(mantissa << 12);
        }
        return exp;
    }

    /**
     * Create a line in the zeta domain for (s, a) where the result is sub-normal.
     * For any {@code (s',a') < (s,a)} the result will be sub-normal or normal.
     * This is the region where the function should support accurate evaluation.
     *
     * <p>This uses the current test implementation to find (s, a). These are recorded
     * to file for independent evaluation.
     */
    @ParameterizedTest
    // TODO:
    // This uniformly steps through s or a. It should increase the number
    // of steps when the line is more convex.
    // Ideally each trapezoid under the curve should be the same area, and
    // the segment between the points (top of trapezoid) should not deviate
    // far from the curve below it.
    @ValueSource(ints = {-1022, -1065, -1075})
    @Disabled("Used to generate test data")
    void testDataSubNormalDomain(int b) throws IOException {
        // Configure for high accuracy with a scaled result.
        // Scaling by 2^106 allows the solver to find roots with a normal number
        final Context context = Context.DOUBLE_REF.addOptions(Context.SCALED_RESULT);
        final DoubleBinaryOperator zeta = (s, a) -> HurwitzZetaTest.zeta(s, a, context);
        final BrentSolver solver = new BrentSolver(0, 0, 0);

        // Map zeta(a, s) = t using n points
        // This can be changed to a lower value, e.g. 2^-(1022 + 23)
        final double t = Math.scalb(1.0, b + 106);
        final int n = 1000;
        final ArrayList<double[]> points = new ArrayList<>();

        // Map zeta(a, s) = t using a
        final double minA = Math.nextUp(1.0);
        final double maxA = 0x1p53;
        // zeta(2^61, 1+2^-52) = 4.377e-223
        // zeta(2^62, 1+2^-52) = 1.916e-445
        final double maxS = solver.findRoot(s -> zeta.applyAsDouble(s, minA) - t, 0x1p61, 0x1p62);
        // zeta(20, 2^53) = 3.65e-289
        // zeta(21, 2^53) = 4.04e-321
        // zeta(22, 2^53) = 4.28e-337
        final double minS = solver.findRoot(s -> zeta.applyAsDouble(s, maxA) - t, 20, 22);

        // a in (1, 2^53]
        // Expecting s in [20.2, 3.19e18]
        // Interpolate a in uniform steps to find s
        final DoubleSupplier nextA = createRange(n, minA, maxA, false);
        for (int i = 0; i < n; i++) {
            double a = nextA.getAsDouble();
            try {
                double s = solver.findRoot(x -> zeta.applyAsDouble(x, a) - t, minS, maxS);
                double z = zeta.applyAsDouble(s, a);
                points.add(new double[] {s, a, z});
            } catch (IllegalArgumentException e) {
                Assertions.fail(e.getMessage() + ": a=" + a);
            }
        }

        // Map zeta(a, s) = t using s
        final DoubleSupplier nextS = createRange(n, minS, maxS, false);
        // Skip the ends as the round-trip of the function is not exact.
        // This then enumerates (minS, maxS)
        nextS.getAsDouble();
        for (int i = 2; i < n; i++) {
            double s = nextS.getAsDouble();
            try {
                double a = solver.findRoot(x -> zeta.applyAsDouble(s, x) - t, minA, maxA);
                double z = zeta.applyAsDouble(s, a);
                points.add(new double[] {s, a, z});
            } catch (IllegalArgumentException e) {
                Assertions.fail(e.getMessage() + ": s=" + s);
            }
        }

        // Note:
        // Mapping s to a can create duplicates for a as a -> 1 and may create zeros.
        // This occurs when zeta(s, a) = 0 is closer to t than zeta(s, nextDown(a)):
        // zeta(2.02294791067527501E18, 1.0000000000000004) = 0
        // zeta(2.02294791067527501E18, 1.0000000000000002) = 8.35e-196

        // Sort by a (descending) then s and check the zeta result is descending for increasing s
        final Comparator<double[]> cs = Comparator.comparingDouble(k -> k[0]);
        final Comparator<double[]> ca = Comparator.comparingDouble(k -> -k[1]);
        points.sort(ca.thenComparing(cs));
        double[] previous = points.get(0);
        for (int i = 1; i < points.size(); i++) {
            final double[] current = points.get(i);
            if (current[1] == previous[1]) {
                Assertions.assertTrue(current[2] <= previous[2]);
            }
            previous = current;
        }
        // Sort by s then a (descending) and check the zeta result is ascending for decreasing a
        points.sort(cs.thenComparing(ca));
        previous = points.get(0);
        for (int i = 1; i < points.size(); i++) {
            final double[] current = points.get(i);
            if (current[0] == previous[0]) {
                Assertions.assertTrue(current[2] >= previous[2]);
            }
            previous = current;
        }

        try (PrintStream out = getPrintStream(String.format("hzeta_2p%d.txt", b))) {
            points.forEach(k -> out.printf("%s, %s%n", k[0], k[1]));
        }
    }

    /**
     * Create a line in the zeta domain for (s, a) where the result is sub-normal.
     * For any {@code (s',a') < (s,a)} the result will be sub-normal or normal.
     * This is the region where the function should support accurate evaluation.
     *
     * <p>This uses the current test implementation to find (s, a). These are recorded
     * to file for independent evaluation.
     */
    @ParameterizedTest
    @CsvSource({
        // Reflection identity must use max s=1837.
        // zeta(1837, 1.5) = 3.31e-324
        "1837, -1075",
    })
    @Disabled("Used to generate test data")
    void testDataSubNormalDomainDD(int maxS, int b) throws IOException {
        // Configure for high accuracy
        final Context context = Context.of(40, 106);
        // This may target numbers that are outside the domain of a double.
        // Use a log representation of the double-double result.
        // Requires mapping log2 exponents to base e logarithm.
        // log(x) = log2(x) * log(2)
        final double ln2 = Math.log(2);
        final int[] e = {0};
        final BiFunction<Integer, Double, DD> zeta = (s, a) -> {
            DD r = HurwitzZetaTest.zeta(s, DD.of(a), null, e, context);
            // Map the double-double result to a log result
            // log(x + xx) = log(x * (1+xx/x))
            //             = log(x) + log1p(xx/x)
            // Note: The round-off may be useless if log(x) is not exact
            // but will be monotonic which allows root finding
            return DD.ofSum(Math.log(r.hi()), Math.log1p(r.lo() / r.hi())).add(e[0] * ln2);
        };
        final BrentSolver solver = new BrentSolver(0, 0, 0);

        // Map zeta(a, s) = t using n points
        // This can be changed to a lower value
        // log(2^b)
        final double t = b * ln2;
        final ArrayList<double[]> points = new ArrayList<>();

        // Map zeta(a, s) = t using s
        final double minA = Math.nextUp(1.0);
        final double maxA = 0x1p53;
        for (int i = 2; i <= maxS; i++) {
            final int s = i;
            double a;
            // Note: This starts at s=2. The first s that achieves a sub-normal result
            // with max a is (2^53)^-s = 2^-1022 => s = -log2(2^1022) / log2(2^53) = 19.28.
            // The target may be lower. Catch the no bracket exception and use maxA.
            try {
                a = solver.findRoot(x -> zeta.apply(s, x).subtract(t).hi(), minA, maxA);
                // Debug the result
                final DD z = HurwitzZetaTest.zeta(s, DD.of(a), null, e, context);
                // Check z ~ 1.0 * 2^b or 0.5 * 2^(b+1)
                Assertions.assertTrue(b <= e[0] && e[0] <= b + 1,
                    () -> "Failed to find target: " + toBigDecimal(z, e[0]));
                final double expected = e[0] == b ? 1.0 : 0.5;
                Assertions.assertEquals(expected, z.hi(), expected * 1e-12,
                    () -> "Failed to find target: " + toBigDecimal(z, e[0]));
            } catch (IllegalArgumentException ignored) {
                a = maxA;
            }
            points.add(new double[] {s, a});
        }

        try (PrintStream out = getPrintStream(String.format("hzeta_is_2p%d.txt", b))) {
            points.forEach(k -> out.printf("%s, %s%n", k[0], k[1]));
        }
    }

    /**
     * Create test data for {@code s in [ls, us)} where the terms of the series are
     * sub-normal.
     *
     * <p>Note: The 1023 and 1040 files when evaluated with mpmath have disagreements
     * with the Java implementation. These were evaluated in MATLAB. 2 cases created
     * out-of-memory errors in MATLAB and were verified using WolframAlpha online.
     * See the test resource CSV file for details.
     */
    @ParameterizedTest
    @CsvSource({
        // (x+n)^-s ~ 2^-b
        // Asymptotic approximation @ 2^53.
        // Use min s ~ b / 53 ~ 19.28
        // Only some sub-normal terms
        "20, 1024, true, 9, 1023",
        // Many terms of the summation are sub-normal
        "20, 1024, true, 1, 1040",
        "20, 1024, true, 1, 1053",
        // Result should be sub-normal
        "20, 1024, true, 0.1, 1073",
    })
    @Disabled("Used to generate test data")
    void testDataSubNormalSample(double ls, double us, boolean uniforms, double n, int b) throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Validate arguments
        Assertions.assertTrue(ls > 1);
        Assertions.assertTrue(us > ls);
        // Create samplers
        final DoubleSupplier s = createSampler(rng, ls, us, uniforms);

        final int size = 3000;
        final double factor = Math.log(2) * b;
        // Limit a to the range where a^-s will not be significantly larger than (a+n)^-s
        final double minA = 1.5;
        try (PrintStream out = getPrintStream(
            String.format("hzeta_s%s_%s_ap%s_subnormal%d.txt",
                shortFormat(ls), shortFormat(us), shortFormat(n), b))) {
            // Allow for rejection
            for (int i = 0, j = 0; i < 10 * size && j < size; i++) {
                final double ss = s.getAsDouble();
                // min normal = 2^1022. Use 1023 to ensure sub-normal (a+n).
                // (x+n)^-s ~ 2^-1023
                // x+n = exp(log(2) * 1023 / s)
                final double apn = Math.exp(factor / ss);
                // Require some range for selection of a
                if (apn - 1 < minA) {
                    // Reject ss as too large
                    continue;
                }
                j++;
                // If (a+n) is small, sample [a, a+n); otherwise use directly
                final double a = apn - n < apn ? rng.nextDouble(Math.max(minA, apn - n), apn) : apn;
                out.printf("%s, %s%n", ss, a);
            }
        }
    }

    /**
     * Create test data for {@code s in [ls, ls + 2, ..., us - 2, us]}
     * and {@code a in -( [la, ua] + offset +/- 2^-b )}.
     *
     * <p>The parameters for {@code s} create either an even or odd series.
     *
     * <p>The parameter {@code a} is a uniform integer sample in a range. This has an
     * offset applied and a uniform wobble. This allows sampling around the poles
     * at integer values, or the location of greatest cancellation when {@code s}
     * is odd at half-integer values.
     *
     * <p>Note: Evaluation of large negative a is not possible using the resource
     * script {@code hzeta.py}. Use of mpmath returns complex results when abs(a) >> dps
     * (dps is the mpmath digits of decimal precision). Using dps=70 is sufficient.
     */
    @ParameterizedTest
    @CsvSource({
        // a = -[1.5, 15.5] +/- 0.5
        "2, 4, 1, 15, 0.5, 1",
        "3, 5, 1, 15, 0.5, 1",
        // a = -[40.5, 99.5] +/- 0.5
        "2, 4, 40, 99, 0.5, 1",
        "3, 5, 40, 99, 0.5, 1",
        // a = -[1, 15] +/- 9.31e-10
        // This approaches the pole at a = -1, -2, -3, ... and is easy to compute as
        // a single term dominates the result
        "2, 4, 1, 15, 0.0, 30",
        "3, 5, 1, 15, 0.0, 30",
        // a = -[1.5, 15.5] +/- 9.31e-10
        // This creates large cancellation when a ~ half-integer and requires
        // an extended precision power function to maintain precision over all
        // terms that cancel.
        "2, 4, 1, 15, 0.5, 30",
        "3, 5, 1, 15, 0.5, 30",
        // a = -[40.5, 99.5] +/- 9.31e-10
        "3, 5, 40, 99, 0.5, 30",
        // larger s and a. Avoid the poles at integer a.
        "6, 50, 400, 999, 0.5, 3",
        "7, 51, 400, 999, 0.5, 3",
    })
    @Disabled("Used to generate test data")
    void testDataNegativeASample(int ls, int us,
                                 int la, int ua, double offset, int b) throws IOException {
        final SplittableRandom rng = new SplittableRandom(SEED);
        // Validate arguments
        Assertions.assertTrue(ls >= 2);
        Assertions.assertTrue(us >= ls);
        Assertions.assertTrue((ls & 1) == (us & 1), "Must be odd or even series");
        Assertions.assertTrue(la >= 0);
        Assertions.assertTrue(ua > la);
        Assertions.assertTrue(b > 0);
        final double scale = Math.scalb(1, -b);
        // Check randomness can be added to the smallest a value
        Assertions.assertTrue(la + offset + 2 * scale != la + offset);
        // Create samplers
        // s in [ls, us]
        final int n = 1 + (us - ls) / 2;
        final IntSupplier s = () -> ls + 2 * rng.nextInt(n);
        // Create a signed random value in [-1, 1) and scale down.
        // signed 53-bits * 2^-53 * scale
        final double f = 0x1.0p-53 * scale;
        final int ra = ua - la + 1;
        final DoubleSupplier a = () -> la + rng.nextInt(ra) + offset + f * (rng.nextLong() >> 10);

        final int size = 3000;
        final StringBuilder name = new StringBuilder("hzeta_s")
            .append(ls).append('_').append(us).append("_na")
            .append(la).append('_').append(ua);
        if (offset != 0) {
            name.append("_p").append(shortFormat(offset));
        }
        name.append("_p0x1p").append(-b).append(".txt");
        try (PrintStream out = getPrintStream(name.toString())) {
            for (int i = 0; i < size; i++) {
                // skip the unlikely integer values of a
                final double x = a.getAsDouble();
                if (Math.rint(x) != x) {
                    out.printf("%d, %s%n", s.getAsInt(), -x);
                }
            }
        }
    }

    /**
     * Creates the sampler. The sample can be log-uniform in the range or uniform.
     *
     * @param rng the source of randomness
     * @param a the lower bound (inclusive)
     * @param b the upper bound (exclusive)
     * @param uniform if true sample from a uniform distribution.
     * @return the double supplier
     */
    private static DoubleSupplier createSampler(SplittableRandom rng, double a, double b, boolean uniform) {
        if (uniform) {
            final double range = b - a;
            return () -> a + rng.nextDouble() * range;
        }
        // limiting log-uniform distribution
        final long bits = Double.doubleToLongBits(a);
        final long range = Double.doubleToLongBits(b) - bits;
        return () -> Double.longBitsToDouble(bits + rng.nextLong(range));
    }

    /**
     * Enumerates the range. The enumeration can be log-uniform in the range or uniform.
     *
     * <p>When log uniform the upper bound may not be reached as the step increment
     * is the range in bits divided by (n-1). If not an exact factor then the upper
     * limit after {@code n} may not be b. The points in bits for i in [0, n-1]:
     * {@code a <= a + i * (range / (n - 1)) <= b}.
     *
     * @param a the lower bound (inclusive)
     * @param b the upper bound (exclusive)
     * @param uniform if true enumerate from a uniform distribution.
     * @return the double supplier
     */
    private static DoubleSupplier createRange(int n, double a, double b, boolean uniform) {
        if (uniform) {
            final int[] i = {-1};
            final double range = b - a;
            final double step = n == 0 ? range : range / (Math.abs((double) n) - 1);
            return () -> {
                i[0]++;
                return a + i[0] * step;
            };
        }
        // limiting log-uniform distribution
        final long bits = Double.doubleToLongBits(a);
        final long range = Double.doubleToLongBits(b) - bits;
        final long step = n == 0 ? range : range / (Math.abs((long) n) - 1);
        final long i[] = {bits - step};
        return () -> {
            i[0] += step;
            return Double.longBitsToDouble(i[0]);
        };
    }

    /**
     * Format the double to a short string. Assumes doubles can be close to integer and
     * returns an integer representation.
     *
     * @param d the double
     * @return the string
     */
    private static String shortFormat(double d) {
        // Close to integer
        if (Math.abs(Math.rint(d) - d) < Math.abs(d) * 1e-5) {
            d = Math.rint(d);
        }
        if ((long) d == d) {
            return Long.toString((long) d);
        }
        final String s = Double.toString(d);
        // Remove trailing zeros for shorter length on integer representations
        return s.replaceFirst("\\.0$", "").replaceFirst("\\.0E", "E");
    }

    /**
     * Gets the prints the stream.
     * Adds a header line to the output indicating how the data was created.
     *
     * @param filename the filename
     * @return the stream
     */
    private static PrintStream getPrintStream(String filename) throws IOException {
        return getPrintStream(filename, HurwitzZetaTest.class.getSimpleName());
    }

    /**
     * Gets the prints the stream.
     * Adds a header line to the output indicating how the data was created.
     *
     * @param filename the filename
     * @param source the source of the data
     * @return the stream
     */
    static PrintStream getPrintStream(String filename, String source) throws IOException {
        final PrintStream out = new PrintStream(Files.newOutputStream(Paths.get("target", filename)));
        Stream.of(
            "# Licensed to the Apache Software Foundation (ASF) under one or more",
            "# contributor license agreements.  See the NOTICE file distributed with",
            "# this work for additional information regarding copyright ownership.",
            "# The ASF licenses this file to You under the Apache License, Version 2.0",
            "# (the \"License\"); you may not use this file except in compliance with",
            "# the License.  You may obtain a copy of the License at",
            "#",
            "#     https://www.apache.org/licenses/LICENSE-2.0",
            "#",
            "# Unless required by applicable law or agreed to in writing, software",
            "# distributed under the License is distributed on an \"AS IS\" BASIS,",
            "# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.",
            "# See the License for the specific language governing permissions and",
            "# limitations under the License.",
            "",
            "# Generated by " + source)
                .forEach(out::println);
        return out;
    }
}
