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

import java.util.SplittableRandom;
import org.apache.commons.numbers.core.DD;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.MethodOrderer;
import org.junit.jupiter.api.Order;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestMethodOrder;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/**
 * Test the extended precision {@code (x+n)^z} power function for use in the Hurwitz zeta function.
 * 
 * <p>This class contains a version of {@code powxpn} for a larger domain than the
 * implementation used in the zeta function. It handles {@code a+n < 1} which can
 * overflow. The final version in the {@link HurwitzZeta} class uses a different method
 * for negative {@code a} and so the {@code powxpn} function supports {@code a + n > 1}.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class PowXPNTest {
    /** Threshold for z in the Taylor series approximation of expm1(z * t). */
    private static final double EXT_POW_LIMIT = 0x1p26;
    /** ln(0.5). */
    private static final double LN_0_5 = Math.log(0.5);
    /** ln(1.5). */
    private static final double LN_1_5 = Math.log(1.5);

    // Collect error statistics

    /** Error of Math.pow when (x+n) is exact. Should be below 1 ulp.*/
    private static TestUtils.ErrorStatistics ERROR_POW = new TestUtils.ErrorStatistics();
    /** Error of Math.pow when (x+n) is inexact. Can be very high. */
    private static TestUtils.ErrorStatistics ERROR_POWP = new TestUtils.ErrorStatistics();
    /** Error of Math.pow when (x+n) is inexact and the result was within 1 ulp. */
    private static TestUtils.ErrorStatistics ERROR_POW1 = new TestUtils.ErrorStatistics();
    /** Error of extended precision power function. */
    private static TestUtils.ErrorStatistics ERROR_POW_PLUS = new TestUtils.ErrorStatistics();
    /** Error of extended precision power function when Math.pow was within 1 ulp. */
    private static TestUtils.ErrorStatistics ERROR_POW_PLUS1 = new TestUtils.ErrorStatistics();
    /** Error of Math.pow when the final implementation uses extended precision. */
    private static TestUtils.ErrorStatistics ERROR_FINAL_STD = new TestUtils.ErrorStatistics();
    /** Error of the final implementation using extended precision. */
    private static TestUtils.ErrorStatistics ERROR_FINAL_EXT = new TestUtils.ErrorStatistics();

    /**
     * Extended precision power function {@code (x+n)^z}.
     *
     * <p>This method is used in the zeta function to compute the power function terms of
     * the infinite series {@code (a+k)^-s}. It is limited to the expected arguments
     * of the zeta function testing use cases:
     *
     * <ul>
     * <li>{@code x > -1}
     * <li>{@code n >= 1}
     * <li>{@code x+n > 0} and finite
     * <li>{@code z < -1} and finite
     * </ul>
     *
     * <p>This method has no advantage over {@code Math.pow(x+n, z)} when there is no
     * round-off component to (x+n). In the event that {@code (x+n)} is exact this method
     * will return the same result as {@code Math.pow(x+n, z)}.
     *
     * <p>Calling with {@code n=0} will consume resources for no benefit. More generally
     * this method should not be called when {@code x + n} does not create a change of
     * exponent between {@code x} and {@code x + n}, and may offer no benefit even if the
     * exponent does change depending on the arguments. For example if {@code |x| < 2^52}
     * and is an integer then addition of any {@code int} value will be exact.
     *
     * <p>Note: This method has been developed for use with integer {@code n} but can
     * be adapted to a {@code double} argument {@code n}; and for a larger domain of
     * input arguments.
     *
     * @param x Argument x ({@code 0 <= |x| < 2^53})
     * @param n Argument n ({@code > 0})
     * @param z Exponent z ({@code < -1})
     * @return the result
     */
    static double powxpn(double x, int n, double z) {
        // (s+ss)^z = s^z * (1+ss/s)^z
        //          = s^z * exp(z*log1p(ss/s))
        //          = s^z * s^z * [ (1+ss/s)^z - 1 ]
        //          = s^z + s^z * expm1(z*log1p(ss/s))
        //
        // ss/s < 1 * machine epsilon : log1p(ss/s) = ss/s
        // ss/s < 2 * machine epsilon : log1p(ss/s) ~ ss/s
        //
        // Taylor series:
        // (1 + t)^z = 1 + zt + z(z-1)t^2/2! + z(z-1)(z-2)t^3/3! + ...
        // for |t| < 1 and real z

        // Compute the exact double-double sum.
        final DD d = DD.ofSum(x, n);
        double s = d.hi();
        double ss = d.lo();

        // Skip impossible rounding
        if (ss == 0) {
            // Handles (x+n) = 0 or 1 for the correct pow result
            return Math.pow(s, z);
        }

        double t;

        // Edge cases:
        // s < 1 can overflow to infinity
        // s = 1 creates 1.0
        // s > 1 can underflow to 0
        // Care must be taken to avoid 0 * infinity == NaN

        // Prevent overflow to infinity if s < 1 && ss > 0 by
        // increasing s by 1 ulp and changing ss to negative by subtracting the ulp.

        // Prevent underflow to zero if s > 1 && ss < 0 by
        // decreasing s by 1 ulp and changing ss to positive by adding the ulp.

        // The correction must be in the same direction so we adjust (s,ss)
        // by 1 ulp to be closer to 1.
        // This does increase accuracy and may prevent infinite or zero results in extreme cases
        // but the branch statements are unpredictable.
        if (s < 1 && ss > 0) {
            t = Math.ulp(s);
            s += t;
            ss -= t;
        } else if (s > 1 && ss < 0) {
            t = Math.ulp(s);
            s -= t;
            ss += t;
        }

        double r = Math.pow(s, z);
        // Do not adjust infinity.
        // NaNs or infinite are not expected from the arguments used by the calling zeta function.
        if (Double.isInfinite(r)) {
            return r;
        }

        // 0 < t < 2^-52
        t = ss / s;

        // Limit of Taylor series where it does not equal expm1(z * t)
        // is at approximately 2^-26.
        // z * t < 2^-26  =>  z < 2^-26 / 2^-52 = 2^26
        if (z < -EXT_POW_LIMIT) {
            // z is expected to be small (practical range for a zeta result is -(1, 1075))
            // and this path is unlikely. Use either:
            //   s^z * exp(z*t)
            //   s^z + s^z * expm1(z*t)
            // Pick expm1 when its result is in the range [-0.5, 0.5].
            // Otherwise use exp (when its result is in [0, 0.5] and [1.5, infinity].
            final double zt = z * t;
            if (zt < LN_0_5 || zt > LN_1_5) {
                // Result can be 0 or infinity.
                // Avoid 0 * infinity (infinite r has been handled).
                t = Math.exp(zt);
                if (t == Double.POSITIVE_INFINITY && r == 0) {
                    // r == 0 when s < 1 and z is large negative
                    // t == infinity when ss < 0 and z is large negative
                    // assume (s,ss) < 1 and return zero
                    return 0;
                }
                // Can correct s=r=1 to infinity is ss < 0; or 0 if ss > 0
                return r * t;
            }
            // Fall-through with the (1 + t)^z - 1 result
            t = Math.expm1(zt);
        } else {
            // Compute (1 + t)^z - 1 using a Taylor series.
            // This is the most common branch path when ss != 0.
            t = z * t * (1 + (z - 1) * t * 0.5);
        }

        // t in [-0.5, 0.5] 
        // Corrects r towards overflow/underflow based on s < 1 or s > 1
        return r + r * t;
    }

    /**
     * Test the limit of the Taylor series used in the extended precision power function.
     *
     * <p>Note if the series is extended to 3 terms the threshold rises from 2^26 to 2^36.
     * The max error is still 2 ulp and the RMS stays around the same 0.7 level.
     * The power argument is not expected to be large so 2 terms is enough.
     */
    @Test
    void testExtPowThreshold() {
        final SplittableRandom rng = new SplittableRandom();
        final TestUtils.ErrorStatistics stats = new TestUtils.ErrorStatistics();
        final int maxB = Math.getExponent(EXT_POW_LIMIT);
        // Test z using powers of 2 to fully explore the range [1, LIMIT)
        for (int b = 0; b < maxB; b++) {
            double threshold = Math.scalb(1.0, b);
            for (int i = 0; i < 100; i++) {
                // The round-off in +/- [0, 2*eps)
                final double n = (0.5 - rng.nextDouble()) * 0x1p-51;
                Assertions.assertTrue(Math.abs(n) < 0x1p-52);
                final double x = rng.nextDouble() * threshold;
                // (1 + x)^n - 1 = expm1(log1p(x) * n) = expm1(x * n) when x is < 2^-52.
                final double expected = Math.expm1(x * n);
                // Taylor series: (1 + x)^n = 1 + nx + n(n-1)x^2/2! + n(n-1)(n-2)x^3/3! + ...
                // for |x| < 1 and real n
                // 2-term version
                final double actual = x * n * (1 + (x - 1) * n * 0.5);
                // 3-term version (see method javadoc why this is not required)
                //final double actual = z * t * (1 + (z - 1) * t * 0.5 * (1 + (z - 2) * t * 0.3333333333333333));
                TestUtils.assertEquals(expected, actual, 2, stats::add,
                    () -> Integer.toString(Math.getExponent(threshold)));
            }
        }
        Assertions.assertTrue(stats.getRMS() < 0.75, () -> Double.toString(stats.getRMS()));
    }

    @ParameterizedTest
    @CsvSource({
        // mpmath version 1.4.1
        // from mpmath import power, mp, mpf
        // mp.pretty = True; mp.dps = 30
        // import random
        // for b in [0, 1, 2, 3, 4]:
        //   for i in range(10):
        //     x, y, z = random.random(), random.randint(1, 9), -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "0.8967729873816368, 5, -1.6835546406097568, 0.0504229831513830078382284502922",
        "0.6813680496625345, 4, -1.4323271109117979, 0.109599073571517284723367371434",
        "0.15957111962929327, 3, -1.5040299185863266, 0.177232896878715036284091730209",
        "0.6274762006071, 2, -1.8900832371708143, 0.16107819854397531562468793418",
        "0.15835093350103202, 8, -1.8105305390080373, 0.0223622782357583436150526476962",
        "0.4020863325429441, 8, -1.9460205265156678, 0.0158899910202388709972855269772",
        "0.027206299967774683, 7, -1.3487737859632958, 0.0720911057995967084088190218744",
        "0.2114182804277419, 9, -1.6170554352797573, 0.0275823068722209610999281694504",
        "0.7544610759032985, 5, -1.315519083586279, 0.100045917270340348530929277985",
        "0.8921409978504445, 8, -1.2132694959818489, 0.0705666589217701488242734398135",
        "0.4654491288553382, 2, -2.289586863934508, 0.1266834731010717058657923331",
        "0.4436892581835531, 7, -3.5262279185435896, 0.000843093977418368485525042447597",
        "0.08676426129231507, 2, -2.6533108245033152, 0.142016169833697526948579073763",
        "0.9817827670033764, 4, -2.8458842299469387, 0.0103591383991516513477974616357",
        "0.873117528079184, 1, -3.326452985581693, 0.123972468620333398133689669714",
        "0.3567206566799368, 9, -2.3635406177955227, 0.00506652016510339424619581920548",
        "0.5838798513149249, 5, -2.084813610295052, 0.0277189979098831710732739418813",
        "0.23098502217341754, 8, -3.912197988780085, 0.000262162919372284302454773427694",
        "0.9153539243770782, 3, -2.1627117171368573, 0.0522404797936984487211733758968",
        "0.20019785769492826, 3, -2.0084425246699986, 0.0966899580599961987380126156015",
        "0.9351316762425306, 1, -4.211509518101879, 0.0620177359699540807659043730553",
        "0.7468381375280576, 3, -4.953130039553916, 0.00144066495635376975729377045573",
        "0.5624227512181915, 7, -6.654791287968949, 1.42132945172462221431994470653e-06",
        "0.5961152345510518, 1, -6.814195540754605, 0.0413314405086557436168062823658",
        "0.8503191145015174, 8, -4.056455385943806, 0.000144113062411118320373390114826",
        "0.8281902580758295, 8, -6.554861857003606, 6.30872760258769986034209042411e-07",
        "0.33961868190419264, 4, -6.130668257509003, 0.000123593574596111610566299338884",
        "0.2323618875805521, 3, -4.708907767269781, 0.00398765356220288172915136349508",
        "0.8553502660262496, 9, -4.099931361857671, 8.43359499875100413249661539051e-05",
        "0.7842919145273811, 2, -6.524361482095172, 0.00125464282914482518133698633117",
        "0.6741388432726885, 9, -13.55898876995823, 4.32607131584136980549794158646e-14",
        "0.92667845521276, 8, -15.27695267340874, 2.99471515717937876797479313646e-15",
        "0.11171434317842655, 1, -8.006278653034787, 0.428317238074598940005185516952",
        "0.1442627935209828, 5, -9.855493492422891, 9.76246209036406493999630963383e-08",
        "0.5815006036710888, 8, -11.094949031943713, 4.38695786563372636969677849668e-11",
        "0.2266772900180679, 9, -13.80626483771978, 4.74600922920426809114450670609e-14",
        "0.32108967780741415, 4, -8.07095053037433, 7.41580127655544308200314595677e-06",
        "0.6975178280336111, 5, -10.413054615313133, 1.35208379857343151464809899542e-08",
        "0.7623634016395487, 8, -9.285346491599615, 1.76780754806086490498540057253e-09",
        "0.7071088650526037, 7, -15.970695628089075, 6.85055726814425650771424877154e-15",
        "0.7544057018807351, 1, -21.25607563717859, 6.46775579784547243769832095937e-06",
        "0.9602830021196862, 5, -25.9865225545419, 7.13690158981107341550766046314e-21",
        "0.44810426047873697, 9, -28.816434364044508, 7.83509253492022613205240948544e-29",
        "0.3428222124829583, 8, -17.70823799908463, 4.84362289444936682642537054766e-17",
        "0.434898437319612, 7, -21.888883474427146, 8.48715567452197979598585948074e-20",
        "0.6604608169814894, 6, -29.829138446741, 2.72632036946382894852622877915e-25",
        "0.5644609072621231, 9, -28.95153638762003, 4.05855984522902712275434148899e-29",
        "0.4768905656467315, 7, -30.596639219471125, 1.84957967284753363187059400419e-27",
        "0.25493581216338623, 2, -28.58794856737284, 8.02819866165411541339754454829e-11",
        "0.9112765515000535, 6, -30.599996030403638, 2.03948364812262740818278983803e-26",
        // for b in [0, 1, 2, 3, 4]:
        //   for i in range(10):
        //     x, y, z = -random.random(), random.randint(1, 9), -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "-0.5188420581899112, 8, -1.0500529218215608, 0.120861083456608767239165165384",
        "-0.0023410578177983643, 4, -1.6389464276880181, 0.103198321253832220159148752105",
        "-0.27403337319731147, 2, -1.9080876739560864, 0.352956684104427431018013950154",
        "-0.15298097848210102, 8, -1.6442534855912778, 0.0337970403046898354295091191801",
        "-0.5194649435141258, 7, -1.206745461499568, 0.104855755326681217808713746299",
        "-0.8504386140231649, 1, -1.950287030470688, 40.6760585542427043757786663944",
        "-0.6955858634199227, 7, -1.2974934957578357, 0.0917204461264590057898813258717",
        "-0.13473082039007456, 9, -1.617294439081989, 0.0293293958302977747870257525432",
        "-0.22185747955445412, 8, -1.013197029706824, 0.125131654760114159161328664626",
        "-0.7790093939891133, 9, -1.977984975930657, 0.0154986465192169632307538189352",
        "-0.6333525154050418, 3, -2.7842314360878095, 0.0908501276837325698797569946968",
        "-0.21668084060371784, 7, -2.1779172486865157, 0.0154592200635505566985378276396",
        "-0.8975504437181289, 3, -3.6947948453561716, 0.0642088916579819425943032542268",
        "-0.49889748096152875, 5, -2.0206789636259552, 0.0478467267649844019327526094389",
        "-0.07398999117624916, 1, -3.128824216603186, 1.27190204726572794485783913342",
        "-0.3678325081532181, 4, -3.570084858170816, 0.0100036921783155965357242306713",
        "-0.5482477289864353, 3, -3.867749160009468, 0.0311603452331619409866210535656",
        "-0.98074196649013, 9, -2.066946342895003, 0.0135270516062240026316809662201",
        "-0.6249272750699911, 8, -3.871104369096746, 0.000437306625041819296391682107705",
        "-0.7454875196000394, 1, -2.5178104704797857, 31.35536727317684425152013266",
        "-0.0324823117452383, 1, -4.1656229751443075, 1.14746529468640288733328562276",
        "-0.7021444755885952, 7, -7.565787729421237, 8.98416513637445939188070701205e-07",
        "-0.8178964452629013, 1, -7.879313422727794, 673260.250542059884618114594589",
        "-0.2922542790617445, 6, -7.989197457257395, 9.04596294287084365568343489553e-07",
        "-0.5482802194307371, 6, -5.901341542295971, 4.50260772487831990673389907768e-05",
        "-0.028184387539250677, 5, -4.062284997569072, 0.00148100835626900419547088696805",
        "-0.3695930248206172, 8, -5.085090764235662, 3.25210703664702609075831856881e-05",
        "-0.6836798764676074, 1, -5.7172675046285235, 720.955865443291975896515663688",
        "-0.3093342027255367, 7, -6.895096513109257, 2.0338047237974370662327189e-06",
        "-0.3036732600372283, 4, -6.646355808966447, 0.00016842098734650267162497874098",
        "-0.8628773923957973, 5, -14.907171263447154, 6.40872232816379320209750529793e-10",
        "-0.1666721107934065, 6, -8.156968143895135, 5.65515375083533659646549442841e-07",
        "-0.7999520727028387, 9, -12.752399097175212, 2.2215121036027045193436875412e-12",
        "-0.7914423148739238, 6, -8.706980194917293, 5.74850760869906596612234562753e-07",
        "-0.7464121732205845, 7, -13.075102365362515, 3.89518880841322545283150065665e-11",
        "-0.5780677801488484, 4, -10.222210519857944, 3.45593244154498341156958279216e-06",
        "-0.5511541650532694, 1, -11.096987389871149, 7255.59428241600993346776437049",
        "-0.39411803197811635, 9, -14.599346366049824, 2.25215557900376651258397276457e-14",
        "-0.22820151012977974, 5, -11.648016318132889, 1.24365487905545817243715814996e-08",
        "-0.4141746041969817, 9, -13.210931527939554, 4.61178704881929470293410093383e-13",
        "-0.7021896292658335, 6, -31.415531417673158, 1.78700805032173845163236124495e-23",
        "-0.7408043449942696, 6, -22.37777033672218, 7.3694121855804336834635742458e-17",
        "-0.04588779442303004, 4, -27.018501691681855, 7.38976925921560037132451943554e-17",
        "-0.94329097086466, 5, -19.79552461004343, 9.13859462817891680372484675451e-13",
        "-0.15010832994612455, 4, -29.219766384707192, 7.82224343869722451820543480013e-18",
        "-0.7739380955029117, 7, -27.50322933315506, 1.43405488177914371954341359517e-22",
        "-0.8564736217418897, 4, -31.739927089617012, 1.62936066403145359522427840024e-16",
        "-0.6489201154598017, 9, -26.86621377846609, 1.72292049411954350510275570369e-25",
        "-0.8789212776641596, 8, -16.769769282578714, 5.04665773360577097399204057959e-15",
        "-0.6960984140690415, 5, -28.742960890645374, 6.03837080548376086791980931957e-19",
        // Hit use of large z above 2^26 when x+y < 1. Requires x to be negative.
        // for b in [40, 48]:
        //   for i in range(10):
        //     x, y, z = -random.random()*2**-55, 1, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "-3.702572446955424e-18, 1, -1674898990276.807, 1.00000620145408176916923331654",
        "-1.2489310044774175e-17, 1, -1112986110849.0952, 1.00001390052522532727466392723",
        "-1.3668877734033251e-17, 1, -1776203173505.0132, 1.00002427899873956123693050753",
        "-1.8669950000218342e-19, 1, -1587019522060.5947, 1.00000029629579515806625859812",
        "-9.23064327323242e-18, 1, -1354097954792.016, 1.00001249927329296421021094687",
        "-2.294534557979883e-17, 1, -1757647728278.8423, 1.00004033064779161592448731205",
        "-1.822349749858894e-18, 1, -1273086997455.8228, 1.00000232001246258696409568169",
        "-1.3706347083007253e-17, 1, -1673243376046.0671, 1.00002293431745386785960591452",
        "-1.3368371468959976e-17, 1, -2047019907358.5977, 1.00002736569695977965621540205",
        "-2.2281411221492123e-17, 1, -1768384250452.335, 1.0000394028729547456244409143",
        "-2.2628822652323836e-17, 1, -528473515082598.25, 1.01203052499534233454238565091",
        "-1.5491918794989766e-17, 1, -404923626868302.5, 1.00629276069223318894563500002",
        "-5.308383788647382e-18, 1, -547728474487075.75, 1.00291178398629534588250326346",
        "-2.421729450202671e-17, 1, -318297125482627.4, 1.00773808061695497265449603121",
        "-1.3196502909214093e-17, 1, -547696346330428.44, 1.00725385912299202526756767346",
        "-3.0096956928695905e-18, 1, -489489435616076.56, 1.00147429995928357377043421652",
        "-1.9261240038423732e-17, 1, -428070729867071.9, 1.00827925813561915484144399599",
        "-1.5969392870282128e-17, 1, -370282637311426.75, 1.00593070632070488188967634583",
        "-1.5199041953018484e-17, 1, -492608817179268.75, 1.00751528111021959686921983214",
        "-1.5595868760284054e-17, 1, -446403666757161.25, 1.00698634443187000084699067345",
        // Expand range of x from the above cases. Requires smaller z.
        // for b in [24, 32]:
        //   for i in range(10):
        //     x, y, z = -random.random()*2**-50, 1, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "-7.974584793469102e-16, 1, -19699105.280573938, 1.00000001570921866493096123982",
        "-8.026455043481995e-17, 1, -24823303.583062746, 1.00000000199243130438647783532",
        "-5.147871492252243e-16, 1, -23933356.187688615, 1.00000001232058427915050446032",
        "-4.2813842429975754e-16, 1, -23967439.3450597, 1.00000001026138176834185263574",
        "-6.114741185384639e-16, 1, -30516011.433849774, 1.00000001865975136691614049818",
        "-4.650329216481511e-16, 1, -29674009.78627184, 1.0000000137993915631374268198",
        "-2.850395425208763e-16, 1, -17755129.083347544, 1.00000000506091388412291496414",
        "-2.9242286891262733e-16, 1, -28928767.93832121, 1.00000000845943335041252696023",
        "-3.4978884696742346e-16, 1, -26106344.021077476, 1.00000000913170801536162672764",
        "-5.69143704732998e-16, 1, -25030300.504482806, 1.00000001424583806117352784594",
        "-2.925153513201412e-16, 1, -5872567063.678579, 1.00000171781749322981016725909",
        "-1.2957789416293734e-17, 1, -4818811562.80563, 1.00000006244114741708502083666",
        "-4.620183002706999e-16, 1, -4625668357.73638, 1.00000213714571594997411154178",
        "-3.6762102420640106e-16, 1, -6354739041.704281, 1.0000023361384038428028203675",
        "-1.710927723049339e-16, 1, -6486103553.55575, 1.00000110972605418023881643394",
        "-4.1579410106948047e-16, 1, -8451853512.751289, 1.00000351423700862176353662436",
        "-4.928804069909112e-16, 1, -7884677418.006172, 1.0000038862105660756482621878",
        "-2.664583753813762e-16, 1, -4759176183.949682, 1.00000126812315819646240352203",
        "-1.7942936553513716e-16, 1, -8378615994.66663, 1.00000150337088204645757426888",
        "-8.4730625733622295e-16, 1, -8227380882.593294, 1.00000697113560156258121837486",
        // Expand range of x from the above cases. Requires smaller z.
        // for b in [24, 27]:
        //   for i in range(10):
        //     x, y, z = -random.random()*2**-20, 1, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "-5.336700990207014e-07, 1, -28529620.880939238, 4095542.16457487743389145857328",
        "-2.7676421202203157e-07, 1, -28130079.167230584, 2405.22870149659197702489807564",
        "-7.613595006311205e-07, 1, -27063715.508285794, 888659556.090259928481641010396",
        "-3.5773026208256195e-07, 1, -29473932.179690637, 37938.4113207743156836825222755",
        "-2.917027873778877e-07, 1, -32149852.573560797, 11827.740926435736038383373406",
        "-8.667881953214492e-07, 1, -23129692.436043844, 509302679.968911875253430552335",
        "-3.172422637597331e-07, 1, -20956467.9262387, 771.455029096370986627210454419",
        "-2.94163329526159e-07, 1, -28836934.477135003, 4830.8122817472372727782397218",
        "-1.6671340878181233e-07, 1, -33468580.995897163, 264.981944583651617763446781623",
        "-1.137060204187555e-07, 1, -26860834.129701048, 21.2050365729885905890455679787",
        "-3.8071855643124373e-07, 1, -138877796.92460406, 91750916643336940591812.368269",
        "-7.144290779980769e-07, 1, -237351433.70019686, 4.40224538072624391145386905457e+73",
        "-5.53495628731812e-08, 1, -206874422.71670645, 93939.776828697774706855543027",
        "-6.554760561797986e-07, 1, -261272727.04967517, 2.37922774964733071853460387381e+74",
        "-4.050519338086029e-07, 1, -244389450.5653385, 9.79493775201041173565577576772e+42",
        "-5.420274166289836e-07, 1, -194759724.15162852, 7.02032539937607881315931112861e+45",
        "-4.7712708433920074e-08, 1, -229813022.59335348, 57814.911331311174734247956045",
        "-9.339478060193397e-07, 1, -228198287.02257952, 3.62393506751463306503493824559e+92",
        "-7.044476791082879e-07, 1, -186811675.65274024, 1.42154938871266502749212296616e+57",
        "-4.937038377064631e-07, 1, -176264618.48152652, 6.21453956097956128821371375931e+37",
        // Hit use of large z above 2^26 when x+y > 1.
        // for b in [40, 48]:
        //   for i in range(10):
        //     x, y, z = random.random()*2**-55, 1, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "8.50000344713563e-18, 1, -1503553918589.503, 0.999987219868175615472678877174",
        "2.3231261435036756e-17, 1, -2147909148059.897, 0.999950102605959091746591465354",
        "1.4549261942930887e-17, 1, -1242943727610.0427, 0.999981916249639663459099702964",
        "1.9042098252583027e-17, 1, -1960140489309.1226, 0.999962675508789864585783010343",
        "7.378927048518629e-18, 1, -1472474580134.4072, 0.999989134776519358264755831606",
        "2.2058036738835262e-17, 1, -1776727944979.6936, 0.999960809637677036840057696291",
        "1.573843750569005e-17, 1, -1261202666189.2393, 0.999980150837653400713454317164",
        "5.366100998998383e-18, 1, -1920149013424.7148, 0.999989696339543910062936514021",
        "1.8538219572890425e-17, 1, -2167259291484.9912, 0.999959823678473087376232195636",
        "1.3834494283321978e-17, 1, -1985132567315.05, 0.999972537071958956566507585896",
        "1.1765340210219243e-17, 1, -350133269074388.1, 0.995889036212155620184788647251",
        "1.4196147304036987e-17, 1, -486348611865121.4, 0.99311950322326986626168972022",
        "1.5954731986869044e-17, 1, -475578779310339.4, 0.992440982204595879894854514368",
        "1.867960662698333e-17, 1, -490560297042445.1, 0.99087838329817476002613126854",
        "2.3238691710036644e-17, 1, -396938131689481.6, 0.990818090683627490490919130376",
        "4.121583487173531e-18, 1, -533984425611730.94, 0.997801558728647196186676640902",
        "2.7388916012363687e-17, 1, -372457178329412.4, 0.98985065736864106661619808884",
        "2.17491725265551e-17, 1, -373165841539423.7, 0.991916797930578266605938585218",
        "1.468579563144174e-17, 1, -404638756782903.7, 0.994075179300333681494935868393",
        "2.3705911911434888e-18, 1, -335499028421195.0, 0.999204985150488945451559998623",
        // Expand range of x from the above cases. Requires smaller z.
        // Hit use of large z above 2^26 when x+y > 1.
        // for b in [24, 32]:
        //   for i in range(10):
        //     x, y, z = random.random()*2**-50, 1, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "6.354876984945542e-16, 1, -31314992.39216908, 0.99999998009970775433682026579",
        "8.618156537619779e-16, 1, -19182954.42025781, 0.999999983467829731875626244765",
        "1.138503879417834e-16, 1, -24904787.00849943, 0.999999997164580341494655225311",
        "2.1351448831057339e-16, 1, -30894112.47901269, 0.999999993403659403990163952071",
        "1.5190501268800374e-16, 1, -21318724.32582101, 0.999999996761578916037791434668",
        "4.707029206011926e-16, 1, -19115778.593199864, 0.999999991002147227095745427837",
        "7.501173322378474e-16, 1, -19508972.53269092, 0.999999985365981676153697557202",
        "2.4682439564913957e-16, 1, -20779860.559955515, 0.99999999487102348876878973864",
        "5.8147098774875e-16, 1, -23137490.976582304, 0.999999986546220358320993177918",
        "4.730929578922396e-16, 1, -23002701.172543094, 0.999999989117584121983790201069",
        "1.4865716575217705e-16, 1, -6060116195.040617, 0.999999099120708108388152086344",
        "3.946346661409329e-16, 1, -4558627012.760921, 0.999998201009368943642301375443",
        "5.871803697692237e-16, 1, -6742109545.591936, 0.999996041173460169183534807774",
        "5.280937164205083e-16, 1, -7738994260.106953, 0.999995913094111245764427048258",
        "2.416750990658052e-16, 1, -5538220892.099985, 0.999998661550812977336643988327",
        "4.231463890023936e-16, 1, -4751810264.957883, 0.999997989290666637286846248478",
        "5.409705405857083e-16, 1, -6005630146.219, 0.999996751136290811212661401309",
        "2.382487849297694e-16, 1, -7194584371.808722, 0.999998285900484412117533319231",
        "4.4371038065451513e-16, 1, -6283695102.712627, 0.999997211863140919265388715379",
        "7.837828960300255e-16, 1, -4986232752.421234, 0.999996091883689733965527056365",
        // Expand range of x from the above cases. Requires smaller z.
        // for b in [24, 27]:
        //   for i in range(10):
        //     x, y, z = random.random()*2**-20, 1, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "1.66346366297373e-07, 1, -25711219.87494736, 0.0138847016795655432126235839339",
        "2.7528269003468857e-07, 1, -20107864.609775133, 0.00394484230550316278669470446481",
        "2.6832589208868387e-07, 1, -30821170.686701678, 0.000256055328423194586068559147143",
        "5.25372753127701e-07, 1, -28044551.568554856, 3.99185587177799478327263621852e-07",
        "5.508912331154624e-07, 1, -32150997.995614335, 2.03192311290848505640750874593e-08",
        "2.2686985210847656e-07, 1, -28484752.549499847, 0.00156115262994557778317982563341",
        "1.7485703317858038e-07, 1, -32813990.364042304, 0.00322198870389374219082526698619",
        "6.01625085916704e-07, 1, -29275581.157029256, 2.24288282367491917114102477619e-08",
        "5.5443680365122744e-08, 1, -21822713.36642589, 0.298217703348955783842778345515",
        "4.634621479389429e-07, 1, -25250341.336813852, 8.27249236324031593552099538999e-06",
        "4.4693507584930594e-08, 1, -246487798.51535332, 1.64299602321725422967085589099e-05",
        "5.147415948496822e-07, 1, -260057118.84939027, 7.31801290688998262498400555403e-59",
        "7.789608216342709e-07, 1, -231965309.9828784, 3.36155702700598176082401954764e-79",
        "3.986569984798169e-07, 1, -193341331.60672498, 3.35695255580539730442370084044e-34",
        "1.620022038558036e-07, 1, -221924386.5060208, 2.43299910697123913334841331888e-16",
        "8.479057402680049e-07, 1, -228850690.23545274, 5.3441433211051358519175110011e-85",
        "9.196504821454351e-07, 1, -248728460.9458745, 4.55108401250702368217369093907e-100",
        "1.298214946032042e-07, 1, -168225510.75086904, 3.27580981894164293733928631651e-10",
        "3.526953995229003e-07, 1, -201489823.3980894, 1.37110454425974507093157286134e-31",
        "5.761339954571439e-07, 1, -141326512.58398587, 4.3495763259420070111573350333e-36",
        // Handle underflow and overflow. No support for -infinity. Use MAX_VALUE.
        "-0.000000000001, 1, -1.7976931348623157E308, Infinity",
        "0, 1, -1.7976931348623157E308, 1.0",
        "0.000000000001, 1, -1.7976931348623157E308, 0",
        "0.999999999999, 1, -1.7976931348623157E308, 0",
        "1.000000000001, 1, -1.7976931348623157E308, 0",
        // Rounding of 1+/-x = 1
        "-1e-20, 1, -1.7976931348623157E308, Infinity",
        "1e-20, 1, -1.7976931348623157E308, 0",
        // Handle x+n = 0
        "-1, 1, -2, Infinity",
        "-1, 1, -1.7976931348623157E308, Infinity",
        // Trigger the use of (x+n)^z * exp(z*t)
        "1e-200, 1, -1.7976931348623157E308, 0",
    })
    @Order(1)
    void testPowxpn(double x, int n, double z, double expected) {
        // Check data is suitable for the zeta use case.
        // We expand the range to include x == -1 and x + n == 0 to check edge cases.
        Assertions.assertTrue(x >= -1, "x should be > -1");
        Assertions.assertTrue(n >= 1, "n should be >= 1");
        Assertions.assertTrue(x + n >= 0 && Double.isFinite(x + n), "(x+n) should be > 0 and finite");
        Assertions.assertTrue(z < -1 && z > Double.NEGATIVE_INFINITY, "z should be < -1 and finite");

        // Math.pow is not expected to be ulp exact for all cases.
        // When it is we collect the error separately to obtain the RMS for this subset.
        final double lo = DD.ofSum(x, n).lo();
        final long stdUlp;
        if (lo == 0) {
            // Expect standard pow to be 1 ulp (it approaches 0 ulp RMS)
            stdUlp = TestUtils.assertEquals(expected, Math.pow(x + n, z), 1,
                ERROR_POW::add, () -> String.format("%s, %s, %s", x, n, z));
        } else {
            // Skip edge cases of x+1 == 1 which have round-off
            final double xpn = x + n;
            if (xpn == 1) {
                // Flag to ignore
                stdUlp = Long.MAX_VALUE;
                // Check the edge case is correct
                // Expected 1 +/= tiny. nextAfter requires the direction > start so scale by 2.
                final double xx = Math.nextAfter(1.0, 2 * Math.signum(lo));
                if (expected == 0 || expected == Double.POSITIVE_INFINITY) {
                    Assertions.assertEquals(expected, Math.pow(xx, z));
                }
            } else {
                // This can be very wrong
                stdUlp = TestUtils.assertEquals(expected, Math.pow(xpn, z), 10000000000L,
                    ERROR_POWP::add, () -> String.format("%s, %s, %s", x, n, z));
            }
        }
        // Always 1 ulp
        final double actual = powxpn(x, n, z);
        final long extUlp = TestUtils.assertEquals(expected, actual, 1,
            ERROR_POW_PLUS::add, () -> String.format("%s, %s, %s", x, n, z));

        // To check the extended precision function is better than Math.pow
        if (Math.abs(stdUlp) <= 1) {
            ERROR_POW1.add(stdUlp);
            ERROR_POW_PLUS1.add(extUlp);
        }

        // Ensure the final implementation in the HurwitzZeta function is the same when
        // the domain is applicable. Collect error statistics to compare to Math.pow.
        // Ignore edge case where x+1 == 1 with round-off.
        if (x > 0 && n >= 1 && z >= -EXT_POW_LIMIT && stdUlp != Long.MAX_VALUE) {
            Assertions.assertEquals(actual, HurwitzZeta.extPowxpn(x, n, z), "Final implementation differs");
            ERROR_FINAL_STD.add(stdUlp);
            ERROR_FINAL_EXT.add(extUlp);
        }
    }

    @Test
    @Order(2)
    void testExtPowRMS() {
        // The Math.pow function is very accurate when there is no round-off
        Assertions.assertTrue(ERROR_POW.getRMS() < 0.1);
        // Some cases are very wrong using Math.pow(x+n, z) when there is round-off
        Assertions.assertTrue(ERROR_POWP.getRMS() > 1e6);
        Assertions.assertTrue(ERROR_FINAL_STD.getRMS() > 1e6);

        // Extended precision is ulp accurate for all x+n cases
        Assertions.assertTrue(ERROR_POW_PLUS.getRMS() < 0.5);
        Assertions.assertTrue(ERROR_FINAL_EXT.getRMS() < 0.5);

        // Where Math.pow was ulp accurate, the extended precision lowered the RMS.
        // Note there are some cases where Math.pow is closer to the correct result
        // but on average the extended precision method is closer.
        Assertions.assertTrue(ERROR_POW_PLUS1.getRMS() < ERROR_POW1.getRMS() * 0.75);
    }
}
