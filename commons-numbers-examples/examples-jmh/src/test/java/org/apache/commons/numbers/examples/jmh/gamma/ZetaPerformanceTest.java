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
package org.apache.commons.numbers.examples.jmh.gamma;

import java.util.stream.Stream;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * Tests the zeta implementation in {@link ZetaPerformance}.
 *
 * <p>Special cases and spot tests have been adapted from the HurwitzZetaTest.
 */
class ZetaPerformanceTest {

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
        Assertions.assertEquals(z, ZetaPerformance.zeta(s, a));
    }

    /**
     * Spot tests for the zeta function to check various points in the domain and extreme values.
     */
    @ParameterizedTest
    @MethodSource(value = "testZetaSpot")
    void testZetaSpot(double s, double a, double z, int ulp) {
        final double actual = ZetaPerformance.zeta(s, a);
        if (!Double.isFinite(z)) {
            Assertions.assertEquals(z, actual);
        } else {
            // Double the ULP for positive a as the performance implementation is less accurate.
            // Not required for many cases but is simpler than adjusting copied ULPs.
            if (a > 0) {
                ulp = ulp * 2 + 1;
            }
            Assertions.assertEquals(z, actual, ulp * Math.ulp(z));
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
            Arguments.of(1.001, 1, 1000.57728847601162684806668989, 0),
            Arguments.of(1.001, 3, 999.077634929516400548962369402, 0),
            Arguments.of(1.001, 156.78, 994.961087152463264230812967825, 0),
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
            Arguments.of(133, -123.5, 3.36756171921967259297351128082e-279, 0),
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

}
