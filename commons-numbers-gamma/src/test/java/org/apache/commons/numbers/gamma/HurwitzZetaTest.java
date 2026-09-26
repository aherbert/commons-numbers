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
import java.util.Arrays;
import java.util.SplittableRandom;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleSupplier;
import java.util.function.IntSupplier;
import java.util.stream.Stream;
import org.apache.commons.numbers.core.DD;
import org.apache.commons.numbers.core.DDMath;
import org.apache.commons.numbers.fraction.BigFraction;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.MethodOrderer;
import org.junit.jupiter.api.Order;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestMethodOrder;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * Test the {@link HurwitzZeta} function.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class HurwitzZetaTest {
    /** Seed for data generation. */
    private static final long SEED = 6516587940839803692L;
    /** Table used to create a histogram of the number of steps to converge the tail series.
     * Used in the {@link #zeta(double, double, int, int)} implementation. */
    private static final int[] M = new int[53];
    /** Optimal N used to test the zeta function. */
    private static final int N = 8;
    /** Minimum N used to test the zeta function.
     * Used for reporting RMS errors with varying N. When MIN_N >= MAX_N no report is printed. */
    private static final int MIN_N = N + 1; // e.g. 5
    /** Maximum N used to test the zeta function. Used for reporting RMS errors with varying N. */
    private static final int MAX_N = N; // e.g. 12
    /** Filenames of resources used for the test zeta function. */
    private static final String[] TEST_RESOURCES = {
        "hzeta_s1_4_a1_8.csv",
        "hzeta_s1_4_a8_32.csv",
        "hzeta_s1_4_a32_2147483648.csv",
        "hzeta_s1_4_a0_1.csv",
        "hzeta_s1_4_a1e-16_1e-14.csv",
        "hzeta_s4_32_a1_8.csv",
    };
    /** Flag set when the JVM version is printed. Used for testing. */
    private static boolean jvm = false;

    /**
     * Numerators of the even Bernoulli numbers {@code B_{2k}}.
     * Taken from:
     * <pre>
     * A000367 Numerators of Bernoulli numbers B_2n.
     * https://oeis.org/A164020/b164020.txt
     * </pre>
     *
     * <p>Contains the sequence up to 2k = 106 required for M=53 in the zeta implementation.
     * Johansson (2015) suggests N ~ M ~ P for P-bits of precision.
     */
    private static final String[] NUM = {
        "0 1",
        "1 1",
        "2 -1",
        "3 1",
        "4 -1",
        "5 5",
        "6 -691",
        "7 7",
        "8 -3617",
        "9 43867",
        "10 -174611",
        "11 854513",
        "12 -236364091",
        "13 8553103",
        "14 -23749461029",
        "15 8615841276005",
        "16 -7709321041217",
        "17 2577687858367",
        "18 -26315271553053477373",
        "19 2929993913841559",
        "20 -261082718496449122051",
        "21 1520097643918070802691",
        "22 -27833269579301024235023",
        "23 596451111593912163277961",
        "24 -5609403368997817686249127547",
        "25 495057205241079648212477525",
        "26 -801165718135489957347924991853",
        "27 29149963634884862421418123812691",
        "28 -2479392929313226753685415739663229",
        "29 84483613348880041862046775994036021",
        "30 -1215233140483755572040304994079820246041491",
        "31 12300585434086858541953039857403386151",
        "32 -106783830147866529886385444979142647942017",
        "33 1472600022126335654051619428551932342241899101",
        "34 -78773130858718728141909149208474606244347001",
        "35 1505381347333367003803076567377857208511438160235",
        "36 -5827954961669944110438277244641067365282488301844260429",
        "37 34152417289221168014330073731472635186688307783087",
        "38 -24655088825935372707687196040585199904365267828865801",
        "39 414846365575400828295179035549542073492199375372400483487",
        "40 -4603784299479457646935574969019046849794257872751288919656867",
        "41 1677014149185145836823154509786269900207736027570253414881613",
        "42 -2024576195935290360231131160111731009989917391198090877281083932477",
        "43 660714619417678653573847847426261496277830686653388931761996983",
        "44 -1311426488674017507995511424019311843345750275572028644296919890574047",
        "45 1179057279021082799884123351249215083775254949669647116231545215727922535",
        "46 -1295585948207537527989427828538576749659341483719435143023316326829946247",
        "47 1220813806579744469607301679413201203958508415202696621436215105284649447",
        "48 -211600449597266513097597728109824233673043954389060234150638733420050668349987259",
        "49 67908260672905495624051117546403605607342195728504487509073961249992947058239",
        "50 -94598037819122125295227433069493721872702841533066936133385696204311395415197247711",
        "51 3204019410860907078243020782116241775491817197152717450679002501086861530836678158791",
        "52 -319533631363830011287103352796174274671189606078272738327103470162849568365549721224053",
        "53 36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901",
    };

    /**
     * Denominators of the even Bernoulli numbers {@code B_{2k}}.
     * Taken from:
     * <pre>
     * A002445 Denominators of Bernoulli numbers B_{2n}.
     * https://oeis.org/A002445/b002445.txt
     * </pre>
     */
    private static final String[] DENOM = {
        "0 1",
        "1 6",
        "2 30",
        "3 42",
        "4 30",
        "5 66",
        "6 2730",
        "7 6",
        "8 510",
        "9 798",
        "10 330",
        "11 138",
        "12 2730",
        "13 6",
        "14 870",
        "15 14322",
        "16 510",
        "17 6",
        "18 1919190",
        "19 6",
        "20 13530",
        "21 1806",
        "22 690",
        "23 282",
        "24 46410",
        "25 66",
        "26 1590",
        "27 798",
        "28 870",
        "29 354",
        "30 56786730",
        "31 6",
        "32 510",
        "33 64722",
        "34 30",
        "35 4686",
        "36 140100870",
        "37 6",
        "38 30",
        "39 3318",
        "40 230010",
        "41 498",
        "42 3404310",
        "43 6",
        "44 61410",
        "45 272118",
        "46 1410",
        "47 6",
        "48 4501770",
        "49 6",
        "50 33330",
        "51 4326",
        "52 1590",
        "53 642",
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses {@code 2k!} divided by Bernoulli number {@code B_2k}.
     * The table size is suitable for N ~ M ~ P for P-bits of precision (53 entries)
     * as stated in Johansson (2015) section 3.1. In practice the result in double
     * precision requires lower N & M values.
     */
    private static final double[] F = {
        12.0, // 2! / (1 / 6)
        -720.0, // 4! / (-1 / 30)
        30240.0, // 6! / (1 / 42)
        -1209600.0, // 8! / (-1 / 30)
        4.790016E7, // 10! / (5 / 66)
        -1.8924375803183792E9, // 12! / (-691 / 2730)
        7.47242496E10, // 14! / (7 / 6)
        -2.950130727918164E12, // 16! / (-3617 / 510)
        1.1646782814350067E14, // 18! / (43867 / 798)
        -4.597978722407473E15, // 20! / (-174611 / 330)
        1.81521054019435456E17, // 22! / (854513 / 138)
        -7.1661652561756672E18, // 24! / (-236364091 / 2730)
        2.82908877253043E20, // 26! / (8553103 / 6)
        -1.1168794925000445E22, // 28! / (-23749461029 / 870)
        4.4092635141854666E23, // 30! / (8615841276005 / 14322)
        -1.7407074646225822E25, // 32! / (-7709321041217 / 510)
        6.872037622739274E26, // 34! / (2577687858367 / 6)
        -2.7129717107520044E28, // 36! / (-26315271553053477373 / 1919190)
        1.0710383014704457E30, // 38! / (2929993913841559 / 6)
        -4.228289733582729E31, // 40! / (-261082718496449122051 / 13530)
        1.669261878547101E33, // 42! / (1520097643918070802691 / 1806)
        -6.589981753232787E34, // 44! / (-27833269579301024235023 / 690)
        2.6016205165923062E36, // 46! / (596451111593912163277961 / 282)
        -1.0270786120209628E38, // 48! / (-5609403368997817686249127547 / 46410)
        4.054743835786749E39, // 50! / (495057205241079648212477525 / 66)
        -1.6007487042788347E41, // 52! / (-801165718135489957347924991853 / 1590)
        6.319502582715391E42, // 54! / (29149963634884862421418123812691 / 798)
        -2.4948396201225357E44, // 56! / (-2479392929313226753685415739663229 / 870)
        9.849232037909393E45, // 58! / (84483613348880041862046775994036021 / 354)
        -3.888320954748035E47, // 60! / (-1215233140483755572040304994079820246041491 / 56786730)
        1.535047584313167E49, // 62! / (12300585434086858541953039857403386151 / 6)
        -6.060124957607529E50, // 64! / (-106783830147866529886385444979142647942017 / 510)
        2.3924414381101893E52, // 66! / (1472600022126335654051619428551932342241899101 / 64722)
        -9.444980218768351E53, // 68! / (-78773130858718728141909149208474606244347001 / 30)
        3.728728733414322E55, // 70! / (1505381347333367003803076567377857208511438160235 / 4686)
        -1.4720431007109737E57, // 72! / (-5827954961669944110438277244641067365282488301844260429 / 140100870)
        5.811393226148101E58, // 74! / (34152417289221168014330073731472635186688307783087 / 6)
        -2.2942460864500874E60, // 76! / (-24655088825935372707687196040585199904365267828865801 / 30)
        9.057320508803927E61, // 78! / (414846365575400828295179035549542073492199375372400483487 / 3318)
        -3.575686814230726E63, // 80! / (-4603784299479457646935574969019046849794257872751288919656867 / 230010)
        1.4116245727459505E65, // 82! / (1677014149185145836823154509786269900207736027570253414881613 / 498)
        -5.5728704383437275E66, // 84! / (-2024576195935290360231131160111731009989917391198090877281083932477 / 3404310)
        2.2000810641991213E68, // 86! / (660714619417678653573847847426261496277830686653388931761996983 / 6)
        -8.685571901589203E69, // 88! / (-1311426488674017507995511424019311843345750275572028644296919890574047 / 61410)
        3.428926346636115E71, // 90! / (1179057279021082799884123351249215083775254949669647116231545215727922535 / 272118)
        -1.3536858624708422E73, // 92! / (-1295585948207537527989427828538576749659341483719435143023316326829946247 / 1410)
        5.344137578373867E74, // 94! / (1220813806579744469607301679413201203958508415202696621436215105284649447 / 6)
        -2.10978095054183E76, // 96! / (-211600449597266513097597728109824233673043954389060234150638733420050668349987259 / 4501770)
        8.329081341920855E77, // 98! / (67908260672905495624051117546403605607342195728504487509073961249992947058239 / 6)
        -3.288189514770133E79, // 100! / (-94598037819122125295227433069493721872702841533066936133385696204311395415197247711 / 33330)
        1.2981251882636475E81, // 102! / (3204019410860907078243020782116241775491817197152717450679002501086861530836678158791 / 4326)
        -5.124792828500739E82, // 104! / (-319533631363830011287103352796174274671189606078272738327103470162849568365549721224053 / 1590)
        2.023187114193683E84, // 106! / (36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901 / 642)
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * This is the inverse of table {@link #M}.
     */
    private static final double[] FM = {
        0.08333333333333333, // (1 / 6) / 2!
        -0.001388888888888889, // (-1 / 30) / 4!
        3.306878306878307E-5, // (1 / 42) / 6!
        -8.267195767195768E-7, // (-1 / 30) / 8!
        2.08767569878681E-8, // (5 / 66) / 10!
        -5.284190138687493E-10, // (-691 / 2730) / 12!
        1.3382536530684679E-11, // (7 / 6) / 14!
        -3.3896802963225827E-13, // (-3617 / 510) / 16!
        8.586062056277845E-15, // (43867 / 798) / 18!
        -2.174868698558062E-16, // (-174611 / 330) / 20!
        5.5090028283602295E-18, // (854513 / 138) / 22!
        -1.3954464685812522E-19, // (-236364091 / 2730) / 24!
        3.534707039629467E-21, // (8553103 / 6) / 26!
        -8.953517427037546E-23, // (-23749461029 / 870) / 28!
        2.267952452337683E-24, // (8615841276005 / 14322) / 30!
        -5.744790668872202E-26, // (-7709321041217 / 510) / 32!
        1.455172475614865E-27, // (2577687858367 / 6) / 34!
        -3.6859949406653103E-29, // (-26315271553053477373 / 1919190) / 36!
        9.336734257095045E-31, // (2929993913841559 / 6) / 38!
        -2.36502241570063E-32, // (-261082718496449122051 / 13530) / 40!
        5.990671762482134E-34, // (1520097643918070802691 / 1806) / 42!
        -1.5174548844682903E-35, // (-27833269579301024235023 / 690) / 44!
        3.843758125454189E-37, // (596451111593912163277961 / 282) / 46!
        -9.736353072646691E-39, // (-5609403368997817686249127547 / 46410) / 48!
        2.466247044200681E-40, // (495057205241079648212477525 / 66) / 50!
        -6.247076741820743E-42, // (-801165718135489957347924991853 / 1590) / 52!
        1.5824030244644914E-43, // (29149963634884862421418123812691 / 798) / 54!
        -4.008273685948936E-45, // (-2479392929313226753685415739663229 / 870) / 56!
        1.0153075855569557E-46, // (84483613348880041862046775994036021 / 354) / 58!
        -2.5718041582418717E-48, // (-1215233140483755572040304994079820246041491 / 56786730) / 60!
        6.514456035233815E-50, // (12300585434086858541953039857403386151 / 6) / 62!
        -1.6501309906896525E-51, // (-106783830147866529886385444979142647942017 / 510) / 64!
        4.179830628539476E-53, // (1472600022126335654051619428551932342241899101 / 64722) / 66!
        -1.058763466770291E-54, // (-78773130858718728141909149208474606244347001 / 30) / 68!
        2.6818791912607708E-56, // (1505381347333367003803076567377857208511438160235 / 4686) / 70!
        -6.793279351107421E-58, // (-5827954961669944110438277244641067365282488301844260429 / 140100870) / 72!
        1.7207577616681404E-59, // (34152417289221168014330073731472635186688307783087 / 6) / 74!
        -4.358730329348894E-61, // (-24655088825935372707687196040585199904365267828865801 / 30) / 76!
        1.1040792903684666E-62, // (414846365575400828295179035549542073492199375372400483487 / 3318) / 78!
        -2.7966655133781345E-64, // (-4603784299479457646935574969019046849794257872751288919656867 / 230010) / 80!
        7.084036501679471E-66, // (1677014149185145836823154509786269900207736027570253414881613 / 498) / 82!
        -1.794407408289224E-67, // (-2024576195935290360231131160111731009989917391198090877281083932477 / 3404310) / 84!
        4.545287063611096E-69, // (660714619417678653573847847426261496277830686653388931761996983 / 6) / 86!
        -1.1513346631982051E-70, // (-1311426488674017507995511424019311843345750275572028644296919890574047 / 61410) / 88!
        2.9163647710923614E-72, // (1179057279021082799884123351249215083775254949669647116231545215727922535 / 272118) / 90!
        -7.387238263497337E-74, // (-1295585948207537527989427828538576749659341483719435143023316326829946247 / 1410) / 92!
        1.8712093117637953E-75, // (1220813806579744469607301679413201203958508415202696621436215105284649447 / 6) / 94!
        -4.739828557761799E-77, // (-211600449597266513097597728109824233673043954389060234150638733420050668349987259 / 4501770) / 96!
        1.2006125993354507E-78, // (67908260672905495624051117546403605607342195728504487509073961249992947058239 / 6) / 98!
        -3.0411872415142924E-80, // (-94598037819122125295227433069493721872702841533066936133385696204311395415197247711 / 33330) / 100!
        7.703417274705106E-82, // (3204019410860907078243020782116241775491817197152717450679002501086861530836678158791 / 4326) / 102!
        -1.951298390909883E-83, // (-319533631363830011287103352796174274671189606078272738327103470162849568365549721224053 / 1590) / 104!
        4.942696565159462E-85, // (36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901 / 642) / 106!
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * This is the same as table {@link #FM} in DECIMAL128 precision.
     */
    private static final BigDecimal[] FD = {
        new BigDecimal("0.08333333333333333333333333333333333"), // (1 / 6) / 2!
        new BigDecimal("-0.001388888888888888888888888888888889"), // (-1 / 30) / 4!
        new BigDecimal("0.00003306878306878306878306878306878307"), // (1 / 42) / 6!
        new BigDecimal("-8.267195767195767195767195767195767E-7"), // (-1 / 30) / 8!
        new BigDecimal("2.087675698786809897921009032120143E-8"), // (5 / 66) / 10!
        new BigDecimal("-5.284190138687493184847682202179557E-10"), // (-691 / 2730) / 12!
        new BigDecimal("1.338253653068467883282698097512912E-11"), // (7 / 6) / 14!
        new BigDecimal("-3.389680296322582866830195391249442E-13"), // (-3617 / 510) / 16!
        new BigDecimal("8.586062056277844564135905450425627E-15"), // (43867 / 798) / 18!
        new BigDecimal("-2.174868698558061873041516423865918E-16"), // (-174611 / 330) / 20!
        new BigDecimal("5.509002828360229515202652608902255E-18"), // (854513 / 138) / 22!
        new BigDecimal("-1.395446468581252334070768626406355E-19"), // (-236364091 / 2730) / 24!
        new BigDecimal("3.534707039629467471693229977803799E-21"), // (8553103 / 6) / 26!
        new BigDecimal("-8.953517427037546850402611318112741E-23"), // (-23749461029 / 870) / 28!
        new BigDecimal("2.267952452337683060310950738868166E-24"), // (8615841276005 / 14322) / 30!
        new BigDecimal("-5.744790668872202445263881987607018E-26"), // (-7709321041217 / 510) / 32!
        new BigDecimal("1.455172475614864901866264867271329E-27"), // (2577687858367 / 6) / 34!
        new BigDecimal("-3.685994940665310178181782479908660E-29"), // (-26315271553053477373 / 1919190) / 36!
        new BigDecimal("9.336734257095044672032555152785623E-31"), // (2929993913841559 / 6) / 38!
        new BigDecimal("-2.365022415700629934559635196369838E-32"), // (-261082718496449122051 / 13530) / 40!
        new BigDecimal("5.990671762482134304659912396819658E-34"), // (1520097643918070802691 / 1806) / 42!
        new BigDecimal("-1.517454884468290261710813135864719E-35"), // (-27833269579301024235023 / 690) / 44!
        new BigDecimal("3.843758125454188232229445290990232E-37"), // (596451111593912163277961 / 282) / 46!
        new BigDecimal("-9.736353072646691035267621279250454E-39"), // (-5609403368997817686249127547 / 46410) / 48!
        new BigDecimal("2.466247044200680957106400280288843E-40"), // (495057205241079648212477525 / 66) / 50!
        new BigDecimal("-6.247076741820743693148756794723369E-42"), // (-801165718135489957347924991853 / 1590) / 52!
        new BigDecimal("1.582403024464491429751081706828764E-43"), // (29149963634884862421418123812691 / 798) / 54!
        new BigDecimal("-4.008273685948935968530012190521983E-45"), // (-2479392929313226753685415739663229 / 870) / 56!
        new BigDecimal("1.015307585556955631163071394537876E-46"), // (84483613348880041862046775994036021 / 354) / 58!
        new BigDecimal("-2.571804158241871749924819409764455E-48"), // (-1215233140483755572040304994079820246041491 / 56786730) / 60!
        new BigDecimal("6.514456035233814931558434858641858E-50"), // (12300585434086858541953039857403386151 / 6) / 62!
        new BigDecimal("-1.650130990689652455506098780479323E-51"), // (-106783830147866529886385444979142647942017 / 510) / 64!
        new BigDecimal("4.179830628539475894850187234709407E-53"), // (1472600022126335654051619428551932342241899101 / 64722) / 66!
        new BigDecimal("-1.058763466770290877027042024279117E-54"), // (-78773130858718728141909149208474606244347001 / 30) / 68!
        new BigDecimal("2.681879191260770666140984858841510E-56"), // (1505381347333367003803076567377857208511438160235 / 4686) / 70!
        new BigDecimal("-6.793279351107421209527180299533895E-58"), // (-5827954961669944110438277244641067365282488301844260429 / 140100870) / 72!
        new BigDecimal("1.720757761668140490536349940758231E-59"), // (34152417289221168014330073731472635186688307783087 / 6) / 74!
        new BigDecimal("-4.358730329348893843400199849773161E-61"), // (-24655088825935372707687196040585199904365267828865801 / 30) / 76!
        new BigDecimal("1.104079290368466675083839597644427E-62"), // (414846365575400828295179035549542073492199375372400483487 / 3318) / 78!
        new BigDecimal("-2.796665513378134507204793753118627E-64"), // (-4603784299479457646935574969019046849794257872751288919656867 / 230010) / 80!
        new BigDecimal("7.084036501679470198509388422380349E-66"), // (1677014149185145836823154509786269900207736027570253414881613 / 498) / 82!
        new BigDecimal("-1.794407408289224066605257309336756E-67"), // (-2024576195935290360231131160111731009989917391198090877281083932477 / 3404310) / 84!
        new BigDecimal("4.545287063611096107085079124639448E-69"), // (660714619417678653573847847426261496277830686653388931761996983 / 6) / 86!
        new BigDecimal("-1.151334663198205181273002900786243E-70"), // (-1311426488674017507995511424019311843345750275572028644296919890574047 / 61410) / 88!
        new BigDecimal("2.916364771092361354703368980052971E-72"), // (1179057279021082799884123351249215083775254949669647116231545215727922535 / 272118) / 90!
        new BigDecimal("-7.387238263497337562573375394709737E-74"), // (-1295585948207537527989427828538576749659341483719435143023316326829946247 / 1410) / 92!
        new BigDecimal("1.871209311763795306225261871408075E-75"), // (1220813806579744469607301679413201203958508415202696621436215105284649447 / 6) / 94!
        new BigDecimal("-4.739828557761799405499563441210887E-77"), // (-211600449597266513097597728109824233673043954389060234150638733420050668349987259 / 4501770) / 96!
        new BigDecimal("1.200612599335450651981710022114129E-78"), // (67908260672905495624051117546403605607342195728504487509073961249992947058239 / 6) / 98!
        new BigDecimal("-3.041187241514292383037122074838896E-80"), // (-94598037819122125295227433069493721872702841533066936133385696204311395415197247711 / 33330) / 100!
        new BigDecimal("7.703417274705106272876503394597219E-82"), // (3204019410860907078243020782116241775491817197152717450679002501086861530836678158791 / 4326) / 102!
        new BigDecimal("-1.951298390909883071112323768145618E-83"), // (-319533631363830011287103352796174274671189606078272738327103470162849568365549721224053 / 1590) / 104!
        new BigDecimal("4.942696565159461474896400153888328E-85"), // (36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901 / 642) / 106!
    };

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

    /**
     * Define the test cases for each resource file for two argument functions.
     * This encapsulates the function to test, the expected maximum and RMS error, and
     * the resource file containing the data.
     */
    private enum BiTestCase implements TestError {
        // Test implementation. Uses combined data from multiple resources in order to
        // find N and M values. Any N above 5 works on this data.
        ZETA_5_15((s, a) -> HurwitzZetaTest.zeta(s, a, 5, 15), TEST_RESOURCES, 25, 3.5),
        ZETA_6_15((s, a) -> HurwitzZetaTest.zeta(s, a, 6, 15), TEST_RESOURCES, 4, 0.56),
        ZETA_7_15((s, a) -> HurwitzZetaTest.zeta(s, a, 7, 15), TEST_RESOURCES, 4, 0.56),
        ZETA_8_15((s, a) -> HurwitzZetaTest.zeta(s, a, 8, 15), TEST_RESOURCES, 4, 0.56),
        ZETA_9_15((s, a) -> HurwitzZetaTest.zeta(s, a, 9, 15), TEST_RESOURCES, 4, 0.56),
        ZETA_10_15((s, a) -> HurwitzZetaTest.zeta(s, a, 10, 15), TEST_RESOURCES, 4, 0.64),
        ZETA_11_15((s, a) -> HurwitzZetaTest.zeta(s, a, 11, 15), TEST_RESOURCES, 4.5, 0.66),
        ZETA_12_15((s, a) -> HurwitzZetaTest.zeta(s, a, 12, 15), TEST_RESOURCES, 4, 0.66),
        ZETA_S1_4_A1_8(HurwitzZeta::value, "hzeta_s1_4_a1_8.csv", 2.9, 0.65),
        ZETA_S1_4_A8_32(HurwitzZeta::value, "hzeta_s1_4_a8_32.csv", 3.6, 0.69),
        ZETA_S1_4_A32_2147483648(HurwitzZeta::value, "hzeta_s1_4_a32_2147483648.csv", 1.8, 0.5),
        ZETA_S1_4_A0_1(HurwitzZeta::value, "hzeta_s1_4_a0_1.csv", 1.9, 0.57),
        ZETA_S1_4_A0(HurwitzZeta::value, "hzeta_s1_4_a1e-16_1e-14.csv", 1.25, 0.22),
        ZETA_S4_32_A1_8(HurwitzZeta::value, "hzeta_s4_32_a1_8.csv", 3.22, 0.62),
        ZETA_S2_4_N_A1_7(HurwitzZeta::value, "hzeta_s2_4_na1_7_p0.5_p0x1p-1.csv", 0.63, 0.1),
        ZETA_S3_5_N_A1_7(HurwitzZeta::value, "hzeta_s3_5_na1_7_p0.5_p0x1p-1.csv", 3, 0.1),
        ZETA_S2_4_N_A8_33(HurwitzZeta::value, "hzeta_s2_4_na8_33_p0.5_p0x1p-1.csv", 0.75, 0.1),
        ZETA_S3_5_N_A8_33(HurwitzZeta::value, "hzeta_s3_5_na8_33_p0.5_p0x1p-1.csv", 0, 0),
        // Extended precision power is exact on the largest term
        ZETA_S2_4_N_A1_7_B30(HurwitzZeta::value, "hzeta_s2_4_na1_7_p0x1p-30.csv", 0, 0),
        ZETA_S3_5_N_A1_7_B30(HurwitzZeta::value, "hzeta_s3_5_na1_7_p0x1p-30.csv", 0, 0),
        ZETA_S2_4_N_A1_7_HALF_B30(HurwitzZeta::value, "hzeta_s2_4_na1_7_p0.5_p0x1p-30.csv", 0.63, 0.1),
        // The method suffers some cancellation here.
        // Further precision gains would require BigDecimal over double-double math.
        ZETA_S3_5_N_A1_7_HALF_B30(HurwitzZeta::value, "hzeta_s3_5_na1_7_p0.5_p0x1p-30.csv", 6, 1.7),
        ZETA_S3_5_N_A8_33_HALF_B30(HurwitzZeta::value, "hzeta_s3_5_na8_33_p0.5_p0x1p-30.csv", 180, 5.2);

//        JDK Temurin 25.492-b09
//        ZETA_5_15                             max    22.7482   RMS    3.32905   mean        1.28589  n 18000
//        ZETA_6_15                             max    3.50083   RMS   0.545590   mean    -0.00931504  n 18000
//        ZETA_7_15                             max    3.19219   RMS   0.543611   mean     -0.0108850  n 18000
//        ZETA_8_15                             max    3.52850   RMS   0.541361   mean    -0.00823326  n 18000
//        ZETA_9_15                             max    3.46324   RMS   0.546674   mean    -0.00892159  n 18000
//        ZETA_10_15                            max    3.80492   RMS   0.540870   mean    -0.00484315  n 18000
//        ZETA_11_15                            max    4.18614   RMS   0.551071   mean   -0.000419412  n 18000
//        ZETA_12_15                            max    3.77276   RMS   0.545796   mean    -0.00703945  n 18000
//        ZETA_S1_4_A1_8                        max    2.81580   RMS   0.641152   mean   -0.000526424  n 3000
//        ZETA_S1_4_A8_32                       max    3.52850   RMS   0.674928   mean     -0.0209228  n 3000
//        ZETA_S1_4_A32_2147483648              max    1.77377   RMS   0.470439   mean    -0.00772827  n 3000
//        ZETA_S1_4_A0_1                        max    1.86288   RMS   0.550082   mean    -0.00582330  n 3000
//        ZETA_S1_4_A0                          max    1.22616   RMS   0.127784   mean    -0.00596360  n 3000
//        ZETA_S4_32_A1_8                       max    3.18614   RMS   0.592956   mean    -0.00843521  n 3000
//        ZETA_S2_4_N_A1_7                      max   0.600309   RMS  0.0644069   mean   -0.000470018  n 3000
//        ZETA_S3_5_N_A1_7                      max    2.89123   RMS  0.0653419   mean   -0.000183116  n 3000
//        ZETA_S2_4_N_A8_33                     max   0.734873   RMS  0.0669002   mean    0.000287677  n 3000
//        ZETA_S3_5_N_A8_33                     max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        ZETA_S2_4_N_A1_7_B30                  max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        ZETA_S3_5_N_A1_7_B30                  max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        ZETA_S2_4_N_A1_7_HALF_B30             max   0.608751   RMS  0.0652802   mean     0.00355178  n 3000
//        ZETA_S3_5_N_A1_7_HALF_B30             max    5.66317   RMS    1.59232   mean     -0.0478892  n 3000
//        ZETA_S3_5_N_A8_33_HALF_B30            max    173.860   RMS    5.02775   mean     -0.0120123  n 3000
//        zeta  N=8   M=1   6596
//        zeta  N=8   M=2   170
//        zeta  N=8   M=3   133
//        zeta  N=8   M=4   165
//        zeta  N=8   M=5   1132
//        zeta  N=8   M=6   1837
//        zeta  N=8   M=7   1678
//        zeta  N=8   M=8   2585
//        zeta  N=8   M=9   3704


        /** The function. */
        private final DoubleBinaryOperator fun;

        /** The filenames containing the test data. */
        private final String[] filename;

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
        BiTestCase(DoubleBinaryOperator fun, String filename, double maxUlp, double rmsUlp) {
            this.fun = fun;
            this.filename = new String[] {filename};
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
        }

        /**
         * Create an instance.
         *
         * @param fun function to test
         * @param filename Filename of the test data
         * @param maxUlp maximum allowed ulp
         * @param rmsUlp maximum allowed RMS ulp
         */
        BiTestCase(DoubleBinaryOperator fun, String[] filename, double maxUlp, double rmsUlp) {
            this.fun = fun;
            this.filename = filename;
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
        }

        /**
         * @return function to test
         */
        DoubleBinaryOperator getFunction() {
            return fun;
        }

        /**
         * @return Filenames of the test data
         */
        String[] getFilenames() {
            return filename;
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
    // See STATISTICS-100 for variations tested during development.
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
    // Summing in extended precision requires a double-double (DD) sum to be used
    // throughout. Use in S and then not in the tail T does not lower the RMS.

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     * See {@link HurwitzZeta} for the formula details.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a >= 1}
     * @param n Argument {@code N}
     * @param m Argument {@code M}
     * @return zeta(s, a)
     */
    static double zeta(double s, double a, int n, int m) {
        // Skip testing
        if (n < MIN_N || n > MAX_N) {
            return Double.NaN;
        }

        // Asymptotic Behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0, the I term and the first term of T.
        if (a > 1e15) {
            return Math.pow(a, 1 - s) / (s - 1) + Math.pow(a, -s) * 0.5;
        }

        final double apn = a + n;
        double p = Math.pow(apn, -s);

        // Initialise sum with the first tail term
        double sum = 0.5 * p;
        // S : k in [0, n-1]
        for (int k = n; --k >= 0;) {
            // Descending k sums in order of magnitude for increased precision.
            // Prevents early exit for large s when the term (a+k)^-s is below
            // machine epsilon of the ascending series sum.
            sum += Math.pow(a + k, -s);
        }

        // I
        // Use of (a+n)^(1-s) = (a+n)^-1 * apn to recycle the power lowers precision.
        sum += Math.pow(apn, 1 - s) / (s - 1);

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T into the sum terms.
        // This sets the first power as (a+n)^-(1+s) not (a+n)^-1.
        // When s is large the loop exits before the rising factorial overflows.

        // Rising factorial term : (s)_{2k-1}
        double f = s;
        // 2k - 1
        double k2 = 1;
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Note: if the factor is too small (e.g. 0x1p-63) then the series continues
        // further and terms may be less accurate (i.e. add noise to the T sum).
        double tsum = 0;
        final double stop = sum * 0x1p-53;
        int i;
        for (i = 0; i < m; i++) {
            // p = (a+n)^-(2k-1+s)
            p /= apn;
            final double t = f * p / F[i];
            tsum += t;
            if (Math.abs(t) <= stop) {
                break;
            }
            p /= apn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
        }
        // Used to histogram convergence when testing
        if (n == N) {
            M[i]++;
        }
        return sum + tsum;
    }

    // TODO - Convert this to use BigDecimal

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be negative.
     *
     * @param s Argument {@code s > 1} and integer
     * @param a Argument {@code a < 0}
     * @param ca Ceil(a)
     * @return zeta(s, a)
     */
    private static double zetaNegativeImp(double s, double a, double ca) {
        // a < 0 (non-integer) and s is a positive integer.
        // If s is odd then the negative series sum will be negative
        // and the addition of zeta(s, x > 0) has cancellation.
        // This is largest when a is close to half-integer.

        // Check case of total cancellation.
        final boolean odd = SpecialMath.isOdd(s);
        double xn = a - ca;
        // Intentional float comparison
        if (odd && xn == -0.5) {
            return zetaImp(s, 1 - a);
        }

        // Compute the two terms either side of zero:
        // -1 < xn < 0 < xn + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        // One term is < 0.5: 0.5^-s overflows when s >= 1024.
        DD sum = DD.of(Double.NaN);
        if (s < 1024) {
            final int n = (int) -s;
            DD pn = DD.of(xn);
            DD pp = DD.ONE.add(xn);
            // Avoid overflow issues using scaling
            final long[] expn = {0};
            final long[] expp = {0};
            if (odd) {
                // Compute accurately in extended precision to handle cancellation.
                // The double-double result is +/- 1 ULP (105 bit precision).
                pn = DDMath.pow(pn, n, expn);
                pp = DDMath.pow(pp, n, expp);
            } else {
                // Power terms accurate to at least double precision.
                pn = pn.pow(n, expn);
                pp = pp.pow(n, expp);
            }
            // If re-scaling and addition create infinity we exit with the IEEE result.
            // Note: if one side overflows then it is unlikely the other side will
            // bring it back to finite and we do not check.
            pn = pn.scalb((int) expn[0]);
            pp = pp.scalb((int) expp[0]);
            sum = pn.add(pp);
        }
        // Check for overflow or return the IEEE result
        if (!sum.isFinite()) {
            return odd && Math.abs(xn) < xn + 1 ?
                Double.NEGATIVE_INFINITY :
                Double.POSITIVE_INFINITY;
        }

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2.
        // Adding to the existing finite sum cannot trigger overflow.

        double xp = 2 + xn;
        if (odd) {
            // If odd compute the terms in extended precision.
            // This is impractical for large |a| so the number of terms
            // is limited and precision will be lost for large |a|.
            // The number of terms depends on how close to
            // half-integer and the size of the exponent. 
            // The sum continues until the cancellation in opposing terms is in
            // the low part of the double-double sum. Computing the remaining
            // terms in double precision will have cancellation, and the difference
            // of their sums will overlap the high part of the double-double sum.
            // The degree of cancellation is dependent on the number of remaining
            // negative terms as the positive zeta evaluation to infinity will be
            // larger and more accurate. Ideally the remaining double precision sums
            // have missing bits that do not affect the result. In practice this
            // strategy works to compute a result with many bit of precision and
            // avoid a useless result with catastrophic cancellation.
            // 
            // sum          |--------|--------|
            // sp1        |--------|
            // -sn1        |------xx|
            // sp1 + sn1          |----xxxx|
            //
            // When a is close to integer this exits very fast otherwise the
            // number of terms can be large.
            // In the extreme this is limited to 5430 terms when 0.5 +/- 2^-40.
            final int n = (int) -s;
            final double x = xn;
            for (int i = 0; xn > a; i++) {
                xn -= 1.0;
                // Note: DDMath here makes no difference as s is small and the
                // standard pow function is accurate to ~100 bits. If s is large
                // then the terms rapidly reduce in magnitude compared to 0.5^-s
                // and trailing imprecise bits in the term do not change the sum.
                final DD pn = DD.of(xn).pow(n);
                final DD pp = DD.ofSum(2 + i, x).pow(n);
                final DD term = pn.add(pp);
                if (Math.abs(term.hi()) < Math.abs(sum.lo())) {
                    // Switch to a double precision tail.
                    // Reset xn to compute as part of the remaining series.
                    xn += 1.0;
                    break;
                }
                sum = sum.add(term);
            }
            // advance positive x by the number of terms computed (x - xn)
            xp = 2 + (x - xn) + x;
        }

        // Compute the remaining terms
        final double sn1 = negativeSeriesSum(a, xn, s);
        final double sp1 = zetaImp(s, xp);

        return sum.add(DD.ofSum(sp1, sn1)).hi();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 0}
     * @return zeta(s, a)
     */
    private static double zetaImp(double s, double a) {
        // Asymptotic Behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0, the I term and the first term of T.
        if (a > 1e16) {
            return Math.pow(a, 1 - s) / (s - 1) + Math.pow(a, -s) * 0.5;
        }

        final double apn = a + N;
        double p = Math.pow(apn, -s);

        // Initialise sum with the first tail term
        double sum = 0.5 * p;
        // S : k in [0, n-1]
        for (int k = N - 1; k >= 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum += Math.pow(a + k, -s);
        }

        // I
        sum += Math.pow(apn, 1 - s) / (s - 1);

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        // When s is large the loop exits before the rising factorial overflows.

        // Rising factorial term : (s)_{2k-1}
        double f = s;
        // 2k - 1
        double k2 = 1;
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        double tsum = 0;
        final double stop = sum * 0x1p-53;
        int i;
        for (i = 0; i < F.length; i++) {
            // p = (a+n)^-(2k-1+s)
            p /= apn;
            final double t = f * p / F[i];
            tsum += t;
            if (Math.abs(t) <= stop) {
                break;
            }
            p /= apn;
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
        }
        return sum + tsum;
    }

    /**
     * Calculates the sum of terms of the power series.
     *
     * <pre>
     *      b-1   1
     *   sum     ---
     *      k=a  k^m
     * </pre>
     *
     * <p>Assumes {@code a} and {@code b} are negative and separated by an integer
     * distance; and {@code exponent >= 2} and integer.
     *
     * <p>Large ranges may be evaluated using a difference of zeta functions.
     *
     * @param a First term in the series to calculate (negative non-integer).
     * @param b Last term in the series to calculate, exclusive (negative non-integer).
     * @param s Exponent (positive integer).
     * @return the sum
     */
    private static double negativeSeriesSum(double a, double b, double s) {
        // This can be computed using a difference of zeta functions.
        // A single call to zeta uses ~10 pow operations; use zeta when the
        // sum will use more.
        // Note: The difference incurs cancellation.
        // When s is even the function is called with b in -[1, 0) and the
        // series is strongly converging and no issue occurs.
        // When s is odd it may be called with large |b|. In this case many
        // terms have been computed to handle most of the cancellation in the
        // result. Worst case is b ~ 5430:
        // zeta(3, 5430) = 1.696e-08
        // zeta(3, 5450) = 1.683e-08
        // Cancellation in lost bits = exponent(max(a, b)) - exponent(a-b) = 7
        // The result is sufficient for reasonable double precision.
        if (b - a > 20) {
            final int sign = SpecialMath.isOdd(s) ? -1 : 1;
            final double zb = zetaImp(s, 1 - b);
            final double za = zetaImp(s, 1 - a);
            return sign * (zb - za);
        }

        // Sum terms in ascending order of magnitude
        double sum = 0;
        double x = a;
        while (x < b) {
            sum += Math.pow(x, -s);
            x += 1.0;
        }
        return sum;
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
        final MathContext mc = MathContext.DECIMAL128;
        // Check factors
        for (int k = 1; k < NUM.length; k++) {
            factorial = factorial.multiply(BigInteger.valueOf(2 * k - 1)).multiply(BigInteger.valueOf(2 * k));
            final BigInteger num = new BigInteger(NUM[k].substring(NUM[k].indexOf(' ') + 1));
            final BigInteger denom = new BigInteger(DENOM[k].substring(DENOM[k].indexOf(' ') + 1));

            // 2k! / B_2k
            final BigFraction factor1 = BigFraction.of(factorial.multiply(denom), num);
            final double d1 = factor1.doubleValue();
            // Cross verify BigFraction vs BigDecimal
            BigDecimal v = factor1.bigDecimalValue(mc);
            Assertions.assertEquals(v.doubleValue(), d1);
            // Find ULP precision
            final double e1 = new BigDecimal(d1).subtract(v)
                .divide(new BigDecimal(Math.ulp(d1)), mc).doubleValue();
            sum1 += Math.abs(e1);

            Assertions.assertEquals(d1, F[k - 1]);

            // Format to print the table:
            // "%s, // %s! / (%s / %s)%n", d2, 2 * k, num, denom

            // B_2k / 2k!
            final BigFraction factor2 = BigFraction.of(num, factorial.multiply(denom));
            final double d2 = factor2.doubleValue();
            // Cross verify BigFraction vs BigDecimal
            v = factor2.bigDecimalValue(mc);
            Assertions.assertEquals(v.doubleValue(), d2);
            // Find ULP precision
            final double e2 = new BigDecimal(d2).subtract(v)
                .divide(new BigDecimal(Math.ulp(d2)), mc).doubleValue();
            sum2 += Math.abs(e2);

            Assertions.assertEquals(d2, FM[k - 1]);
            Assertions.assertEquals(v, FD[k - 1]);

            // Format to print the table:
            // "%s, // (%s / %s) / %s!%n", d2, num, denom, 2 * k

            // For any M, the cumulative error is lower using 2k! / B_2k
            // as the first 6/7 factors are exact and errors in the later factors
            // are comparable.
            Assertions.assertTrue(sum1 < sum2, "2k! / B_2k does not have lower combined error");
        }
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

            // Requires M=12 when N=9.
            // This is largest M noted during development when the RMS error is close to optimal.
            Arguments.of(31.76, 8.23, 8.6811830191090714059294777379e-30, 1),
            Arguments.of(31.76, 11.23, 4.68180272588528321856530594913e-34, 0),
            Arguments.of(31.76, 12.23, 3.17258072050987981477830566397e-35, 0),
            Arguments.of(31.76, 13.23, 2.66312289027731475813712555339e-36, 1),
            Arguments.of(31.76, 14.23, 2.68377485329827353841985109724e-37, 2),
            Arguments.of(31.76, 15.23, 3.1669792006042583105605159352e-38, 1),
            Arguments.of(31.76, 17.23, 6.55526711653597399706064463689e-40, 0),
            Arguments.of(31.76, 19.23, 2.08909657211945256234757218604e-41, 0),

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
            Arguments.of(1.5, 0.1, 34.0529755150756003469433380579, 0),
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
            Arguments.of(3, -21.499, -0.09637775418460447271221850760411324846873791133448, 37),
            Arguments.of(3, -21.501, 0.098442803974649377073746600966748291488497220978961, 37),
            Arguments.of(5, -21.499, -0.64094298652061283507484770384318604647832428321319, 16),
            Arguments.of(5, -21.501, 0.64094511727200934057870930923285174211166495014722, 16),
            Arguments.of(3, -7.5000000001, 0.007782265638668063994730919817396152634504258389332, 11),
            Arguments.of(3, -7.4999999999, 0.0077822461572358593294768713462220939956242365981849, 28),
            Arguments.of(3, -7.499999999999999, 0.0077822558978654467309133842121074379966549149936996, 20),
            Arguments.of(3, -7.500000000000001, 0.0077822558980384765932871763254914638339941060968054, 18),

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
            // 0.5^-s overflows with odd or even s
            Arguments.of(1025, -1.25, -inf, 0), // -1.29e617
            Arguments.of(1024, -1.75, inf, 0), // 3.23e616
            // very large s: will be odd if cast to a long which triggers the cancellation path with s half-integer
            Arguments.of(1e+19, -1.25, inf, 0), // 1.88e+6020599913279623904
            Arguments.of(1e+19, -1.5, inf, 0), // 2.74e+3010299956639811952
            Arguments.of(1e+19, -1.75, inf, 0), // 1.88e+6020599913279623904

            // Note: large negative a can be evaluated as complex using mpmath.
            // To obtain a real result requires the digits of precision (dps)
            // to be set above |a|. Larger a are tested using Matlab.
            // mp.dps = 70
            Arguments.of(2, -66.26738, 17.78437226907544744370996767080669047664124490204306137519070727420819, 0),
            Arguments.of(3, -66.26738, -50.12250472489562851240391661869509127911964240532672215389564430564016, 0),

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
            Arguments.of(3, -0.5 + 0x1p-54, 0.41439832211715462961751618626938, 1),
            Arguments.of(3, -0.5 + 0x1p-53, 0.41439832211714926143686524195861, 1),
            Arguments.of(3, -0.5 - 0x1p-53, 0.41439832211717073415946901920171, 2),
            Arguments.of(3, -1.5 + 0x1p-52, 0.11810202582084209719727895331851, 2),
            Arguments.of(3, -1.5 - 0x1p-52, 0.11810202582088530580646271524921, 2),
            Arguments.of(3, -3.5 + 0x1p-51, 0.0307784106604705956811455131349, 4),
            Arguments.of(3, -3.5 - 0x1p-51, 0.030778410660557098867785659805987, 2),
            Arguments.of(3, -7.5 + 0x1p-50, 0.0077822558978654467309133842121074, 1),
            Arguments.of(3, -7.5 - 0x1p-50, 0.0077822558980384765932871763254915, 3),
            Arguments.of(3, -15.5 + 0x1p-49, 0.0019512219787038490075394821332913, 4),
            Arguments.of(3, -15.5 - 0x1p-49, 0.0019512219790499147520220765569762, 4),
            Arguments.of(3, -31.5 + 0x1p-48, 0.00048816210820002952659777173497611, 2),
            Arguments.of(3, -31.5 - 0x1p-48, 0.00048816210889216253017517339304625, 3),
            Arguments.of(3, -63.5 + 0x1p-47, 0.00012206286228806035641514115942213, 10),
            Arguments.of(3, -63.5 - 0x1p-47, 0.00012206286367232674283574332594035, 10),
            Arguments.of(3, -127.5 + 0x1p-46, 0.000030517111096024468697569156501944, 10),
            Arguments.of(3, -127.5 - 0x1p-46, 0.000030517113864557336393645086698446, 10),
            Arguments.of(3, -255.5 + 0x1p-45, 0.0000076293626591457113750252501943586, 10),
            Arguments.of(3, -255.5 - 0x1p-45, 0.0000076293681962114704832983463164781, 10),
            Arguments.of(3, -511.5 + 0x1p-44, 0.0000019073412767613820522958206941628, 10),
            Arguments.of(3, -511.5 - 0x1p-44, 0.0000019073523508929061980225612021874, 10),
            Arguments.of(3, -1023.5 + 0x1p-43, 0.00000047682597038482563656491683673066, 10),
            Arguments.of(3, -1023.5 - 0x1p-43, 0.00000047684811864787541032292535846052, 10),
            Arguments.of(3, -2047.5 + 0x1p-42, 0.00000011918713418230492155747450649894, 10),
            Arguments.of(3, -2047.5 - 0x1p-42, 0.00000011923143070840483965021033635372, 10),
            Arguments.of(3, -4095.5 + 0x1p-41, 0.000000029758025417506153675797113283396, 10),
            Arguments.of(3, -4095.5 - 0x1p-41, 0.000000029846618469706082505485151864752, 10),
            // Limit of a double-double summation of opposing terms from 0.5 +/- 2^-40
            Arguments.of(3, -8191.5 + 0x1p-40, 0.0000000073619875169683123404158617742165, 10),
            Arguments.of(3, -8191.5 - 0x1p-40, 0.0000000075391736213681931608483274726682, 10),
            Arguments.of(3, -16383.5 + 0x1p-39, 0.0000000016854590430963498434783046783686, 10),
            Arguments.of(3, -16383.5 - 0x1p-39, 0.0000000020398312518961172746074881289298, 10),
            // 50 ulp = Loss of ~ 6-bits of precision
            Arguments.of(3, -65535.5 + 0x1p-37, -0.00000000059232909577937793989464549320949, 50),
            Arguments.of(3, -65535.5 - 0x1p-37, 0.00000000082515973941969504164666837084463, 50),
            Arguments.of(3, -1048575.5 + 0x1p-33, -0.000000011339455934241697907112956454121, 10),
            Arguments.of(3, -1048575.5 - 0x1p-33, 0.000000011340365428943470628555741917531, 10),

            // Check some larger powers
            Arguments.of(7, -31.5 + 0x1p-48, 0.00000000014222080058151244961535769914742, 5),
            Arguments.of(9, -127.5 - 0x1p-46, 0.00000000026193893956547650755982737479515, 0),
            Arguments.of(5, -1023.5 - 0x1p-43, 0.00000000007309223831961703808738193264568, 1),

            Arguments.of(1.5, 4789, 0.0289021565574206125831622859402, 0),
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
    @EnumSource(value = BiTestCase.class)
    @Order(1)
    void testZeta(BiTestCase tc) {
        assertFunction(tc);
    }

    @Test
    @Order(2)
    void testZetaM() {
        // Check the M used to converge the tail series in the chosen implementation.
        // This depends on the epsilon used to stop the sum.
        int m = 0;
        for (int i = 0; i < M.length; i++) {
            if (M[i] != 0) {
                m = i + 1;
                // This is used for testing.
                // CHECKSTYLE: stop regex
                System.out.printf("zeta  N=%-2d  M=%-2d  %d%n", N, m, M[i]);
                // CHECKSTYLE: resume regex
            }
        }
        Assertions.assertTrue(m < 15);
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
     */
    private static void assertFunction(BiTestCase tc) {
        final TestUtils.ErrorStatistics stats = new TestUtils.ErrorStatistics();
        for (final String filename : tc.getFilenames()) {
            try (DataReader in = new DataReader(filename)) {
                while (in.next()) {
                    try {
                        final double x = in.getDouble(0);
                        final double y = in.getDouble(1);
                        final double actual = tc.getFunction().applyAsDouble(x, y);
                        // Skip results from test zeta function
                        if (Double.isNaN(actual)) {
                            return;
                        }
                        final BigDecimal expected = in.getBigDecimal(2);
                        TestUtils.assertEquals(expected, actual, tc.getTolerance(), stats::add,
                            () -> tc + " x=" + x + ", y=" + y);
                    } catch (final NumberFormatException ex) {
                        Assertions.fail("Failed to load data: " + Arrays.toString(in.getFields()), ex);
                    }
                }
            } catch (final IOException ex) {
                Assertions.fail("Failed to load data: " + filename, ex);
            }
        }

        assertRms(tc, stats);
    }

    /**
     * Assert the Root Mean Square (RMS) error of the function is below the allowed
     * maximum for the specified TestError.
     *
     * @param te Test error
     * @param stats Error statistics
     */
    private static void assertRms(TestError te, TestUtils.ErrorStatistics stats) {
        final double rms = stats.getRMS();
        debugRms(te.toString(), stats.getMaxAbs(), rms, stats.getMean(), stats.size());
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
     */
    private static void debugRms(String name, double maxAbsUlp, double rmsUlp, double meanUlp, int size) {
        // CHECKSTYLE: stop regexp
        if (!jvm) {
            jvm = true;
            System.out.printf("JDK %s %s%n",
                System.getProperty("java.vm.vendor"),
                System.getProperty("java.vm.version")
            );
        }
        System.out.printf("%-35s   max %10.6g   RMS %10.6g   mean %14.6g  n %4d%n",
            name, maxAbsUlp, rmsUlp, meanUlp, size);
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
     * Create test data for {@code s in [ls, ls + 2, ..., us - 2, us]}
     * and {@code a in -( [la, ua] + offset +/- 2^-b )}.
     *
     * <p>The parameters for {@code s} create either an even or odd series.
     *
     * <p>The parameter {@code a} is a uniform integer sample in a range. This has an
     * offset applied and a uniform wobble. This allows sampling around the poles
     * for at integer values, or the location of greatest cancellation when {@code s}
     * is odd at half-integer values.
     *
     * <p>Note: Evaluation of large negative a is not possible using the resource
     * script {@code hzeta.py}. Use of mpmath returns complex results when abs(a) >> dps
     * (dps is the mpmath digits of decimal precision).
     */
    @ParameterizedTest
    @CsvSource({
        // a = -[1.5, 7.5] +/- 0.5
        "2, 4, 1, 7, 0.5, 1",
        "3, 5, 1, 7, 0.5, 1",
        // a = -[8.5, 33.5] +/- 0.5
        "2, 4, 8, 33, 0.5, 1",
        "3, 5, 8, 33, 0.5, 1",
        // a = -[1, 7] +/- 9.31e-10
        // This approaches the pole at a = -1, -2, -3, ... and is easy to compute as
        // a single term dominates the result
        "2, 4, 1, 7, 0.0, 30",
        "3, 5, 1, 7, 0.0, 30",
        // a = -[1.5, 7.5] +/- 9.31e-10
        // This creates large cancellation when a ~ half-integer and requires
        // an extended precision power function to maintain precision over all
        // terms that cancel. The implementation only uses extended precision on
        // the most significant terms.
        "2, 4, 1, 7, 0.5, 30",
        "3, 5, 1, 7, 0.5, 30",
        // a = -[8.5, 33.5] +/- 9.31e-10
        "3, 5, 8, 33, 0.5, 30",
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
    private PrintStream getPrintStream(String filename) throws IOException {
        return getPrintStream(filename, getClass().getSimpleName());
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
