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
import java.util.function.BiFunction;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleSupplier;
import java.util.function.DoubleUnaryOperator;
import java.util.function.IntSupplier;
import java.util.stream.Stream;
import org.apache.commons.numbers.core.DD;
import org.apache.commons.numbers.core.DDMath;
import org.apache.commons.numbers.fraction.BigFraction;
import org.apache.commons.numbers.rootfinder.BrentSolver;
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
     * Used in the {@link #zeta(double, double, Context)} implementation. */
    private static final int[] M = new int[53];
    /** Filenames of resources used for the test zeta function. */
    private static final String[] TEST_RESOURCES = {
        "hzeta_s1_4_a32_2147483648.csv",
        "hzeta_s1_4_a0_1.csv",
        "hzeta_s1_4_a1e-16_1e-14.csv",
        // Higher error on these data
        "hzeta_s1_4_a1_8.csv",
        "hzeta_s1_4_a8_32.csv",
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
     * This is the same as table {@link #FM} in 70-digits of precision.
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
    };

    /**
     * Precomputed factors for {@code k}-th element of the tail function {@code T}.
     * Uses Bernoulli number {@code B_2k} divided by {@code 2k!}.
     * This is the same as table {@link #FM} in double-double precision.
     */
    private static final DD[] FDD = {
        DD.ofSum(0.08333333333333333, 4.625929269271485E-18), // (1 / 6) / 2!
        DD.ofSum(-0.001388888888888889, 5.300543954373577E-20), // (-1 / 30) / 4!
        DD.ofSum(3.306878306878307E-5, -2.2300719288557665E-21), // (1 / 42) / 6!
        DD.ofSum(-8.267195767195768E-7, 3.457597454003665E-23), // (-1 / 30) / 8!
        DD.ofSum(2.08767569878681E-8, -1.2073450591132599E-24), // (5 / 66) / 10!
        DD.ofSum(-5.284190138687493E-10, 3.517096671929869E-27), // (-691 / 2730) / 12!
        DD.ofSum(1.3382536530684679E-11, -2.828354019907999E-29), // (7 / 6) / 14!
        DD.ofSum(-3.3896802963225827E-13, -1.4986928409964295E-29), // (-3617 / 510) / 16!
        DD.ofSum(8.586062056277845E-15, -6.05252374381974E-31), // (43867 / 798) / 18!
        DD.ofSum(-2.174868698558062E-16, 4.961617782549996E-33), // (-174611 / 330) / 20!
        DD.ofSum(5.5090028283602295E-18, -1.49827152194499E-35), // (854513 / 138) / 22!
        DD.ofSum(-1.3954464685812522E-19, -1.0350590497256251E-35), // (-236364091 / 2730) / 24!
        DD.ofSum(3.534707039629467E-21, 1.894231142684204E-37), // (8553103 / 6) / 26!
        DD.ofSum(-8.953517427037546E-23, -5.728752743153026E-39), // (-23749461029 / 870) / 28!
        DD.ofSum(2.267952452337683E-24, 1.3043458462619563E-40), // (8615841276005 / 14322) / 30!
        DD.ofSum(-5.744790668872202E-26, 1.663242973708004E-43), // (-7709321041217 / 510) / 32!
        DD.ofSum(1.455172475614865E-27, -5.613265715443096E-44), // (2577687858367 / 6) / 34!
        DD.ofSum(-3.6859949406653103E-29, 1.0778256413554197E-45), // (-26315271553053477373 / 1919190) / 36!
        DD.ofSum(9.336734257095045E-31, -3.9347970210731877E-47), // (2929993913841559 / 6) / 38!
        DD.ofSum(-2.36502241570063E-32, 2.0347170931532494E-49), // (-261082718496449122051 / 13530) / 40!
        DD.ofSum(5.990671762482134E-34, 1.6265467158179092E-50), // (1520097643918070802691 / 1806) / 42!
        DD.ofSum(-1.5174548844682903E-35, 5.493014407946745E-52), // (-27833269579301024235023 / 690) / 44!
        DD.ofSum(3.843758125454189E-37, -3.685053096067968E-53), // (596451111593912163277961 / 282) / 46!
        DD.ofSum(-9.736353072646691E-39, 2.258059165188444E-55), // (-5609403368997817686249127547 / 46410) / 48!
        DD.ofSum(2.466247044200681E-40, -1.505641802268162E-56), // (495057205241079648212477525 / 66) / 50!
        DD.ofSum(-6.247076741820743E-42, -2.7106815859687654E-58), // (-801165718135489957347924991853 / 1590) / 52!
        DD.ofSum(1.5824030244644914E-43, 2.545428531496969E-60), // (29149963634884862421418123812691 / 798) / 54!
        DD.ofSum(-4.008273685948936E-45, -2.2124211668946826E-61), // (-2479392929313226753685415739663229 / 870) / 56!
        DD.ofSum(1.0153075855569557E-46, -9.404269751258486E-63), // (84483613348880041862046775994036021 / 354) / 58!
        DD.ofSum(-2.5718041582418717E-48, -6.537655454012542E-65), // (-1215233140483755572040304994079820246041491 / 56786730) / 60!
        DD.ofSum(6.514456035233815E-50, -2.763626172529861E-66), // (12300585434086858541953039857403386151 / 6) / 62!
        DD.ofSum(-1.6501309906896525E-51, 3.1794529475063687E-68), // (-106783830147866529886385444979142647942017 / 510) / 64!
        DD.ofSum(4.179830628539476E-53, 2.617556823159939E-69), // (1472600022126335654051619428551932342241899101 / 64722) / 66!
        DD.ofSum(-1.058763466770291E-54, 6.6915528436035195E-71), // (-78773130858718728141909149208474606244347001 / 30) / 68!
        DD.ofSum(2.6818791912607708E-56, -8.70695425146146E-73), // (1505381347333367003803076567377857208511438160235 / 4686) / 70!
        DD.ofSum(-6.793279351107421E-58, 2.795667911354165E-74), // (-5827954961669944110438277244641067365282488301844260429 / 140100870) / 72!
        DD.ofSum(1.7207577616681404E-59, 4.65433497191727E-76), // (34152417289221168014330073731472635186688307783087 / 6) / 74!
        DD.ofSum(-4.358730329348894E-61, 2.8840522874209336E-77), // (-24655088825935372707687196040585199904365267828865801 / 30) / 76!
        DD.ofSum(1.1040792903684666E-62, 6.624841731022409E-79), // (414846365575400828295179035549542073492199375372400483487 / 3318) / 78!
        DD.ofSum(-2.7966655133781345E-64, 2.628041826403209E-81), // (-4603784299479457646935574969019046849794257872751288919656867 / 230010) / 80!
        DD.ofSum(7.084036501679471E-66, -5.026235239023924E-82), // (1677014149185145836823154509786269900207736027570253414881613 / 498) / 82!
        DD.ofSum(-1.794407408289224E-67, 1.5372719769275798E-84), // (-2024576195935290360231131160111731009989917391198090877281083932477 / 3404310) / 84!
        DD.ofSum(4.545287063611096E-69, 9.87696151726261E-87), // (660714619417678653573847847426261496277830686653388931761996983 / 6) / 86!
        DD.ofSum(-1.1513346631982051E-70, -7.192856523313341E-87), // (-1311426488674017507995511424019311843345750275572028644296919890574047 / 61410) / 88!
        DD.ofSum(2.9163647710923614E-72, -3.8911087510195904E-89), // (1179057279021082799884123351249215083775254949669647116231545215727922535 / 272118) / 90!
        DD.ofSum(-7.387238263497337E-74, -6.923136687699924E-90), // (-1295585948207537527989427828538576749659341483719435143023316326829946247 / 1410) / 92!
        DD.ofSum(1.8712093117637953E-75, 1.5886680102062367E-92), // (1220813806579744469607301679413201203958508415202696621436215105284649447 / 6) / 94!
        DD.ofSum(-4.739828557761799E-77, -9.517121002177184E-94), // (-211600449597266513097597728109824233673043954389060234150638733420050668349987259 / 4501770) / 96!
        DD.ofSum(1.2006125993354507E-78, -1.109850335891779E-95), // (67908260672905495624051117546403605607342195728504487509073961249992947058239 / 6) / 98!
        DD.ofSum(-3.0411872415142924E-80, 5.125117133572647E-97), // (-94598037819122125295227433069493721872702841533066936133385696204311395415197247711 / 33330) / 100!
        DD.ofSum(7.703417274705106E-82, 3.948211996024456E-99), // (3204019410860907078243020782116241775491817197152717450679002501086861530836678158791 / 4326) / 102!
        DD.ofSum(-1.951298390909883E-83, -1.2533409416284754E-99), // (-319533631363830011287103352796174274671189606078272738327103470162849568365549721224053 / 1590) / 104!
        DD.ofSum(4.942696565159462E-85, -2.8094990080509668E-101), // (36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901 / 642) / 106!
    };

    // TODO - refactor Context to double, double-double and BigDecimal implementations

    /** Context for the zeta implementation.
     * This is used to control how the test zeta implementation is computed.
     * Contains options for double, double-double and BigDecimal implementations. */
    private static class Context {
        // Default contexts optimised using the precision tests

        /** Context for double precision zeta implementation.
         * N=9 has max M used as 9 in the test data. */
        public static final Context DOUBLE = Context.of(9, 15);
        // TODO: To be optimised
        /** Context for double-double precision using the DD zeta implementation. */
        public static final Context DOUBLE_DOUBLE = Context.of(20, 30).withEpsilon(0x1p-106);
        /** Context for double-double precision using the BigDecimal zeta implementation. */
        public static final Context BD = Context.of(20, 30).withMathContext(new MathContext(37));
        /** Context for double precision using the BigDecimal zeta implementation. */
        public static final Context BD_DOUBLE = Context.of(9, 15).withMathContext(new MathContext(20));

        /** Default epsilon. */
        private static final double EPS = 0x1p-53;
        /** Default math context. */
        private static final MathContext MC = MathContext.DECIMAL128;
        /** Default pow(n+x, y) function. */
        private static final boolean POWNP = true;
        /** Default option for tail series. */
        private static final int TAIL = 0;
        /** Default option for extended precision sum. */
        private static final boolean EP_SUM = true;

        /** N. */
        private final int n;
        /** M. */
        private final int m;
        /** Epsion for convergence of the tail series. */
        private final double eps;
        /** Math context for extended precision evaluations. */
        private final MathContext mc;
        /** Use extended precision Math.pow(n + x, y). */
        private final boolean usePowNp;
        /** Tail series option. */
        private final int tailOption;
        /** Use extended precision sum. */
        private final boolean useExtendedPrecisionSum;

        /**
         * Create an instance.
         *
         * @param n the n
         * @param m the m
         * @param eps the eps
         * @param mc the mc
         * @param usePowNp flag for extended precision pow(n+x, y)
         * @param tailOption the tail option
         * @param useExtendedPrecisionSum Use an extended precision sum
         */
        Context(int n, int m, double eps, MathContext mc,
            boolean usePowNp, int tailOption, boolean useExtendedPrecisionSum) {
            if (m >= F.length) {
                throw new IllegalArgumentException("Unsupported M: " + m);
            }
            this.n = n;
            this.m = m;
            this.eps = eps;
            this.mc = mc;
            this.usePowNp = usePowNp;
            this.tailOption = tailOption;
            this.useExtendedPrecisionSum = useExtendedPrecisionSum;
        }

        /**
         * Create a context.
         *
         * @param n the number of terms N.
         * @param m the number of terms M.
         * @return the context
         */
        static Context of(int n, int m) {
            return new Context(n, m, EPS, MC, POWNP, TAIL, EP_SUM);
        }

        /**
         * Return a context that uses the provided convergence epsilon in the tail series.
         *
         * @param value the value
         * @return the context
         */
        Context withEpsilon(double value) {
            return new Context(n, m, value, mc, usePowNp, tailOption, useExtendedPrecisionSum);
        }

        /**
         * Return a context that uses the provided MathContext.
         *
         * @param value the value
         * @return the context
         */
        Context withMathContext(MathContext value) {
            return new Context(n, m, eps, value, usePowNp, tailOption, useExtendedPrecisionSum);
        }

        /**
         * Return a context that uses extended precision for Math.pow(n + x, y).
         *
         * @param value the value
         * @return the context
         */
        Context withUseNp(boolean value) {
            return new Context(n, m, eps, mc, value, tailOption, useExtendedPrecisionSum);
        }

        /**
         * Return a context that uses the provided option in the tail series.
         *
         * @param value the value
         * @return the context
         */
        Context withTailOption(int value) {
            return new Context(n, m, eps, mc, usePowNp, value, useExtendedPrecisionSum);
        }

        /**
         * Return a context that uses an extended precision sum.
         *
         * @param value the value
         * @return the context
         */
        Context withExtendedPrecisionSum(boolean value) {
            return new Context(n, m, eps, mc, usePowNp, tailOption, value);
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
         * Gets the convergence epsilon for the tail series T.
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
         * Get the function to compute {@code (x+y)^z}.
         * @return function
         */
        DoubleTernaryOperator getPowNp() {
            return usePowNp ? Context::powNp : (x, y, z) -> Math.pow(x + y, z);
        }

        /**
         * Gets the option to use in the tail series.
         * Implementations may use this as a bit flag to change multiple options.
         * @return the option
         */
        int getTailOption() {
            return tailOption;
        }

        /**
         * Gets the option to use an extended precision sum.
         * @return the option
         */
        boolean getUseExtendedPrecisionSum() {
            return useExtendedPrecisionSum;
        }
        /**
         * Extended precision {@code (x+y)^z}.
         *
         * <p>Warning: assumes x+y is finite.
         * This method is used for testing where the arguments should be {@code n + y}
         * with y finite and less than 2^53.
         *
         * @param x the x
         * @param y the y
         * @param z the z
         * @return the result
         */
        static double powNp(double x, double y, double z) {
            // (s+ss)^y = s^y * (1+ss/s)^y
            //          = s^y * exp(y*log1p(ss/s))
            // ss/s < machine epsilon : log1p(ss/s) ~ ss/s
            final DD s = DD.ofSum(x, y);
            double r = Math.pow(s.hi(), z);
            // This does not check all pow edge cases and assumes the round-off is finite
            if (s.lo() != 0) {
                r *= Math.exp(z * s.lo() / s.hi());
            }
            return r;
        }
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

    /** Define a test for a double-precision zeta function. */
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
     * Define the test cases for each resource file for two argument functions.
     * This encapsulates the function to test, the expected maximum and RMS error, and
     * the resource file containing the data.
     */
    private enum ZetaTestCase implements DoublePrecisionTestCase {
        ZETA_S1_4_A1_8(HurwitzZeta::value, "hzeta_s1_4_a1_8.csv", 2.9, 0.65),
        ZETA_S1_4_A8_32(HurwitzZeta::value, "hzeta_s1_4_a8_32.csv", 3.6, 0.69),
        ZETA_S1_4_A32_2147483648(HurwitzZeta::value, "hzeta_s1_4_a32_2147483648.csv", 1.8, 0.5),
        ZETA_S1_4_A0_1(HurwitzZeta::value, "hzeta_s1_4_a0_1.csv", 1.9, 0.57),
        ZETA_S1_4_A0(HurwitzZeta::value, "hzeta_s1_4_a1e-16_1e-14.csv", 1.25, 0.22),
        ZETA_S4_32_A1_8(HurwitzZeta::value, "hzeta_s4_32_a1_8.csv", 3.22, 0.62),
        ZETA_S2_4_N_A1_7(HurwitzZeta::value, "hzeta_s2_4_na1_7_p0.5_p0x1p-1.csv", 0.7, 0.1),
        ZETA_S3_5_N_A1_7(HurwitzZeta::value, "hzeta_s3_5_na1_7_p0.5_p0x1p-1.csv", 3, 0.1),
        ZETA_S2_4_N_A8_33(HurwitzZeta::value, "hzeta_s2_4_na8_33_p0.5_p0x1p-1.csv", 0.75, 0.1),
        ZETA_S3_5_N_A8_33(HurwitzZeta::value, "hzeta_s3_5_na8_33_p0.5_p0x1p-1.csv", 0, 0),
        // Extended precision power is exact on the largest term
        ZETA_S2_4_N_A1_7_B30(HurwitzZeta::value, "hzeta_s2_4_na1_7_p0x1p-30.csv", 0, 0),
        ZETA_S3_5_N_A1_7_B30(HurwitzZeta::value, "hzeta_s3_5_na1_7_p0x1p-30.csv", 0, 0),
        ZETA_S2_4_N_A1_7_HALF_B30(HurwitzZeta::value, "hzeta_s2_4_na1_7_p0.5_p0x1p-30.csv", 0.63, 0.1),
        // The method suffers some cancellation here.
        // Further precision gains would require BigDecimal over double-double math.
        ZETA_S3_5_N_A1_7_HALF_B30(HurwitzZeta::value, "hzeta_s3_5_na1_7_p0.5_p0x1p-30.csv", 6.5, 1.7),
        ZETA_S3_5_N_A8_33_HALF_B30(HurwitzZeta::value, "hzeta_s3_5_na8_33_p0.5_p0x1p-30.csv", 180, 5.2),
        BD_ZETA_S3_5_N_A1_7_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_s3_5_na8_33_p0.5_p0x1p-30.csv", 0, 0),
        BD_ZETA_S3_5_N_A8_33_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_s3_5_na8_33_p0.5_p0x1p-30.csv", 0, 0),
        BD_ZETA_ROOT_S3_9_N_A0_100((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_root_s3_9_na0_100.csv", 0, 0),
        BD_ZETA_ROOT_S3_9_N_A1000_1100((s, a) -> HurwitzZetaTest.zetaNegativeBD((int) s, a), "hzeta_root_s3_9_na1000_1100.csv", 0, 0),
        DD_ZETA_S3_5_N_A1_7_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_s3_5_na8_33_p0.5_p0x1p-30.csv", 0, 0),
        DD_ZETA_S3_5_N_A8_33_HALF_B30((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_s3_5_na8_33_p0.5_p0x1p-30.csv", 0, 0),
        // TODO - Can this be fixed?
        DD_ZETA_ROOT_S3_9_N_A0_100((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_root_s3_9_na0_100.csv", 1, 0.05),
        DD_ZETA_ROOT_S3_9_N_A1000_1100((s, a) -> HurwitzZetaTest.zetaNegativeDD((int) s, a), "hzeta_root_s3_9_na1000_1100.csv", 1, 0.05),
        // TODO - Must be optimised for double-double precision
        BD_ZETA_IS2_2_A40_41((s, a) -> HurwitzZetaTest.zeta((int) s, new BigDecimal(a), null, Context.BD).doubleValue(), "hzeta_ia2_2_a1_1.csv", 0, 0),
        DD_ZETA_IS2_2_A40_41((s, a) -> HurwitzZetaTest.zeta((int) s, DD.of(a), null, Context.DOUBLE_DOUBLE).doubleValue(), "hzeta_ia2_2_a1_1.csv", 0, 0),
        ;

//        JDK Temurin 25.492-b09
//        ZETA_5_15                             max    22.7482   RMS    3.32271   mean        1.28523  n 18000
//        ZETA_6_15                             max    4.03990   RMS   0.517851   mean    -0.00641581  n 18000
//        ZETA_7_15                             max    3.19219   RMS   0.515154   mean    -0.00790843  n 18000
//        ZETA_8_15                             max    3.52850   RMS   0.512206   mean    -0.00776736  n 18000
//        ZETA_9_15                             max    3.03990   RMS   0.523973   mean    -0.00502301  n 18000
//        ZETA_10_15                            max    3.80492   RMS   0.515456   mean    -0.00676918  n 18000
//        ZETA_11_15                            max    3.50083   RMS   0.525740   mean    -0.00333901  n 18000
//        ZETA_12_15                            max    3.16162   RMS   0.521647   mean    -0.00653556  n 18000
//        ZETA_S1_4_A1_8                        max    2.81580   RMS   0.609868   mean     0.00407202  n 3000
//        ZETA_S1_4_A8_32                       max    3.52850   RMS   0.675110   mean     -0.0240288  n 3000
//        ZETA_S1_4_A32_2147483648              max    1.77377   RMS   0.470439   mean    -0.00772827  n 3000
//        ZETA_S1_4_A0_1                        max    1.83029   RMS   0.483327   mean     -0.0128829  n 3000
//        ZETA_S1_4_A0                          max   0.952404   RMS   0.107769   mean    -0.00278978  n 3000
//        ZETA_S4_32_A1_8                       max    3.05405   RMS   0.529045   mean    -0.00324641  n 3000
//        ZETA_S2_4_N_A1_7                      max   0.653535   RMS  0.0672688   mean    0.000652193  n 3000
//        ZETA_S3_5_N_A1_7                      max    2.89123   RMS  0.0658265   mean   -0.000869074  n 3000
//        ZETA_S2_4_N_A8_33                     max   0.633617   RMS  0.0622978   mean    0.000343235  n 3000
//        ZETA_S3_5_N_A8_33                     max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        ZETA_S2_4_N_A1_7_B30                  max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        ZETA_S3_5_N_A1_7_B30                  max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        ZETA_S2_4_N_A1_7_HALF_B30             max   0.589429   RMS  0.0613828   mean     0.00299001  n 3000
//        ZETA_S3_5_N_A1_7_HALF_B30             max    6.26561   RMS    1.58252   mean     -0.0695818  n 3000
//        ZETA_S3_5_N_A8_33_HALF_B30            max    173.860   RMS    5.02735   mean     -0.0130644  n 3000
//        BD_ZETA_S3_5_N_A1_7_HALF_B30          max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        BD_ZETA_S3_5_N_A8_33_HALF_B30         max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        BD_ZETA_ROOT_S3_9_N_A0_100            max    0.00000   RMS    0.00000   mean        0.00000  n 1212
//        DD_ZETA_S3_5_N_A1_7_HALF_B30          max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        DD_ZETA_S3_5_N_A8_33_HALF_B30         max    0.00000   RMS    0.00000   mean        0.00000  n 3000
//        DD_ZETA_ROOT_S3_9_N_A0_100            max   0.638297   RMS  0.0362792   mean   -0.000613956  n 1212
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
        ZetaTestCase(DoubleBinaryOperator fun, String filename, double maxUlp, double rmsUlp) {
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
        ZetaTestCase(DoubleBinaryOperator fun, String[] filename, double maxUlp, double rmsUlp) {
            this.fun = fun;
            this.filename = filename;
            this.maxUlp = maxUlp;
            this.rmsUlp = rmsUlp;
        }

        @Override
        public
        DoubleBinaryOperator getFunction() {
            return fun;
        }

        @Override
        public
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
    //
    // Note: The method is sensitive to the initial loop over N to create S.
    // Under certain conditions the N cannot be too high if using an ascending
    // sum of k as the sum does not converge and later terms are added with
    // low precision.
    //
    // Better results are obtained using descending k. However this prevents
    // an early exit if the series is rapidly converging and the term (a+k)^-s
    // drops below machine epsilon of the sum.

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
        final int n = c.getN();

        // Asymptotic Behavior as a -> inf
        // https://dlmf.nist.gov/25.11#E43
        // When a is large the series cannot use a+k.
        // This reduces to N=0, the I term and the first term of T.
        if (a > 1e15) {
            return Math.pow(a, 1 - s) / (s - 1) + Math.pow(a, -s) * 0.5;
        }

        // Can overflow if 0 < a < 1.
        final double t0 = Math.pow(a, -s);
        if (!Double.isFinite(t0)) {
            return t0;
        }
        // Now any (a+n)^-s cannot overflow and the sum cannot overflow.

        DoubleTernaryOperator pow = c.getPowNp();

        // Check the extra precision power will make a difference
        if (!(a > 1 && Math.abs(s * DD.ofSum(a, n).lo()) >= 0x1p-53)) {
            pow = (x, y, z) -> Math.pow(x + y, z);
        }

        double p = pow.applyAsDouble(a, n, -s);

        // We always use a DD sum. If not using extended precision
        // we add in double precision and create a new DD.
        BiFunction<DD, Double, DD> add = c.getUseExtendedPrecisionSum() ?
            DD::add :
            (x, y) -> DD.of(x.hi() + y);

        // Initialise sum with the first tail term
        DD sum = DD.of(0.5 * p);
        // S : k in [0, n-1]
        for (int k = n; --k > 0;) {
            // Descending k sums in order of magnitude for increased precision.
            // Prevents early exit for large s when the term (a+k)^-s is below
            // machine epsilon of the ascending series sum.
            sum = add.apply(sum, pow.applyAsDouble(a, k, -s));
        }

        // I : (a+p)^(1-s) / (s-1)
        // Use of (a+n)^(1-s) = (a+n)^-1 * apn to recycle the power lowers precision.
        final double ti = pow.applyAsDouble(a, n, 1 - s) / (s - 1);

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (t0 > ti) {
            sum = add.apply(sum, ti);
            sum = add.apply(sum, t0);
        } else {
            sum = add.apply(sum, t0);
            sum = add.apply(sum, ti);
        }

        // T
        // The following incorporates the factor for T into the sum terms
        // as (a+n)^-(2k-1+s).
        // This sets the first power as (a+n)^-(1+s) not (a+n)^-1.
        // When s is large the loop exits before the rising factorial overflows.
        // This factor can be computed using alternative implementations.

        // Rising factorial term : (s)_{2k-1}
        double f = s;
        // 2k - 1
        double k2 = 1;
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        // Note: if the factor is too small (e.g. 0x1p-63) then the series continues
        // further and terms may be less accurate (i.e. add noise to the T sum).

        // tsum can use extended precision.
        DD tsum = DD.ZERO;
        add = (c.getTailOption() & 4) != 0 ?
            DD::add :
            (x, y) -> DD.of(x.hi() + y);

        final double stop = sum.hi() * c.getEps();

        // Alternative implementations for (a+n)^-(2k-1+s).
        // Set using the first two bits of the tail option.
        double apn = 0;
        int powerTermOption = c.getTailOption() & 0x3;
        if (powerTermOption == 0) {
            // Initialise (a+n)^-(2k-1+s) to (a+n)^-(1+s)
            // Divide by (a+n)^2 using multiplication
            p = pow.applyAsDouble(a, n, -s - 1);
            apn = pow.applyAsDouble(a, n, -2);
        } else if (powerTermOption == 1) {
            // Initialise (a+n)^-(2k-1+s) to (a+n)^-s
            // Divide by (a+n) twice inside the loop
            apn = a + n;
        } else {
            // Compute using the power function
            p = pow.applyAsDouble(a, n, -(k2 + s));
        }

        int i;
        for (i = 0; i < c.getM(); i++) {
            // p = (a+n)^-(2k-1+s)
            if (powerTermOption == 1) {
                p /= apn;
            }
            // Note that this uses divide by F rather than multiply by FM.
            // Testing shows negligible difference. The first 6/7 terms of
            // M are exact so divide is used.
            final double t = f * p / F[i];
            tsum = add.apply(tsum, t);
            if (Math.abs(t) <= stop) {
                break;
            }
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            f *= s + k2;
            k2 += 1.0;
            f *= s + k2;
            k2 += 1.0;
            // Update (a+n)^-(2k-1+s)
            if (powerTermOption == 0) {
                p *= apn;
            } else if (powerTermOption == 1) {
                p /= apn;
            } else {
                p = pow.applyAsDouble(a, n, -(k2 + s));
            }
        }

        // Used to histogram convergence when testing
        M[i]++;
        return sum.add(tsum).hi();
    }

    // TODO
    // Is the current double-double implementation good enough?
    // How slow is the DD vs BigDecimal (optimised) for odd s.
    // The BigDecimal can compute the zeta difference accurately.
    // Simple benchmark of random values with mix of odd and even s
    // calling each implementation.

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
        double xn = a - ca;
        // Intentional float comparison
        if (odd && xn == -0.5) {
            // Use extended precision but evaluated with precision for a double result
            return zeta(s, BigDecimal.ONE.subtract(new BigDecimal(a)), null,
                Context.DOUBLE.withMathContext(new MathContext(20))).doubleValue();
        }

        // Compute dominant term using closest to zero.
        double d = Math.pow(xn > -0.5 ? xn : 1 + xn, -s);
        if (!Double.isFinite(d)) {
            return d;
        }

        // Odd computation requires twice the precision of a double (17 digits).
        // Even requires some extra.
        final Context c = odd ? Context.BD : Context.BD_DOUBLE;

        // Compute the two terms either side of zero:
        // -1 < xn < 0 < xn + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        BigDecimal x = new BigDecimal(xn);
        BigDecimal xp = BigDecimal.ONE.add(x);
        MathContext mc = c.getMathContext();
        final BigDecimal pn = x.pow(-s, mc);
        final BigDecimal pp = xp.pow(-s, mc);

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2.
        // Exit early if remaining terms cannot be added.
        d = pn.add(pp, mc).doubleValue();
        if (Math.abs(d) > 0x1p54) {
            return d;
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) +/- [ zeta(s, -x) - zeta(s, 1 - a) ]

        // Compute the remaining terms passing in the known values:
        final BigDecimal sn1 = negativeSeriesSum(a, xn, s, pn, c);
        final BigDecimal sp1 = zeta(s, xp, pp, c);

        return sp1.add(sn1, mc).doubleValue();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}
     * @param a Argument {@code a > 1}
     * @param a0 {@code a^-s}
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    private static BigDecimal zeta(int s, BigDecimal a, BigDecimal a0, Context c) {
        final int n = c.getN();
        final MathContext mc = c.getMathContext();
        final BigDecimal apn = a.add(BigDecimal.valueOf(n));
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
        // When s is large the loop exits before the rising factorial overflows.

        // Rising factorial term : (s)_{2k-1}
        BigDecimal f = BigDecimal.valueOf(s);
        p = p.divide(apn, mc);
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        BigDecimal tsum = BigDecimal.ZERO;
        // Used to divide by (a+n)^2
        BigDecimal apn2 = apn.pow(-2, mc);
        final int stop = sum.scale() + mc.getPrecision();
        int i;
        for (i = 0; i < c.getM(); i++) {
            final BigDecimal t = f.multiply(p, mc).multiply(FBD[i], mc);
            tsum = tsum.add(t, mc);
            if (t.scale() >= stop) {
//                System.out.printf("%d %s  %d%n", s, a.doubleValue(), i + 1);
                break;
            }
            // p = (a+n)^-(2k-1+s)
            p = p.multiply(apn2, mc);
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            f = f.multiply(BigDecimal.valueOf((s + (2L * i) + 1) * (s + (2L * i) + 2)), mc);
        }
        return sum.add(tsum, mc);
    }

    /**
     * Calculates the sum of terms of the power series.
     *
     * <pre>
     *      b     1
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
     * @param b Last term in the series inclusive (negative non-integer); result provided.
     * @param s Exponent (positive integer).
     * @param bn {@code b^-s}.
     * @param c Evaluation context.
     * @return the sum
     */
    private static BigDecimal negativeSeriesSum(double a, double b, int s,
            BigDecimal bn, Context c) {
        // This can be computed using a difference of zeta functions.
        // A single call to zeta uses many pow operations; use zeta when the
        // sum will use more.
        final MathContext mc = c.getMathContext();
        if (b - a > 2 * (c.getN() + 2)) {
            // Note: The difference incurs cancellation.
            // This should not be an issue as function is called with b in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 0.5)  = 1.6449
            // zeta(2, 31.5) = 0.03225
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            final BigDecimal zb = zeta(s, new BigDecimal(-b), bn.abs(), c);
            final BigDecimal za = zeta(s, BigDecimal.ONE.subtract(new BigDecimal(a)), null, c);
            final BigDecimal r = zb.subtract(za, mc);
            return (s & 1) == 1 ? r.negate() : r;
        }

        // Sum terms in ascending order of magnitude.
        // Use a double to track the iterations, and mirror with a BigDecimal.
        BigDecimal sum = BigDecimal.ZERO;
        double x = a;
        BigDecimal bx = new BigDecimal(a);
        while (x < b) {
            sum = sum.add(bx.pow(-s, mc), mc);
            x += 1.0;
            bx = bx.add(BigDecimal.ONE);
        }
        return sum.add(bn, mc);
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}
     * when {@code a} is negative. Uses double-double ({@link DD}) arithmetic.
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
        double x = a - ca;
        // Intentional float comparison
        if (odd && x == -0.5) {
            // Use extended precision but evaluated with precision for a double result
            return zeta(s, DD.ONE.subtract(a), null, Context.DOUBLE).doubleValue();
        }

        // Compute dominant term using closest to zero.
        // Note 0.5^-1024 overflows.
        double d = Math.pow(x > -0.5 ? x : 1 + x, -s);
        if (!Double.isFinite(d)) {
            return d;
        }

        // Compute the two terms either side of zero:
        // -1 < xn < 0 < xn + 1 < 1
        // These are the largest terms and contain most of the error of the function.
        DD xn = DD.of(x);
        DD xp = DD.ONE.add(x);
        final long[] expn = {0};
        final long[] expp = {0};
        DD pn = DDMath.pow(xn, -s, expn);
        DD pp = DDMath.pow(xp, -s, expp);
        pn = pn.scalb((int) expn[0]);
        pp = pp.scalb((int) expp[0]);

        // Here the remaining series above and below zero are effectively both zeta
        // evaluations with zeta(s >= 2, a > 1). This is always < 2.
        // Exit early if remaining terms cannot be added.
        d = pn.add(pp).doubleValue();
        if (Math.abs(d) > 0x1p54) {
            return d;
        }

        // x = a - ceil(a) : x in -(1, 0)
        // zeta(s, x + 1) +/- [ zeta(s, -x) - zeta(s, 1 - a) ]

        // Worst case single term cancellation:
        // pow(0.5 + 2^-53, -3) - pow(0.5 - 2^-53, -3)
        // priority = exponent(max(a, b)) - exponent(a - b) = 50 bits

        // We have the two largest power terms for each side.
        // The remaining terms are increasing smaller. Computing with
        // double-double (DD) precision for the zeta evaluations should handle cancellation.

        // Evaluate zeta with extra precision
        final Context c = Context.DOUBLE_DOUBLE;

        // Compute the remaining terms passing in the known values:
        final DD sn1 = negativeSeriesSum(a, x, s, pn, c);
        final DD sp1 = zeta(s, xp, pp, c);

        return sp1.add(sn1).hi();
    }

    /**
     * Compute the value of the Hurwitz zeta function {@code zeta(s, a)}.
     *
     * <p><strong>Warning</strong>: No parameter validation is performed.
     * The domain of {@code a} is expected to be positive.
     *
     * @param s Argument {@code s > 1}, expected {@code s < 1024}
     * @param a Argument {@code a > 1}
     * @param a0 {@code a^-s}
     * @param c Evaluation context.
     * @return zeta(s, a)
     */
    private static DD zeta(int s, DD a, DD a0, Context c) {
        // Note:
        // Closest sum of terms when cancellation occurs and we need full DD accuracy:
        // a in [0, 1]
        // (1-2^-53)^-3 + 2^-3 + 3^-3 + 4^-3 + ... ~ 1.0 + 0.125 + 0.0370 + 0.01562
        // zeta(3, 2) 0.202056
        // The initial series of terms are > 2-fold smaller. Computing with the standard
        // DD pow function has enough accuracy to not accumulate error to the zeta result.
        // When a >> 1 then the terms are all similar magnitude and we use DDmath.

        // TODO : get pow function from the context to allow testing

        final int n = c.getN();
        final DD apn = a.add(n);
        DD p = apn.pow(-s);

        // Initialise sum with the first tail term: 0.5 * (a+n)^-s
        DD sum = p.scalb(-1);
        // S : k in [0, n-1]
        for (int k = n - 1; k > 0; k--) {
            // Descending k sums in order of magnitude for increased precision
            sum = sum.add(a.add(k).pow(-s));
        }
        // First term may be provided
        final DD t0 = a0 == null ? a.pow(-s) : a0;

        // I : (a+p)^(1-s) / (s-1)
        final DD ti = apn.pow(1 - s).divide(s - 1);

        // Add in magnitude order. When a in [0, 1] it may be the dominant term
        if (t0.hi() > ti.hi()) {
            sum = sum.add(ti).add(t0);
        } else {
            sum = sum.add(t0).add(ti);
        }

        // T
        // The following recycles the power term p: (a+n)^-(2k-1+s).
        // This incorporates the factor for T, (a+n)^-s, into the sum terms.
        // The first power is (a+n)^-(1+s) not (a+n)^-1.
        // When s is large the loop exits before the rising factorial overflows.
        // Max expected s is < 1024. This overflows at k=51:
        // pochammer(1024, 101) = 1.30e+306
        // pochammer(1024, 102) = 1.46e+309

        // Rising factorial term : (s)_{2k-1}
        DD f = DD.of(s);
        p = p.divide(apn);
        // Sum of an alternating series as each F changes sign.
        // Sum until terms will not impact the result.
        DD tsum = DD.ZERO;
        // Used to divide by (a+n)^2
        DD apn2 = apn.pow(-2);
        final double stop = sum.hi() * c.getEps();
        int i;
        for (i = 0; i < c.getM(); i++) {
            final DD t = f.multiply(p).multiply(FDD[i]);
            tsum = tsum.add(t);
            if (Math.abs(t.hi()) <= stop) {
                break;
            }
            // p = (a+n)^-(2k-1+s)
            p = p.multiply(apn2);
            // f = s * (s+1) * (s+2) * ... * (s+2k-2)
            // compute the multiplicand as a long as it cannot overflow when M is small
            f = f.multiply((s + (2L * i) + 1) * (s + (2L * i) + 2));
        }
        return sum.add(tsum);
    }

    /**
     * Calculates the sum of terms of the power series.
     *
     * <pre>
     *      b     1
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
     * @param b Last term in the series inclusive (negative non-integer); result provided.
     * @param s Exponent (positive integer).
     * @param bn {@code b^-s}.
     * @param c Evaluation context.
     * @return the sum
     */
    private static DD negativeSeriesSum(double a, double b, int s,
            DD bn, Context c) {
        // This can be computed using a difference of zeta functions.
        // A single call to zeta uses many pow operations; use zeta when the
        // sum will use more.
        if (b - a > 2 * (c.getN() + 2)) {
            // Note: The difference incurs cancellation.
            // This should not be an issue as function is called with b in -(1, 0)
            // and the series is strongly converging, e.g.
            // zeta(2, 0.5)  = 1.6449
            // zeta(2, 31.5) = 0.03225
            // Significant cancellation (leading digits the same) is not possible.
            // Take care to change the sign of a provided result for the zeta method.
            final DD zb = zeta(s, DD.of(-b), bn.abs(), c);
            final DD za = zeta(s, DD.ONE.subtract(a), null, c);
            final DD r = zb.subtract(za);
            return (s & 1) == 1 ? r.negate() : r;
        }

        // Sum terms in ascending order of magnitude
        DD sum = DD.ZERO;
        double x = a;
        while (x < b) {
            sum = sum.add(DD.of(x).pow(-s));
            x += 1.0;
        }
        return sum.add(bn);
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

    /**
     * Test the double-precision zeta test function.
     */
    @ParameterizedTest
    @CsvSource({
        // Notes:
        // - N >= 7 has similar error. A lower N uses more M in the tail series.
        // - Using a high precision pow(n + x, y) lowers max error by ~1 ULP and rms by 20%.
        //   Investigation shows it makes most difference when s is large which makes
        //   sense given the implementation (a+n)^-s using -s * round-off of (a+n)
        //   (verified using test resources with s above 4).
        //   It does not impact a in [0, 1] where (a+n) is inexact as a^-s is magnitudes
        //   larger than all other terms.
        //   It could be applied when: a>1; s>2; and (a+n) is inexact. Here the round-off
        //   multiplied by s will be above 2^-53 and exp(-s * round-off) != 1.
        // - Using divide makes negligible difference from using multiply in the tail.
        // - Using the power function to compute (a+n)^-(2k-1+s) makes no difference
        //   and will be expensive.
        // - Extended precision sum makes a small difference. Max error is similar and RMS drops.
        //   Terms are added in ascending magnitude which is fine in double precision
        //   unless the largest terms have error. Can be helped using high precision pow.
        // - Extended precision tail sum makes no difference. The error is in the other part
        //   of the computation.
        // - Using a smaller epsilon requires more terms in the tail but does not impact
        //   error unless epsilon is too high.
        // - RMS drops with increasing N then plateaus. The plateau is at larger N when
        //   using the higher precision power **and** sum, e.g. N=8 vs N=10. Using one or
        //   the other N=8 is OK.

        // Extended precision power function (difference)
        "5, 15, 15, -53, false, 0, false",
        "5, 15, 15, -53, true, 0, false",
        // Extended precision sum (small difference with/without the power function)
        // RMS drops as N increases. Max error is variable.
        "5, 15, 15, -53, false, 0, true",
        "5, 15, 15, -53, true, 0, true",
//        // Using divide in the tail series (no difference)
//        "5, 12, 15, -53, false, 1, false",
//        "5, 12, 15, -53, true, 1, false",
//        // Use power in the tail series (no difference)
//        "5, 12, 15, -53, false, 2, true",
//        "5, 12, 15, -53, true, 2, true",
//        // Use extended precision sum in the tail series (no difference)
//        "5, 12, 15, -53, false, 4, true",
//        "5, 12, 15, -53, true, 4, true",
//        "5, 12, 15, -53, false, 6, true",
//        "5, 12, 15, -53, true, 6, true",
//        // Convergence (negligible error change unless to high, does increase required M)
//        "8, 10, 15, -49, true, 0, true",
//        "8, 10, 15, -50, true, 0, true",
//        "8, 10, 15, -51, true, 0, true",
//        "8, 10, 15, -52, true, 0, true",
//        "8, 10, 15, -53, true, 0, true",
//        "8, 10, 15, -54, true, 0, true",
    })
    @Disabled("Used to parameterize the zeta function")
    void testDoublePrecisionZeta(int ln, int un, int m, int b,
        boolean powNp, int tail, boolean epSum)
        throws IOException {
        final double eps = Math.scalb(1.0, b);
        for (int n = ln; n <= un; n++) {
            // Reset M
            Arrays.fill(M, 0);

            final Context c = Context.of(n, m)
                .withEpsilon(eps)
                .withUseNp(powNp)
                .withTailOption(tail)
                .withExtendedPrecisionSum(epSum);
            final String name = String.format("ZETA %2d %2d 2^%d %6s %6s %6s",
                n, m, b,
                powNp ? "powNp" : "",
                tail,
                epSum ? "EP sum" : "");
            final DoublePrecisionTestCase test = new DoublePrecisionTestCase() {
                @Override
                public double getTolerance() {
                    return 100;
                }

                @Override
                public double getRmsTolerance() {
                    return 10;
                }

                @Override
                public DoubleBinaryOperator getFunction() {
                    return (s, a) -> HurwitzZetaTest.zeta(s, a, c);
                }

                @Override
                public String[] getFilenames() {
                    // Use combined data from multiple resources in order to
                    // find N and M values.
                    return TEST_RESOURCES;
                }

                @Override
                public String toString() {
                    return name;
                }
            };
            assertFunction(test);

            // Check usage of M.
            // As N increases the convergence of the tail is faster (M decreases).
            // Johansson (2015) recommends N ~ M. However the test implementation does
            // not use Math.pow in the tail function and can evaluate higher M efficiently.
            int max = 0;
            for (int i = 0; i < M.length; i++) {
                if (M[i] != 0) {
                    max = i + 1;
                    // This is used for testing.
                    // CHECKSTYLE: stop regex
                    // System.out.printf("%s  M=%-2d  %d%n", name, max, M[i]);
                    // CHECKSTYLE: resume regex
                }
            }
            // n = 5 is too low and the tail does not converge
            Assertions.assertTrue(n < 6 || max <= m);
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
//        if (a < 0 && s < Integer.MAX_VALUE) {
//            //assertClose((x, y) -> HurwitzZetaTest.zetaNegative((int) x, y), s, a, z, 0);
//            assertClose((x, y) -> HurwitzZetaTest.zetaNegativeDD((int) x, y), s, a, z, 0);
//        }
    }

    // TODO - remove
    @Test
    void test() {
        assertClose(HurwitzZeta::value, 50, 2000, 3.6698119957034991027055454908981e-164, 0);
//        assertClose((x, y) -> HurwitzZetaTest.zetaNegativeDD((int) x, y),
//            5, -22.500000000921442, 0.00000148253693746789985363830295586415232, 0);
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
            Arguments.of(50, 2000000, 3.625262448698370064130045968154440065060042043605e-311, 1),
            Arguments.of(50, 3000000, 8.5283690506160333026898323949253691719144047850959e-320, 3),
            Arguments.of(50, 3500000, 4.4717013098147129077153665995697691837355569545233e-323, 0),
            Arguments.of(50, 4000000, 0, 0), // 6.43e-326
            Arguments.of(5, 10000000000.0, 2.5000000005000000000416666666666666666666666666529e-41, 0),
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
            Arguments.of(3, -7.499999999999999, 0.0077822558978654467309133842121074379966549149936996, 1),
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
            Arguments.of(3, -0.5 - 0x1p-53, 0.41439832211717073415946901920171, 3),
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
    @EnumSource(value = ZetaTestCase.class)
    @Order(1)
    void testZeta(ZetaTestCase tc) {
        assertFunction(tc);
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
    private static void assertFunction(DoublePrecisionTestCase tc) {
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
     * Find roots of the the zeta function.
     * This uses positive a parameters for convenience.
     */
    @ParameterizedTest
    @CsvSource({
        "3, 9, 0, 100",
        // Require evaluation in matlab
        // Note the root f(x) gets further from zero as |a| or s increases
        // and the root x -> half-integer
        "3, 9, 1000, 1100",
    })
    @Disabled("Used to generate test data")
    void testDataZetaRoots(int ls, int us, double la, double ua) throws IOException {
        // Validate arguments
        Assertions.assertTrue(ls > 2);
        Assertions.assertTrue(la >= 0);
        Assertions.assertNotEquals(la, la + 1, "a must be iterable with +1");
        Assertions.assertEquals(la, Math.floor(la), "lower a must be an integer");
        Assertions.assertEquals(ua, Math.floor(ua), "upper a must be an integer");
        Assertions.assertEquals(1, ls & 1, "s must be odd");
        // Lowest tolerance allowed
        final BrentSolver solver = new BrentSolver(0, 0, 0);
        try (PrintStream out = getPrintStream(
            String.format("hzeta_root_s%d_%d_na%s_%s.txt", ls, us, shortFormat(la), shortFormat(ua)))) {
            for (int s = ls; s <= us; s += 2) {
                final int ss = s;
                // Assume the function is optimised for accuracy
                final DoubleUnaryOperator f = x -> HurwitzZetaTest.zetaNegativeBD(ss, x);
                for (double a = la; a <= ua; a += 1) {
                    // a is integer: bracket -(a, a+1)
                    double min = -a - 1;
                    double max = -a;
                    double x = solver.findRoot(f, Math.nextUp(min), Math.nextDown(max));
                    // Check the solver found a bracket
                    double x0 = Math.nextDown(x);
                    double x1 = Math.nextUp(x);
                    double f0 = f.applyAsDouble(x0);
                    double f1 = f.applyAsDouble(x1);
                    Assertions.assertTrue(f0 * f1 <= 0);
                    out.printf("%s, %s%n", s, x0);
                    out.printf("%s, %s%n", s, x);
                    out.printf("%s, %s%n", s, x1);
                }
            }
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
        "2, 2, 0.9999999999988898, 1, false",
        // zeta called with a > 1.
        // Occurs when a > 2N with N the number of power terms in a zeta evaluation.
        "2, 2, 40, 41, true",
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
