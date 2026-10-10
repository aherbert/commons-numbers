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
 * Test the extended precision {@code (x+y)^z} power function for use in the Hurwitz zeta function.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class ExtPowTest {
    /** Threshold for z in the Taylor series approximation of expm1(z * t). */
    private static final double EXT_POW_LIMIT = 0x1p26;

    // Collect error statistics
    private static TestUtils.ErrorStatistics ERROR_POW = new TestUtils.ErrorStatistics();
    private static TestUtils.ErrorStatistics ERROR_POW1 = new TestUtils.ErrorStatistics();
    private static TestUtils.ErrorStatistics ERROR_EXT_POW = new TestUtils.ErrorStatistics();
    private static TestUtils.ErrorStatistics ERROR_EXT_POW1 = new TestUtils.ErrorStatistics();

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
                // 3-term version
                //final double actual = z * t * (1 + (z - 1) * t * 0.5 * (1 + (z - 2) * t * 0.333333333333333333333));
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
        // Hit use of large z above 2^26 when x+y < 1. Requires y to be zero.
        // for b in [24, 32]:
        //   for i in range(10):
        //     x, y, z = 1 - random.random()*2**-50, 0, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "0.9999999999999999, 0, -27904046.0762659, 1.00000000309797144820587975985",
        "0.9999999999999999, 0, -26189061.829705685, 1.00000000290756994789408988814",
        "0.9999999999999994, 0, -26816529.446500212, 1.00000001488616432682046736835",
        "0.9999999999999997, 0, -20183696.129322454, 1.00000000672252127203957913463",
        "0.9999999999999993, 0, -27542439.285582773, 1.00000001834695031782196243315",
        "0.9999999999999994, 0, -29270996.58120038, 1.0000000162486673110960489012",
        "0.9999999999999993, 0, -22306103.37358912, 1.00000001485884984340922629951",
        "0.9999999999999994, 0, -20465456.057365105, 1.00000001136061032670229494172",
        "0.9999999999999998, 0, -29863107.516385473, 1.00000000663094193229424143576",
        "0.9999999999999994, 0, -17674416.295603756, 1.00000000981127200722520995009",
        "0.9999999999999999, 0, -6267212798.305672, 1.00000069580063695959234472686",
        "0.9999999999999998, 0, -4347828409.366737, 1.00000096541230744982635703003",
        "0.9999999999999997, 0, -5351053069.490152, 1.00000178226028534570104611447",
        "0.9999999999999993, 0, -6538739205.55032, 1.00000435568477678009498410791",
        "0.9999999999999997, 0, -5128534294.039316, 1.00000170814651562724984966074",
        "0.9999999999999998, 0, -6410172652.885048, 1.0000014233452671660140657825",
        "0.9999999999999999, 0, -6387470093.090347, 1.0000007091518880934309139254",
        "1.0, 0, -6200000244.389243, 1.0",
        "0.9999999999999999, 0, -7137073172.4818535, 1.00000079237461038098229296614",
        "0.9999999999999998, 0, -6888409257.1306505, 1.00000152953528179940023108887",
        // for b in [24, 27]:
        //   for i in range(10):
        //     x, y, z = 1 - random.random()*2**-20, 0, -random.uniform(1, 2)*2**b
        //     print(f'"{x}, {y}, {z}, {mp.power(mpf(x) + mpf(y), z)}",')
        "0.9999998176777922, 0, -24503492.36240287, 87.1413256102570175334100378276",
        "0.9999997192695803, 0, -26213578.825065367, 1570.18702905628443030634027023",
        "0.9999995391611006, 0, -17185809.10998093, 2751.47161429117464715445932581",
        "0.9999992920294289, 0, -27875117.813870788, 372136046.871827762829410649539",
        "0.9999993894259576, 0, -22812664.00076673, 1119983.74524076262533638078138",
        "0.9999991672085344, 0, -17281094.72699289, 1778986.17759816267713600130245",
        "0.999999225950329, 0, -20195840.71029392, 6153858.91032657693878404725672",
        "0.9999995126732188, 0, -26794198.447530862, 468613.523779248629574003823484",
        "0.9999996806954553, 0, -30502460.003091864, 16976.3360842443684112091728549",
        "0.9999997476939285, 0, -29780964.18247216, 1833.38456308991800358436389874",
        "0.9999995476091361, 0, -147942488.9394067, 1.16518260583501182403014659931e+29",
        "0.9999998972731136, 0, -228186919.62698698, 15144950060.8963227988865318094",
        "0.9999991125878169, 0, -176589043.80780017, 1.14059730872805806276858437621e+68",
        "0.9999996007572814, 0, -264741821.61460328, 8.00396056866403549352225545496e+45",
        "0.9999997466684181, 0, -167936136.3068788, 2995169395663905872.15728469694",
        "0.9999997869531986, 0, -201266536.4068937, 4189849056864775697.99106627307",
        "0.9999993890593745, 0, -197228338.05651915, 2.13916635023629946194570550086e+52",
        "0.9999991197263971, 0, -208585172.23997986, 5.51724999544252598823594276745e+79",
        "0.9999990604556022, 0, -153412742.56403875, 3.96646305491769763354735089792e+62",
        "0.9999999172736319, 0, -146464465.42655826, 182859.57362861956675558067228",
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
        // Handle underflow and overflow. No support for -infinity. Use MAX_VALUE
        "1.0, 0, -1.7976931348623157E308, 1.0",
        "0.999999999999, 0, -1.7976931348623157E308, Infinity",
        "1.000000000001, 0, -1.7976931348623157E308, 0",
    })
    @Order(1)
    void testExtPow(double x, double y, double z, double expected) {
        // Check data is suitable for the zeta use case
        Assertions.assertTrue(z < -1, "Expected domain of zeta(s, a) is s > 1");

        // Math.pow is not expected to be ulp exact for all cases.
        // When it is we collect the error separately to obtain the RMS for this subset.
        final long stdUlp = TestUtils.assertEquals(expected, Math.pow(x + y, z), 10000000000L,
            ERROR_POW::add, () -> String.format("%s, %s, %s", x, y, z));
        final long extUlp = TestUtils.assertEquals(expected, ExtPowTest.extPow(x, y, z), 1,
            ERROR_EXT_POW::add, () -> String.format("%s, %s, %s", x, y, z));
        if (Math.abs(stdUlp) <= 1) {
            ERROR_POW1.add(stdUlp);
            ERROR_EXT_POW1.add(extUlp);
        }
    }

    @Test
    @Order(2)
    void testExtPowRMS() {
        // Some cases are very wrong using Math.pow(x+y, z)
        Assertions.assertTrue(ERROR_POW.getRMS() > 1e6);
        // When x+y is exact then Math.pow is ulp accurate
        Assertions.assertTrue(ERROR_POW1.getRMS() < 0.5);
        // Extended precision is ulp accurate for all x+y cases
        Assertions.assertTrue(ERROR_EXT_POW.getRMS() < 0.5);
        // Where Math.pow was ulp accurate, the extended precision lowered the RMS.
        // Note there are some cases where Math.pow is closer to the correct result
        // but on average the extended precision method is closer.
        Assertions.assertTrue(ERROR_EXT_POW1.getRMS() < ERROR_POW1.getRMS() * 0.75);
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
    static double extPow(double x, double y, double z) {
        // (s+ss)^z = s^z * (1+ss/s)^z
        //          = s^z * s^z * [ (1+ss/s)^z - 1 ]
        //          = s^z * exp(z*log1p(ss/s))
        //          = s^z + s^z * expm1(z*log1p(ss/s))
        //
        // ss/s < 2 * machine epsilon : log1p(ss/s) ~ ss/s
        //
        // Taylor series: (1 + x)^n = 1 + nx + n(n-1)x^2/2! + n(n-1)(n-2)x^3/3! + ...
        // for |x| < 1 and real n

        final DD s = DD.ofSum(x, y);
        if (s.lo() == 0) {
            // Skip attempted rounding
            return Math.pow(s.hi(), z);
        }

        // Remove the lowest set bit from (x+y) if present.
        // This removes some magnitude from the pow result and allows rounding to be
        // performed by this method with the round-off from (x+y).
        // This makes the most difference when s.lo() is much smaller than ulp(s.hi()) / 2.
        final double hx = highPart(s.hi());
        final double lx = s.hi() - hx;

        final double r = Math.pow(hx, z);
        // Do not adjust infinity (avoids 0 * infinity == NaN)
        if (!Double.isFinite(r)) {
            return r;
        }

        // t < 2^-52
        double t = (lx + s.lo()) / hx;

        // Limit of Taylor series where it does not equal expm1(z * t)
        // is at approximately 2^-26.
        // z * t < 2^-26  =>  z < 2^-26 / 2^-52 = 2^26
        if (z < -EXT_POW_LIMIT) {
            // z is expected to be small (practical range for a zeta result is -[1, 1075))
            // and this path is unlikely
            t = Math.expm1(z * t);
        } else {
            t = z * t * (1 + (z - 1) * t * 0.5);
        }

        // Note: This could return a double-double result for extended summation
        return r + r * t;
    }

    /**
     * Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) creates
     * a big value from which to derive the two split parts.
     * <pre>
     * c = (2^s + 1) * a
     * a_big = c - a
     * a_hi = c - a_big
     * a_lo = a - a_hi
     * a = a_hi + a_lo
     * </pre>
     *
     * <p>The multiplicand allows a p-bit value to be split into
     * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
     * Combined they have (p-1) bits of significand but the sign bit of {@code a_lo}
     * contains a bit of information. This uses s = 1 to create a 52-bit value and
     * the least significant bit.
     *
     * <p>This conversion does not use scaling and the result of overflow is NaN. Overflow
     * may occur when the exponent of the input value is above 1021.
     *
     * <p>Splitting a NaN or infinite value will return NaN.
     *
     * @param value Value.
     * @return the high part of the value.
     * @see Math#getExponent(double)
     */
    static double highPart(double value) {
        final double c = 3 * value;
        return c - (c - value);
    }
}
