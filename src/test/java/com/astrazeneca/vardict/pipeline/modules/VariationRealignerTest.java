package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.GlobalReadOnlyScope;
import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.pipeline.modules.VariationRealigner;
import com.astrazeneca.vardict.variations.Variation;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static com.astrazeneca.vardict.pipeline.modules.VariationRealigner.*;

public class VariationRealignerTest {

    @Test
    public void testIsLowComplexSeq() throws Exception {
        Assert.assertTrue(VariationRealigner.isLowComplexSeq("AAAAAAAAA"));
        Assert.assertFalse(VariationRealigner.isLowComplexSeq("ATATATATATAT"));
        Assert.assertFalse(VariationRealigner.isLowComplexSeq("ACGTACGTACGT"));
        Assert.assertTrue(VariationRealigner.isLowComplexSeq("CCCCCCCCGA"));
    }

    @Test
    public void testIsMatch() throws Exception {
        Assert.assertTrue(VariationRealigner.isMatch("AAAAAAAAA","AAAAAAAAA", 1, false));
        Assert.assertTrue(VariationRealigner.isMatch("AAAAAAAAA","AAAAAAAAA", -1, false));
        Assert.assertTrue(VariationRealigner.isMatch("ACGTACGTACGTACGT","AAGTACTTACGTACGT", 1, false));
        Assert.assertTrue(VariationRealigner.isMatch("ACGTACGTACGTACGT", "AAGTACTTACGTACGT", 1, false));
        Assert.assertTrue(VariationRealigner.isMatch("ACGTACGTACGTACGT","AAGTACTTACGT", 1, false));

        Assert.assertFalse(VariationRealigner.isMatch("ACGTACGTACGT","AAGTACTTACGT", 1, false));
        Assert.assertFalse(VariationRealigner.isMatch("ACGTCAGCAT","ACGACTGACT", 1, false));
    }

    @DataProvider(name = "isRelevantVariantTestDataProvider")
    public Object[][] isRelevantVariantTestDataProvider() {
        Configuration conf = new Configuration();
        conf.goodq = 9;
        GlobalReadOnlyScope.init(conf, null, null, null, "");
        return new Object[][] {
                // all condition satisfied
                {constructVariation(1, 1, 10), true},
                // variation is null
                {null, false},
                // varsCount in variation is 0, so return false
                {constructVariation(0, 1, 10), false},
                // tv.meanQuality / tv.varsCount < conf.goodq, so return false
                {constructVariation(10, 1, 10), false}
        };

    }

    @Test(dataProvider = "isRelevantVariantTestDataProvider")
    public void isRelevantVariantTest(Variation variation, boolean expected) {
        Assert.assertTrue(
                expected == isRelevantVariant(
                        10,
                        8,
                        3,
                        3,
                        3,
                        variation
                )
        );
    }

    @DataProvider(name = "find35MatchTestDataProvider")
    public Object[][] find35MatchTestDataProvoder() {
        return new Object[][] {
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCATGCA", new ExpectedMatchInfo(1, 0, 24)},
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCCTGCATGCATGCA", new ExpectedMatchInfo(1, 0, 23)},
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCATGCA", new ExpectedMatchInfo(1, 0, 22)},
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCAAAAA", new ExpectedMatchInfo(13, 0, 12)},
                {"ACGTACGTACGTACGTACGTACGT", "AAAATGCATGCATGGGTGCAAAAA", new ExpectedMatchInfo(11, 14, 10)},
        };
    }

    @Test(dataProvider = "find35MatchTestDataProvider")
    public void find35MatchTest(String seq5, String seq3, ExpectedMatchInfo expected) {
        Tuple.Tuple3<Integer, Integer, Integer> match = find35match(seq5, seq3);
        Assert.assertEquals(match._1, expected.p3);
        Assert.assertEquals(match._2, expected.p5);
        Assert.assertEquals(match._3, expected.match);

    }

    private Variation constructVariation(int count, int pmean, int qmean) {
        Variation variation = new Variation();
        variation.varsCount = count;
        variation.meanPosition = pmean;
        variation.meanQuality = qmean;
        return variation;
    }

    private class ExpectedMatchInfo {
        private final Integer p5;
        private final Integer p3;
        private final Integer match;

        public ExpectedMatchInfo(int p5, int p3, int match) {
            this.p5 = p5;
            this.p3 = p3;
            this.match = match;
        }
    }
}
