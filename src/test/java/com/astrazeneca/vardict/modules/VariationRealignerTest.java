package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class VariationRealignerTest {

    @AfterMethod
    public void cleanUp() {
        GlobalReadOnlyScope.clear();
    }

    @Test
    public void testIsLowComplexSeq() {
        Assert.assertTrue(VariationRealigner.islowcomplexseq("AAAAAAAAA"));
        Assert.assertTrue(VariationRealigner.islowcomplexseq("ATATATATATAT"));
        Assert.assertTrue(VariationRealigner.islowcomplexseq("CCCCCCCCGA"));
        Assert.assertFalse(VariationRealigner.islowcomplexseq("ACGTACGTACGT"));
        Assert.assertFalse(VariationRealigner.islowcomplexseq("CCGTAACGGGGT"));
    }

    @Test
    public void testIsMatch() {
        GlobalReadOnlyScope.init(new Configuration(), null, null, null, "");

        Assert.assertTrue(VariationRealigner.ismatch("AAAAAAAAA","AAAAAAAAA", 1));
        Assert.assertTrue(VariationRealigner.ismatch("AAAAAAAAA","AAAAAAAAA", -1));
        Assert.assertTrue(VariationRealigner.ismatch("ACGTACGTACGTACGT","AAGTACTTACGTACGT", 1));
        Assert.assertTrue(VariationRealigner.ismatch("ACGTACGTACGTACGT", "AAGTACTTACGTACGT", 1));
        Assert.assertTrue(VariationRealigner.ismatch("ACGTACGTACGTACGT","AAGTACTTACGT", 1));

        Assert.assertFalse(VariationRealigner.ismatch("ACGTACGTACGT","AAGTACTTACGT", 1));
        Assert.assertFalse(VariationRealigner.ismatch("ACGTCAGCAT","ACGACTGACT", 1));
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
        GlobalReadOnlyScope.init(new Configuration(), null, null, null, "");
        Tuple.Tuple3<Integer, Integer, Integer> match = new VariationRealigner().find35match(seq5, seq3);
        Assert.assertEquals(match._1, expected.p3);
        Assert.assertEquals(match._2, expected.p5);
        Assert.assertEquals(match._3, expected.match);

    }

    private class ExpectedMatchInfo {
        private final Integer p5;
        private final Integer p3;
        private final Integer match;

        ExpectedMatchInfo(int p5, int p3, int match) {
            this.p5 = p5;
            this.p3 = p3;
            this.match = match;
        }
    }
}
