package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.Match35;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.HashMap;

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
        GlobalReadOnlyScope.init(new Configuration(), null, null, null, "", new HashMap<>(), new HashMap<>());

        Assert.assertTrue(new VariationRealigner().ismatch("AAAAAAAAA","AAAAAAAAA", 1));
        Assert.assertTrue(new VariationRealigner().ismatch("AAAAAAAAA","AAAAAAAAA", -1));
        Assert.assertTrue(new VariationRealigner().ismatch("ACGTACGTACGTACGT","AAGTACTTACGTACGT", 1));
        Assert.assertTrue(new VariationRealigner().ismatch("ACGTACGTACGTACGT", "AAGTACTTACGTACGT", 1));
        Assert.assertTrue(new VariationRealigner().ismatch("ACGTACGTACGTACGT","AAGTACTTACGT", 1));

        Assert.assertFalse(new VariationRealigner().ismatch("ACGTACGTACGT","AAGTACTTACGT", 1));
        Assert.assertFalse(new VariationRealigner().ismatch("ACGTCAGCAT","ACGACTGACT", 1));
    }

        @DataProvider(name = "find35MatchTestDataProvider")
    public Object[][] find35MatchTestDataProvoder() {
        return new Object[][] {
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCATGCA", new Match35(0, 1, 24)},
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCCTGCATGCATGCA", new Match35(0, 1, 23)},
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCATGCA", new Match35(0, 1, 22)},
                {"ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCAAAAA", new Match35(0, 13, 12)},
                {"ACGTACGTACGTACGTACGTACGT", "AAAATGCATGCATGGGTGCAAAAA", new Match35(14, 11, 10)},
        };
    }

    @Test(dataProvider = "find35MatchTestDataProvider")
    public void find35MatchTest(String seq5, String seq3, Match35 expected) {
        GlobalReadOnlyScope.init(new Configuration(), null, null, null, "", new HashMap<>(), new HashMap<>());
        Match35 match35 = new VariationRealigner().find35match(seq5, seq3);
        Assert.assertEquals(match35.matched5end, expected.matched5end);
        Assert.assertEquals(match35.matched3End, expected.matched3End);
        Assert.assertEquals(match35.maxMatchedLength, expected.maxMatchedLength);

    }
}
