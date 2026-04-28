package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static com.astrazeneca.vardict.data.Patterns.threeIndelsPattern;
import static com.astrazeneca.vardict.modules.CigarParser.getInsertionDeletionLength;

public class CigarModifierTest {
    @AfterMethod
    public void cleanUp () {
        GlobalReadOnlyScope.clear();
    }

    @Test
    public void findOffsetTest() {
        Configuration config = new Configuration();
        config.goodq = 23;
        config.vext = 3;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());

        int referencePosition = 1;
        int readPosition = 2;
        int cigarLength = 3;
        String querySequence = "ACGTACGT";
        String queryQuality = "<<<<<<<<";
        Map<Integer, Integer> refCoverage = new HashMap<>();
        Map<Integer, Character> reference = new HashMap<Integer, Character>() {{
            put(1, 'A');
            put(2, 'A');
        }};
        CigarParser cigarParser = new CigarParser(false);
        CigarParser.Offset result = cigarParser.findOffset(
                referencePosition,
                readPosition,
                cigarLength,
                querySequence,
                queryQuality,
                reference,
                refCoverage);
        CigarParser.Offset expectedResult = cigarParser.new Offset(2, "GT", "<<", 2);
        Assert.assertEquals(result.offset, expectedResult.offset);
        Assert.assertEquals(result.sequence, expectedResult.sequence);
        Assert.assertEquals(result.qualitySequence, expectedResult.qualitySequence);
        Assert.assertEquals(result.offsetNumberOfMismatches, expectedResult.offsetNumberOfMismatches);
    }

    @DataProvider(name = "variationsForSub")
    public Object[][] variationsForSub() {
        return new Object[][] {
                {
                        new Variation() {{
                            varsCount = 4;
                            varsCountOnReverse = 4;
                            varsCountOnForward = 4;
                            meanPosition = 9;
                            meanQuality = 10.5;
                            meanMappingQuality = 31;
                            numberOfMismatches = 8;
                            highQualityReadsCount = 45;
                            lowQualityReadsCount = 34;
                        }},
                        true,
                        22.5
                },
                {
                        new Variation() {{
                            varsCount = 4;
                            varsCountOnReverse = 5;
                            varsCountOnForward = 3;
                            meanPosition = 9;
                            meanQuality = 10.5;
                            meanMappingQuality = 31;
                            numberOfMismatches = 8;
                            highQualityReadsCount = 45;
                            lowQualityReadsCount = 34;
                        }},
                        false,
                        22.5
                },
                {
                        new Variation() {{
                            varsCount = 4;
                            varsCountOnReverse = 5;
                            varsCountOnForward = 3;
                            meanPosition = 9;
                            meanQuality = 10.5;
                            meanMappingQuality = 31;
                            numberOfMismatches = 8;
                            highQualityReadsCount = 44;
                            lowQualityReadsCount = 35;
                        }},
                        false,
                        0.0
                }
        };
    }

    @Test(dataProvider = "variationsForSub")
    public void subCntTest(Object expectedVariationObject, boolean direction, double goodQuality) {
        Configuration config = new Configuration();
        config.goodq = goodQuality;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());

        Variation variation = new Variation() {{
            varsCount = 5;
            varsCountOnReverse = 5;
            varsCountOnForward = 4;
            meanPosition = 12;
            meanQuality = 12.0;
            meanMappingQuality = 36;
            numberOfMismatches = 10;
            highQualityReadsCount = 45;
            lowQualityReadsCount = 35;
        }};
        new CigarParser(false).subCnt(variation, direction, 3, 1.5, 5, 2);
        Variation expectedVariation = (Variation) expectedVariationObject;
        Assert.assertEquals(variation.varsCount, expectedVariation.varsCount);
        Assert.assertEquals(variation.varsCountOnReverse, expectedVariation.varsCountOnReverse);
        Assert.assertEquals(variation.varsCountOnForward, expectedVariation.varsCountOnForward);
        Assert.assertEquals(variation.meanPosition, expectedVariation.meanPosition);
        Assert.assertEquals(variation.meanQuality, expectedVariation.meanQuality);
        Assert.assertEquals(variation.meanMappingQuality, expectedVariation.meanMappingQuality);
        Assert.assertEquals(variation.numberOfMismatches, expectedVariation.numberOfMismatches);
        Assert.assertEquals(variation.highQualityReadsCount, expectedVariation.highQualityReadsCount);
        Assert.assertEquals(variation.lowQualityReadsCount, expectedVariation.lowQualityReadsCount);
    }

    @DataProvider(name = "variationsForAdd")
    public Object[][] variationsForAdd() {
        return new Object[][] {
                {
                        new Variation() {{
                            varsCount = 6;
                            varsCountOnReverse = 6;
                            varsCountOnForward = 4;
                            meanPosition = 15;
                            meanQuality = 13.5;
                            meanMappingQuality = 41;
                            numberOfMismatches = 12;
                            highQualityReadsCount = 45;
                            lowQualityReadsCount = 36;
                        }},
                        true,
                        22.5
                },
                {
                        new Variation() {{
                            varsCount = 6;
                            varsCountOnReverse = 5;
                            varsCountOnForward = 5;
                            meanPosition = 15;
                            meanQuality = 13.5;
                            meanMappingQuality = 41;
                            numberOfMismatches = 12;
                            highQualityReadsCount = 45;
                            lowQualityReadsCount = 36;
                        }},
                        false,
                        22.5
                },
                {
                        new Variation() {{
                            varsCount = 6;
                            varsCountOnReverse = 5;
                            varsCountOnForward = 5;
                            meanPosition = 15;
                            meanQuality = 13.5;
                            meanMappingQuality = 41;
                            numberOfMismatches = 12;
                            highQualityReadsCount = 46;
                            lowQualityReadsCount = 35;
                        }},
                        false,
                        0.0
                }
        };
    }

    @Test(dataProvider = "variationsForAdd")
    public void addCntTest(Object expectedVariationObject, boolean direction, double goodQuality) {
        Configuration config = new Configuration();
        config.goodq = goodQuality;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());

        Variation variation = new Variation() {{
            varsCount = 5;
            varsCountOnReverse = 5;
            varsCountOnForward = 4;
            meanPosition = 12;
            meanQuality = 12.0;
            meanMappingQuality = 36;
            numberOfMismatches = 10;
            highQualityReadsCount = 45;
            lowQualityReadsCount = 35;
        }};
        new CigarParser(false).addCnt(variation, direction, 3, 1.5, 5, 2);
        Variation expectedVariation = (Variation) expectedVariationObject;
        Assert.assertEquals(variation.varsCount, expectedVariation.varsCount);
        Assert.assertEquals(variation.varsCountOnReverse, expectedVariation.varsCountOnReverse);
        Assert.assertEquals(variation.varsCountOnForward, expectedVariation.varsCountOnForward);
        Assert.assertEquals(variation.meanPosition, expectedVariation.meanPosition);
        Assert.assertEquals(variation.meanQuality, expectedVariation.meanQuality);
        Assert.assertEquals(variation.meanMappingQuality, expectedVariation.meanMappingQuality);
        Assert.assertEquals(variation.numberOfMismatches, expectedVariation.numberOfMismatches);
        Assert.assertEquals(variation.highQualityReadsCount, expectedVariation.highQualityReadsCount);
        Assert.assertEquals(variation.lowQualityReadsCount, expectedVariation.lowQualityReadsCount);
    }

    @DataProvider(name = "cigar")
    public Object[][] cigar() {
        return new Object[][] {
                {
                        new Cigar(new ArrayList<CigarElement>() {{
                            add(new CigarElement(1, CigarOperator.M));
                            add(new CigarElement(2, CigarOperator.S));
                            add(new CigarElement(4, CigarOperator.I));
                            add(new CigarElement(8, CigarOperator.D));
                            add(new CigarElement(16, CigarOperator.N));
                            add(new CigarElement(32, CigarOperator.H));
                        }})
                }
        };
    }

    @Test(dataProvider = "cigar")
    public void getInsertionDeletionLengthTest(Object cigar) {
        Assert.assertEquals(getInsertionDeletionLength((Cigar) cigar), 12);
    }

    @Test(dataProvider = "cigar")
    public void getMatchInsertionLengthTest(Object cigar) {
        Assert.assertEquals(CigarParser.getMatchInsertionLength((Cigar) cigar), 5);
    }

    @Test(dataProvider = "cigar")
    public void getSoftClippedLengthTest(Object cigar) {
        Assert.assertEquals(CigarParser.getSoftClippedLength((Cigar) cigar), 7);
    }

    @Test(dataProvider = "cigar")
    public void getCigarOperatorTest(Object cigarObject) {
        Cigar cigar = (Cigar) cigarObject;
        CigarOperator[] operators = new CigarOperator[] {
                CigarOperator.M,
                CigarOperator.S,
                CigarOperator.I,
                CigarOperator.D,
                CigarOperator.N,
                CigarOperator.H,
        };
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            Assert.assertEquals(CigarParser.getCigarOperator(cigar, i), operators[i]);
        }
    }

    @DataProvider(name = "dataForIsBEGIN_ATGC_AMP_ATGCs_END_Test")
    public Object[][] dataForIsBEGIN_ATGC_AMP_ATGCs_END_Test() {
        return new Object[][] {
                {"A&ACGT", true},
                {"A&", false},
                {"A&ASGT", false}
        };
    }

    @Test(dataProvider = "dataForIsBEGIN_ATGC_AMP_ATGCs_END_Test")
    public void isBEGIN_ATGC_AMP_ATGCs_END(final String s, final boolean result) {
        Assert.assertEquals(new CigarParser(false).isBEGIN_ATGC_AMP_ATGCs_END(s), result);
    }

    @DataProvider(name = "threeIndelsCigar")
    public Object[][] threeIndelsCigar() {
        return new Object[][] {
                {"15M1I1M1I2M1I27M", "27M3I18M", 3, 1,
                        "GTGAAAACTGGAAACAAAAAAAAAAAAAAACAGCCTCTTCTTGGAAAG",
                        "C@@FFFFFHHGDHIJIIJGJJJIIGJJIHD6=3@A;;>;(;ACCCDCB",
                        "GTGAAAACTGGAAACAAAAAAAAAAAACAGCCTCTTCTTGG" },
                {"14M3I2M2I2M2D26M1834N23M", "18M2D5I26M1834N23M", 7, 1,
                        "GACTCTTTTGTAGGTGCTGATCTTCTTAAGAAATTTAAATTTCTAAAAGGTGCTACACTGTGTGTCATCTGC",
                        "AAAAAEEEEEEEE6EEEEAEEEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEAE",
                        "GACTCTTTTGTAGGTGCTNNTCTTAAGAAATTTAAATTTCT" },
                {"2S2M2I2M2D26M1834N10M", "6S28M1834N10M", 4, 1,
                        "GCTGATCTTCTTAAGAAATTTAAATTTCTAAAAGGTGCTACACT",
                        "AAAA6/EAEEEEEEEEEAEEAEE/6AEEE6EEA/EE/EEAA/EE",
                        "GGTGCTGATCTTAAGAAATTTAAATTTCTAAAAGNNNNNNN" },
                {"34M3I2M2I2M2D26M1834N4M", "44M3I22M1834N4M", 7,  1,
                        "GACTCCTTTTGTACCTCTGTGACTCTTTTGTAGGTGCTGATCTTCTTAAGAAATTTAAATTTCTAAAAGGTGC",
                        "EEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEAAAAA",
                        "GACTCCTTTTGTACCTCTGTGACTCTTTTGTAGGTGCTGATCTTAAGAAATTTAAATTTCTAAAAGNNNNNNNNNNNNNNNN" },
                {"11M1I2M1I8M1I51M", "21M3I51M", 3, 1,
                        "GTTTTTTGTTTCCTTTTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTT",
                        "AAAAAEEEEEEEEEEEEEEEEEEEEEEEEA<AAEEEAEEEEEE6EE/EEEEE6AEEEAEEEEEA<EEEEEEEEEA",
                        "GTTTTTTGTTTCCTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTTNNNNNNNNNN" },
                {"10M2I2M1I8M1I51M", "21M3I50M", 4, 1,
                        "GTTTTTTGTTTCCTTTTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTT",
                        "AAAAAEEEEEEEEEEEEEEEEEEEEEEEEA<AAEEEAEEEEEE6EE/EEEEE6AEEEAEEEEEA<EEEEEEEEEA",
                        "GTTTTTTGTTTCCTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTTNNNNNNNNNN" },
                {"16M1I1M1I2M1D3M1I50M", "23M2I50M", 4, 1,
                        "CCTGTTCAAAGAAACCAAGAGGGAAAAAAAAAAGTGCAGAAGAAAAGTAATTTGCACTTAAAAATGAAAACAAAC",
                        "EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/AEEEEEEEEEEAAAAA",
                        "CCTGTTCAAAGAAACCAAGAGGGAAAAAAAAGTGCAGAAGAAAAGTAATTTGCACTTAAAAATGAAAACAAACCACTTTTCTC" },

        };
    }

    @Test(dataProvider = "threeIndelsCigar")
    public void threeIndels(final String initialCigar, final String expected, final int indel, final int position,
                            String querySequence, final String queryQuality, String refSequence) {
        Map<Integer, Character> reference = new HashMap<>();
        int i = 1;
        for (String base : refSequence.split("")) {
            reference.put(i, base.charAt(0));
            i++;
        }
        Reference ref = new Reference();
        ref.referenceSequences = reference;
        Region region = new Region("chr1", 1, 100, "");
        CigarModifier cigarModifier = new CigarModifier(position, initialCigar, querySequence,
                queryQuality, ref, indel, region, 150);
        String actual = cigarModifier.modifyCigar().cigar;
        Assert.assertEquals(actual, expected);
    }

    /**
     * Regression test for a malformed CIGAR produced by the merge logic when a read aligns
     * inside a long homopolymer with three closely-spaced deletions.
     *
     * Reproducer is taken verbatim from a duplex consensus read observed in production:
     *   read name : duplex_GTGTC-GTCAC_140434267-140434623_817012 (mate2)
     *   chr/pos   : 7:140434471 (hg19)
     *   cigar     : 104M1I2M2D1M1D5M1D2M9D27M
     *
     * Without a fix, threeDeletions()/captureMisSoftly3Mismatches() let the matching-base
     * scan run past the read end, producing the malformed intermediate CIGAR
     *   104M15D12I-1M9D27M
     * (note the negative -1M operator) and then throwing
     *   StringIndexOutOfBoundsException: Index 143 out of bounds for length 142
     * inside VariationUtils.isHasAndNotEquals.
     *
     * The expected behaviour is: modifyCigar() must NOT throw, and must return a CIGAR
     * that htsjdk's TextCigarCodec can decode (i.e. all operator lengths positive).
     */
    @Test
    public void modifyCigarOnHomopolymerThreeDeletionsDoesNotProduceMalformedCigar() {
        // Read sequence and qualities are 142 nt (the bug reports "Index 143 out of bounds for length 142").
        final String querySequence =
                "AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTAC"
              + "AAAAAAAAAAAAAAAAAAAAGAAAAAAGAAAAAAAAAG";
        final String queryQuality =
                "BBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
              + "=74??44444444FFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        Assert.assertEquals(querySequence.length(), 142);
        Assert.assertEquals(queryQuality.length(), 142);

        // Reference window covering position .. position + cigarRefSpan(154) with a small margin,
        // taken from human_g1k_v37.fa at chr7:140434471-140434670. Reference is 1-indexed in the test,
        // so reference position 1 corresponds to the read alignment start.
        final int position = 1;
        final String refSequence =
                "AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTAC"
              + "AAAAAAAAAAAGAAAAAAAAAAGAAAAAAAAAGAAAAAAGAAAAAAAAAGAAAGAAAGAAAAAGAAAAAACAGAAAGAAGAATGAAAATCTGGGTG";

        Map<Integer, Character> reference = new HashMap<>();
        for (int i = 0; i < refSequence.length(); i++) {
            reference.put(i + 1, refSequence.charAt(i));
        }
        Reference ref = new Reference();
        ref.referenceSequences = reference;

        // GlobalReadOnlyScope is required by some merge-helper paths.
        Configuration config = new Configuration();
        config.goodq = 23;
        config.vext = 3;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());

        Region region = new Region("7", 1, 700, "");

        final String initialCigar = "104M1I2M2D1M1D5M1D2M9D27M";

        CigarModifier cigarModifier = new CigarModifier(
                position, initialCigar, querySequence, queryQuality,
                ref, /*indel=*/3, region, /*maxReadLength=*/150);

        String modifiedCigar = cigarModifier.modifyCigar().cigar;

        // Hard requirement #1: htsjdk must be able to decode the result. This currently fails
        // with "Malformed CIGAR string: 104M15D12I-1M9D27M".
        try {
            TextCigarCodec.decode(modifiedCigar);
        } catch (IllegalArgumentException e) {
            Assert.fail("modifyCigar() returned a CIGAR htsjdk cannot decode: '"
                    + modifiedCigar + "' (" + e.getMessage() + ")");
        }

        // Hard requirement #2: every operator length must be positive.
        for (CigarElement el : TextCigarCodec.decode(modifiedCigar).getCigarElements()) {
            Assert.assertTrue(el.getLength() > 0,
                    "modifyCigar() produced non-positive length for operator "
                            + el.getOperator() + " in CIGAR '" + modifiedCigar + "'");
        }
    }
}