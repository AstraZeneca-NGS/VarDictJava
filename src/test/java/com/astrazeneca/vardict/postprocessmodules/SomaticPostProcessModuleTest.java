package com.astrazeneca.vardict.postprocessmodules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.printers.SystemOutVariantPrinter;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.Mockito;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.AfterTest;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import static org.testng.Assert.assertEquals;

public class SomaticPostProcessModuleTest {

    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private Region region = new Region("1", 1, 10, "gene_name");
    private String sample = "test_bam";

    @AfterMethod
    public void clearOutContent() {
        outContent.reset();
        GlobalReadOnlyScope.clear();
    }

    @AfterTest
    public void closeOutContent() throws IOException {
        outContent.close();
    }

    Variant badVariant = new Variant() {{
        descriptionString = "T";
        totalPosCoverage = 0;
        positionCoverage = 4;
        varsCountOnReverse = 5;
        varsCountOnForward = 3;
        strandBiasFlag = "2";
        frequency = 0.4;
        meanPosition = 2.2;
        meanQuality = 2.6;
        meanMappingQuality = 7.8;
        numberOfMismatches = 2.0;
        hicnt = 44;
        highQualityToLowQualityRatio = 1.2571428571428571;
    }};

    Variant goodVariant = new Variant() {{
        refallele = "T";
        varallele = "A";
        descriptionString = "T";
        genotype = "T/A";
        totalPosCoverage = 0;
        positionCoverage = 4;
        varsCountOnReverse = 5;
        varsCountOnForward = 3;
        strandBiasFlag = "2";
        frequency = 0.4;
        meanPosition = 9;
        meanQuality = 26;
        meanMappingQuality = 7.8;
        numberOfMismatches = 2.0;
        hicnt = 44;
        highQualityToLowQualityRatio = 2.5;
    }};

    Variant refVariant = new Variant() {{
        refallele = "T";
        varallele = "T";
        descriptionString = "T";
        genotype = "T/T";
        totalPosCoverage = 0;
        positionCoverage = 4;
        varsCountOnReverse = 5;
        varsCountOnForward = 3;
        strandBiasFlag = "2";
        frequency = 0.4;
        meanPosition = 9;
        meanQuality = 26;
        meanMappingQuality = 7.8;
        numberOfMismatches = 2.0;
        hicnt = 44;
        highQualityToLowQualityRatio = 2.5;
    }};

    @Test
    public void testNoVariantsForBothSamples() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{ }},
                new Vars(){{ }},
                region, new HashSet<>());
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testBadVariantsInBothSamples() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{variants.add(badVariant); referenceVariant = refVariant; }},
                new Vars(){{variants.add(badVariant); referenceVariant = refVariant; }},
                region, new HashSet<>());
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testBadVariantInFirstNoInSecondSample() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{ variants.add(badVariant); referenceVariant = refVariant; }},
                new Vars(){{ }},
                region, new HashSet<>());
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testBadVariantInSecondNoInFirstSample() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{  }},
                new Vars(){{ variants.add(badVariant); referenceVariant = refVariant; }},
                region, new HashSet<>());
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testGoodVariantInFirstNoInSecondSample() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{ variants.add(goodVariant); referenceVariant = refVariant; }},
                new Vars(){{ }},
                region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t" +
                "7.8\t2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "StrongSomatic\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    @Test
    public void testGoodVariantInSecondNoInFirstSample() {
        SomaticPostProcessModule module = runMethod();
        Map<String, Variant> varDesc = new HashMap<>();
        varDesc.put("T", goodVariant);

        module.callingForBothSamples(1,
                new Vars(){{ }},
                new Vars(){{ variants.add(goodVariant); referenceVariant = refVariant; varDescriptionStringToVariants = varDesc; }},
                region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" +
                "0\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "StrongLOH\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    @Test
    public void testGoodVariantCallingForOneSample() {
        SomaticPostProcessModule module = runMethod();
        Map<String, Variant> varDesc = new HashMap<>();
        varDesc.put("T", goodVariant);

        module.callingForOneSample(new Vars(){{ variants.add(goodVariant); referenceVariant = refVariant; varDescriptionStringToVariants = varDesc; }},
                true, "Deletion", region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" +
                "0\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "Deletion\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    @Test
    public void testGoodVariantsInBothSamples() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{variants.add(goodVariant); referenceVariant = refVariant; }},
                new Vars(){{variants.add(goodVariant); referenceVariant = refVariant; }},
                region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8\t" +
                "2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "StrongSomatic\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    @Test
    public void testGoodVariantsInBothSamplesWithDesc() {
        SomaticPostProcessModule module = runMethod();
        Map<String, Variant> varDesc = new HashMap<>();
        varDesc.put("T", goodVariant);

        module.callingForBothSamples(1,
                new Vars(){{ variants.add(goodVariant); referenceVariant = refVariant; varDescriptionStringToVariants = varDesc; }},
                new Vars(){{ variants.add(goodVariant); referenceVariant = refVariant; varDescriptionStringToVariants = varDesc; }},
                region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t" +
                "7.8\t2.500\t0\t0\t2.0\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "Germline\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    @Test
    public void testGoodVariantsInFirstSample() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{variants.add(goodVariant); referenceVariant = refVariant; }},
                new Vars(){{variants.add(badVariant); referenceVariant = refVariant; }},
                region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8\t" +
                "2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "StrongSomatic\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    @Test
    public void testGoodVariantsInSecondSample() {
        SomaticPostProcessModule module = runMethod();
        module.callingForBothSamples(1,
                new Vars(){{variants.add(badVariant); referenceVariant = refVariant; }},
                new Vars(){{variants.add(goodVariant); referenceVariant = refVariant; }},
                region, new HashSet<>());
        String expectedOutput = "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t0\t3\t5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" +
                "0\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t2.0\t0\t0\t0\t0\t0\t1:1-10\t" +
                "StrongLOH\tSNV\t0\t0\t0\t0\n";
        assertEquals(outContent.toString(), expectedOutput);
    }

    private SomaticPostProcessModule runMethod() {
        VariantPrinter variantPrinter = new SystemOutVariantPrinter();
        variantPrinter.setOut(new PrintStream(outContent));
        SomaticPostProcessModule somaticPostProcessModule = new SomaticPostProcessModule(
                new ReferenceResource(), variantPrinter);
        Configuration config = new Configuration();
        config.goodq = 23;
        GlobalReadOnlyScope.init(config, null, sample, null, "", new HashMap<>(), new HashMap<>());
        return Mockito.spy(somaticPostProcessModule);
    }
}
