package com.astrazeneca.vardict.postprocessmodules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.data.scopedata.Scope;
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
import java.util.*;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static org.testng.Assert.assertEquals;

public class SimplePostProcessModuleTest {

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
        varsCountOnReverse = 0;
        varsCountOnForward = 0;
        strandBiasFlag = "0";
        frequency = 0.9;
        meanPosition = 8;
        meanQuality = 26;
        meanMappingQuality = 7.8;
        numberOfMismatches = 2.0;
        hicnt = 44;
        highQualityToLowQualityRatio = 2.5;
    }};

    @Test
    public void testNoVariantNoReferenceVariant() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, region, null,
                null, 0, null, null, tuple);
        SimplePostProcessModule module = runMethod();
        module.accept(newScope);
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testNoVariantButReferenceVariantAndPileup() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            referenceVariant = refVariant;
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, region, null,
                null, 0, null, null, tuple);

        SimplePostProcessModule module = runMethod();
        instance().conf.doPileup = true;
        module.accept(newScope);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tT\t0\t4\t0\t0\t0\t0\tT/T\t0.9000" +
                "\t0\t8.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\t\t0\t0\n");
    }

    @Test
    public void testNoVariantButReferenceVariantWithoutPileup() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            referenceVariant = refVariant;
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, region, null,
                null, 0, null, null, tuple);

        SimplePostProcessModule module = runMethod();
        module.accept(newScope);
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testBadVariant() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            variants.add(badVariant); referenceVariant = refVariant;
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, region, null,
                null, 0, null, null, tuple);

        SimplePostProcessModule module = runMethod();
        module.accept(newScope);
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testGoodVariant() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, region, null,
                null, 0, null, null, tuple);
        SimplePostProcessModule module = runMethod();
        module.accept(newScope);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000" +
                "\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\tSNV\t0\t0\n");
    }

    @Test
    public void testFewGoodVariants() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            variants.add(goodVariant);
            variants.add(goodVariant);
            referenceVariant = refVariant;
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, region, null,
                null, 0, null, null, tuple);
        SimplePostProcessModule module = runMethod();
        module.accept(newScope);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000" +
                "\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\tSNV\t0\t0\n" +
                "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000\t2\t9.0\t0\t26.0\t0\t7.8" +
                "\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\tSNV\t0\t0\n");
    }

    private SimplePostProcessModule runMethod() {
        VariantPrinter variantPrinter = new SystemOutVariantPrinter();
        variantPrinter.setOut(new PrintStream(outContent));
        SimplePostProcessModule simplePostProcessModule = new SimplePostProcessModule(variantPrinter);
        Configuration config = new Configuration();
        config.goodq = 23;
        GlobalReadOnlyScope.init(config, null, sample, null, "", new HashMap<>(), new HashMap<>());
        return Mockito.spy(simplePostProcessModule);
    }
}
