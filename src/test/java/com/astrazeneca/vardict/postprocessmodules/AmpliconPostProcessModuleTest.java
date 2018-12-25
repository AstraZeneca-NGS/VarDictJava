package com.astrazeneca.vardict.postprocessmodules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.printers.SystemOutVariantPrinter;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.Mockito;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import static org.testng.Assert.assertEquals;

public class AmpliconPostProcessModuleTest {

    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private Region region = new Region("1", 1, 10, "gene_name");
    private String sample = "test_bam";
    private VariantPrinter variantPrinter = new SystemOutVariantPrinter();
    private Map<Integer, List<Tuple.Tuple2<Integer, Region>>> positionsMap = new HashMap<>();
    List<Tuple.Tuple2<Integer, Region>> positions = new ArrayList<>();

    @BeforeMethod
    public void initStructures() {
        positions.add(new Tuple.Tuple2<>(0,region));
        positionsMap.put(1, positions);
    }

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

    Variant anotherGoodVariant = new Variant() {{
        refallele = "T";
        varallele = "A";
        descriptionString = "A";
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
    public void testGoodVariantsOnBothAmps() {
        Map<Integer, Vars> variantMap1Amp = new HashMap<>();
        variantMap1Amp.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        Map<Integer, Vars> variantMap2Amp = new HashMap<>();
        variantMap2Amp.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        List<Map<Integer, Vars>> mapsOnAmplicon = new ArrayList<>();
        mapsOnAmplicon.add(variantMap1Amp);
        mapsOnAmplicon.add(variantMap2Amp);

        AmpliconPostProcessModule module = runMethod();
        module.process(region, mapsOnAmplicon, positionsMap, new HashSet<>(), variantPrinter);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000" +
                "\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\tSNV\t2\t2\t0\t0\n");
    }

    @Test
    public void testBadVariantsOnBoth() {
        Map<Integer, Vars> variantMap1Amp = new HashMap<>();
        variantMap1Amp.put(1, new Vars() {{
            variants.add(badVariant); referenceVariant = refVariant;
        }});
        Map<Integer, Vars> variantMap2Amp = new HashMap<>();
        variantMap2Amp.put(1, new Vars() {{
            variants.add(badVariant); referenceVariant = refVariant;
        }});
        List<Map<Integer, Vars>> mapsOnAmplicon = new ArrayList<>();
        mapsOnAmplicon.add(variantMap1Amp);
        mapsOnAmplicon.add(variantMap2Amp);

        AmpliconPostProcessModule module = runMethod();
        module.process(region, mapsOnAmplicon, positionsMap, new HashSet<>(), variantPrinter);
        assertEquals(outContent.toString(), "");
    }

    @Test
    public void testGoodVariantsOnTwoAmps() {
        Map<Integer, Vars> variantMap1Amp = new HashMap<>();
        variantMap1Amp.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        Map<Integer, Vars> variantMap2Amp = new HashMap<>();
        variantMap2Amp.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        List<Map<Integer, Vars>> mapsOnAmplicon = new ArrayList<>();
        mapsOnAmplicon.add(variantMap1Amp);
        mapsOnAmplicon.add(variantMap2Amp);
        positions.add(new Tuple.Tuple2<>(1,region));
        positionsMap.put(1,positions);

        AmpliconPostProcessModule module = runMethod();
        module.process(region, mapsOnAmplicon, positionsMap, new HashSet<>(), variantPrinter);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000" +
                "\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\tSNV\t2\t2\t0\t0\n");

    }

    @Test
    public void testGoodVariantOnOneAmpBadOnAnother() {
        Map<Integer, Vars> variantMap1Amp = new HashMap<>();
        variantMap1Amp.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        Map<Integer, Vars> variantMap2Amp = new HashMap<>();
        variantMap2Amp.put(1, new Vars() {{
            variants.add(badVariant); referenceVariant = refVariant;
        }});
        List<Map<Integer, Vars>> mapsOnAmplicon = new ArrayList<>();
        mapsOnAmplicon.add(variantMap1Amp);
        mapsOnAmplicon.add(variantMap2Amp);

        AmpliconPostProcessModule module = runMethod();
        module.process(region, mapsOnAmplicon, positionsMap, new HashSet<>(), variantPrinter);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000" +
                "\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-10\tSNV\t1\t2\t0\t0\n");
    }

    @Test
    public void testDifferentGoodVariantsOnDifferentAmps() {
        Map<Integer, Vars> variantMap1Amp = new HashMap<>();
        variantMap1Amp.put(1, new Vars() {{
            variants.add(goodVariant); referenceVariant = refVariant;
        }});
        Map<Integer, Vars> variantMap2Amp = new HashMap<>();

        variantMap2Amp.put(1, new Vars() {{
            variants.add(anotherGoodVariant); referenceVariant = refVariant;
        }});
        List<Map<Integer, Vars>> mapsOnAmplicon = new ArrayList<>();
        mapsOnAmplicon.add(variantMap1Amp);
        mapsOnAmplicon.add(variantMap2Amp);

        positions.add(new Tuple.Tuple2<>(1,region));
        positionsMap.put(1, positions);

        AmpliconPostProcessModule module = runMethod();
        module.process(region, mapsOnAmplicon, positionsMap, new HashSet<>(), variantPrinter);
        assertEquals(outContent.toString(), "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000" +
                "\t2\t9.0\t0\t26.0\t0\t7.8\t2.500\t0\t0\t0\t0\t0\t2.0\t44\t0\t0\t0\t1:1-1\tSNV\t2\t2\t0\t1\n");
    }

    private AmpliconPostProcessModule runMethod() {
        variantPrinter.setOut(new PrintStream(outContent));
        AmpliconPostProcessModule ampliconPostProcessModule = new AmpliconPostProcessModule();

        Configuration config = new Configuration();
        config.goodq = 23;
        GlobalReadOnlyScope.init(config, null, sample, null, "", new HashMap<>(), new HashMap<>());
        return Mockito.spy(ampliconPostProcessModule);
    }
}
