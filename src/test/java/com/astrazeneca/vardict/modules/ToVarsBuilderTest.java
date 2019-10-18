package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.*;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.assertEquals;

public class ToVarsBuilderTest {
    @Spy
    private ToVarsBuilder toVarsBuilder = new ToVarsBuilder();

    @BeforeMethod
    public void setUpStreams() {
        MockitoAnnotations.initMocks(this);
    }

    @AfterMethod
    public void cleanUp () {
        GlobalReadOnlyScope.clear();
    }

    @Test
    public void createInsertion() {
        Configuration config = new Configuration();
        config.goodq = 23;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());

        Map<Integer, VariationMap<String, Variation>> variants = new HashMap<>();
        VariationMap<String, Variation> variantsMap = new VariationMap<>();
        Variation initialVariation = new Variation() {{
            varsCount = 4;
            varsCountOnReverse = 5;
            varsCountOnForward = 3;
            meanPosition = 9;
            meanQuality = 10.5;
            meanMappingQuality = 31;
            numberOfMismatches = 8;
            highQualityReadsCount = 44;
            lowQualityReadsCount = 35;
        }};
        variantsMap.put("T", initialVariation);
        variants.put(1234567, variantsMap);

        List<Variant> var = new ArrayList<>();
        List<String> tmp = new ArrayList<>();

        Mockito.when(toVarsBuilder.getInsertionVariants()).thenReturn(variants);
        int insertionTcov = toVarsBuilder.createInsertion(0.0, 1234567, 10, var, tmp, 0);

        Variant expectedVariant = new Variant() {{
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
            hicov = 44;
            highQualityToLowQualityRatio = 1.2571428571428571;
            highQualityReadsFrequency = 1.0;
        }};

        assertEquals(var.get(0).toString(), expectedVariant.toString());
        assertEquals(insertionTcov, 10);
    }

    @Test
    public void createVariant() {
        Configuration config = new Configuration();
        config.goodq = 23;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());

        Map<Integer, VariationMap<String, Variation>> variants = new HashMap<>();
        VariationMap<String, Variation> variantsMap = new VariationMap<>();
        Map<Integer, Vars> vars = new HashMap<>();

        Variation initialVariation = new Variation() {{
            varsCount = 4;
            varsCountOnReverse = 5;
            varsCountOnForward = 3;
            meanPosition = 9;
            meanQuality = 10.5;
            meanMappingQuality = 31;
            numberOfMismatches = 8;
            highQualityReadsCount = 44;
            lowQualityReadsCount = 35;
        }};
        variantsMap.put("T", initialVariation);
        variants.put(1234567, variantsMap);

        List<Variant> var = new ArrayList<>();
        List<String> tmp = new ArrayList<>();

        Mockito.when(toVarsBuilder.getNonInsertionVariants()).thenReturn(variants);
        toVarsBuilder.createVariant(0.0, vars,1234567,
                toVarsBuilder.getNonInsertionVariants().get(1234567),10, var, tmp,
                new ArrayList<>(toVarsBuilder.getNonInsertionVariants().get(1234567).keySet()),0);

        Variant expectedVariant = new Variant() {{
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

        assertEquals(var.get(0).toString(), expectedVariant.toString());
    }

    @Test
    public void testValidateRefAllele(){
        List<String> alleles = Arrays.asList("A", "C", "G", "T", "N", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B");
        List<String> expected = Arrays.asList("A", "C", "G", "T", "N", "A", "A", "A", "C", "C", "G", "A", "A", "A", "C");
        List<String> validatedAlleles = new ArrayList<>();
        for (String allele: alleles) {
            validatedAlleles.add(toVarsBuilder.validateRefallele(allele));
        }
        assertEquals(validatedAlleles, expected);

        List<String> alleles_complex = Arrays.asList("ANYCGT", "MRACT", "CCGKBG");
        List<String> expected_complex = Arrays.asList("ANCCGT", "AAACT", "CCGGCG");

        validatedAlleles = new ArrayList<>();
        for (String allele: alleles_complex) {
            validatedAlleles.add(toVarsBuilder.validateRefallele(allele));
        }
        assertEquals(validatedAlleles, expected_complex);
    }
}
