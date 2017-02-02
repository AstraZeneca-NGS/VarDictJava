package com.astrazeneca.vardict;


import com.astrazeneca.vardict.pipeline.postprocessmodules.SomaticPostProcessModule;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.VariationUtils;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.Mockito;
import org.powermock.api.mockito.PowerMockito;
import org.powermock.core.classloader.annotations.PrepareForTest;
import org.powermock.modules.testng.PowerMockObjectFactory;
import org.testng.IObjectFactory;
import org.testng.annotations.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.powermock.api.mockito.PowerMockito.doReturn;
import static org.powermock.api.mockito.PowerMockito.mockStatic;
import static org.powermock.api.mockito.PowerMockito.when;
import static org.testng.Assert.assertEquals;

@PrepareForTest({SomaticPostProcessModule.class, VariationUtils.class})
public class SomaticPostProcessModuleTest {

    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    Region region = new Region("1", 1, 10, "gene_name");
    String sample = "test_bam";

    @AfterMethod
    public void clearOutContent() {
        outContent.reset();
    }

    @AfterTest
    public void closeOutContent() throws IOException {
        outContent.close();
    }

    @ObjectFactory
    public IObjectFactory setObjectFactory() {
        return new PowerMockObjectFactory();
    }

    @Test
    public void testNoCoverageForBothSamples() throws IOException {
        runMethod(new HashMap<Integer, Vars>(){{
            put(1, null);
        }}, new HashMap<Integer, Vars>(){{
            put(1, null);
        }});
        assertEquals(outContent.toString(), "");
    }

    private void mockGetType(String type) {
        PowerMockito.stub(PowerMockito.method(SomaticPostProcessModule.class, "determinateType")).toReturn(type);
    }

    private void mockGetVarMaybe(Variant result) {
        when(VariationUtils.getVarMaybe(
                any(Vars.class),
                eq(VariationUtils.VarsType.varn),
                any(String.class)
        )).thenReturn(result);
    }

    private void checkPrintResult(String type) {
        String[] content = outContent.toString().split("\t");
        int length = content.length;
        assertEquals(sample, content[0]);
        assertEquals(region.gene, content[1]);
        assertEquals(type, content[length-2]);
    }

    private Map<Integer, Vars> fillTestMap(int count) {
        final Map<Integer, Vars> testVars = new HashMap<>();
        for (int i = 0; i < count; i++) {
            testVars.put(i, new Vars(){{
                referenceVariant = new Variant();
                Variant spy = Mockito.spy(new Variant());
                mockVariant(spy);
                variants.add(spy);
            }});
        }
        return testVars;
    }

    private void mockVariant(Variant v) {
        doReturn(Variant.Type.noInfo).when(v).getType();
    }

    private void runMethod(Map<Integer, Vars> testVars1, Map<Integer, Vars> testVars2) throws IOException {
        new SomaticPostProcessModule();
    }

    private void mockIsGoodVar() {
        mockStatic(VariationUtils.class);
        when(VariationUtils.isGoodVar(
                any(Variant.class),
                any(Variant.class),
                any(Variant.Type.class),
                any(Set.class)
        )).thenReturn(true);
    }
}
