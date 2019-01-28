package com.astrazeneca.vardict;

import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.postprocessmodules.SimplePostProcessModule;
import com.astrazeneca.vardict.printers.SystemOutVariantPrinter;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.Mockito;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import static org.testng.Assert.assertTrue;

public class ExceptionCounterTest {
    private ByteArrayOutputStream outContent;

    @BeforeMethod
    public void setUpStreams() {
        outContent = new ByteArrayOutputStream();
    }

    @AfterMethod
    public void cleanUpStreams() throws IOException {
        outContent.reset();
        GlobalReadOnlyScope.clear();
    }

    @AfterTest
    public void closeOutContent() throws IOException {
        outContent.close();
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testExceptionMoreThenLimitMustStopWorkAndThrowException() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            variants.add(null);
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, null, null,
                null, 0, null, null, tuple);
        SimplePostProcessModule module = runMethod();
        Configuration.MAX_EXCEPTION_COUNT = 1;
        System.setErr(new PrintStream(outContent));
        module.accept(newScope);
        module.accept(newScope);
    }

    @Test
    public void testExceptionLessThenLimitMustContinueWork() {
        Map<Integer, Vars> variantMap = new HashMap<>();
        variantMap.put(1, new Vars() {{
            variants.add(null);
        }});
        AlignedVarsData tuple = new AlignedVarsData(1, variantMap);
        Scope<AlignedVarsData> newScope = new Scope<>(null, null, null,
                null, 0, null, null, tuple);
        SimplePostProcessModule module = runMethod();
        Configuration.MAX_EXCEPTION_COUNT = 1;
        System.setErr(new PrintStream(outContent));
        module.accept(newScope);
        assertTrue(outContent.toString().contains("NullPointerException"));
    }

    private SimplePostProcessModule runMethod() {
        VariantPrinter variantPrinter = new SystemOutVariantPrinter();
        variantPrinter.setOut(new PrintStream(outContent));
        SimplePostProcessModule simplePostProcessModule = new SimplePostProcessModule(variantPrinter);
        Configuration config = new Configuration();
        config.goodq = 23;
        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());
        return Mockito.spy(simplePostProcessModule);
    }
}
