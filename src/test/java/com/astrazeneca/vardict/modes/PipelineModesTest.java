package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.integrationtests.IntegrationTest;
import com.astrazeneca.vardict.integrationtests.utils.CSVReferenceManager;
import com.astrazeneca.vardict.printers.SystemOutVariantPrinter;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.*;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

import static com.astrazeneca.vardict.VarDictLauncher.readChr;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.ArgumentMatchers.anyString;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;

public class PipelineModesTest {
    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private Region region = new Region("chr7", 55259400, 55259600, "gene_name");
    private String sample = "test_bam";
    private String bam;
    private VariantPrinter variantPrinter = new SystemOutVariantPrinter();
    private Map<String, Integer> chrLengths;
    private List<List<Region>> segments = new ArrayList<>();
    private List<Region> regions = new ArrayList<>();

    private Scope<InitialData> initialScope;

    @Spy
    private ReferenceResource referenceResource;

    private SimpleMode simpleMode;
    private SomaticMode somaticMode;
    private AmpliconMode ampliconMode;

    @BeforeClass
    public void initResources() {
        CSVReferenceManager.init();
    }

    @BeforeMethod
    public void setUpStreams() throws IOException {
        bam = IntegrationTest.class.getResource("L861Q.bam").getPath();
        chrLengths = readChr(bam);
        GlobalReadOnlyScope.init(new Configuration(), chrLengths, sample, null, null, new HashMap<>(), new HashMap<>());
        instance().conf.bam = new Configuration.BamNames(bam);
        instance().conf.freq = 0.001;
        regions.add(region);
        segments.add(regions);
        System.setErr(new PrintStream(outContent));

        MockitoAnnotations.initMocks(this);
        mockReferenceResource("hg19.fa");
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

    @Test
    public void testSimpleModePipelineNormal() throws ExecutionException, InterruptedException {
        simpleMode = new SimpleMode(segments, referenceResource);
        createScope();
        CompletableFuture<Scope<AlignedVarsData>> pipeline = simpleMode.pipeline(initialScope, new DirectThreadExecutor());

        Scope<AlignedVarsData> varsDataScope = pipeline.get();
        Set<Map.Entry<Integer, Vars>> entryVars = varsDataScope.data.alignedVariants.entrySet();
        List<Vars> actualListOfVars = entryVars.stream().map(Map.Entry::getValue).collect(Collectors.toList());

        assertEquals(actualListOfVars.size(), 1);
        Vars actualVars = actualListOfVars.get(0);
        assertEquals(actualVars.variants.size(), 1);
        assertNotNull(actualVars.referenceVariant);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testSimpleModePipelineException()  {
        simpleMode = new SimpleMode(segments, referenceResource);
        SimpleMode simpleModeSpy = Mockito.spy(simpleMode);
        bam = null;
        createScope();

        Mockito.doNothing().when(simpleModeSpy).stopVardictWithException(any(), any());
        CompletableFuture<Scope<AlignedVarsData>> pipeline = simpleModeSpy.pipeline(initialScope, new DirectThreadExecutor());
        Assert.assertTrue(pipeline.isCompletedExceptionally());
        pipeline.join();
    }

    @Test
    public void testSomaticModePipelineNormal() throws ExecutionException, InterruptedException {
        somaticMode = new SomaticMode(segments, referenceResource);
        createScope();
        CompletableFuture<Scope<AlignedVarsData>> pipeline = somaticMode.pipeline(initialScope, new DirectThreadExecutor());

        Scope<AlignedVarsData> tuple2Scope = pipeline.get();
        Set<Map.Entry<Integer, Vars>> entryVars = tuple2Scope.data.alignedVariants.entrySet();
        List<Vars> actualListOfVars = entryVars.stream().map(Map.Entry::getValue).collect(Collectors.toList());

        assertEquals(actualListOfVars.size(), 1);
        Vars actualVars = actualListOfVars.get(0);
        assertEquals(actualVars.variants.size(), 1);
        assertNotNull(actualVars.referenceVariant);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testSomaticModePipelineException() {
        somaticMode = new SomaticMode(segments, referenceResource);
        SomaticMode somaticModeSpy = Mockito.spy(somaticMode);
        bam = null;
        createScope();

        Mockito.doNothing().when(somaticModeSpy).stopVardictWithException(any(), any());
        CompletableFuture<Scope<AlignedVarsData>> pipeline = somaticModeSpy.pipeline(initialScope, new DirectThreadExecutor());
        Assert.assertTrue(pipeline.isCompletedExceptionally());
        pipeline.join();
    }

    @Test
    public void testAmpliconModePipelineNormal() throws ExecutionException, InterruptedException {
        GlobalReadOnlyScope.clear();
        GlobalReadOnlyScope.init(new Configuration(), chrLengths, sample, null, "100:0.95", new HashMap<>(), new HashMap<>());
        instance().conf.bam = new Configuration.BamNames(bam);
        instance().conf.freq = 0.001;
        ampliconMode = new AmpliconMode(segments, referenceResource);
        createScope();

        CompletableFuture<Scope<AlignedVarsData>> pipeline = ampliconMode.pipeline(initialScope, new DirectThreadExecutor());
        Scope<AlignedVarsData> tuple2Scope = pipeline.get();
        Set<Map.Entry<Integer, Vars>> entryVars = tuple2Scope.data.alignedVariants.entrySet();
        List<Vars> actualListOfVars = entryVars.stream().map(Map.Entry::getValue).collect(Collectors.toList());
        assertEquals(actualListOfVars.size(), 178);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testAmpliconModePipelineException()  {
        GlobalReadOnlyScope.clear();
        GlobalReadOnlyScope.init(new Configuration(), chrLengths, sample, null, "100:0.95", new HashMap<>(), new HashMap<>());
        ampliconMode = new AmpliconMode(segments, referenceResource);
        AmpliconMode ampliconModeSpy = Mockito.spy(ampliconMode);
        createScope();

        Mockito.doNothing().when(ampliconModeSpy).stopVardictWithException(any(), any());
        CompletableFuture<Scope<AlignedVarsData>> pipeline = ampliconModeSpy.pipeline(initialScope, new DirectThreadExecutor());
        Assert.assertTrue(pipeline.isCompletedExceptionally());
        pipeline.join();
    }

    private void mockReferenceResource(String fastaFileName) {
        Mockito.doAnswer(invocation -> {
            Object[] args = invocation.getArguments();
            return CSVReferenceManager.getReader(fastaFileName).queryFasta((String) args[1], (Integer) args[2], (Integer) args[3]);
        }).when(referenceResource).retrieveSubSeq(any(), anyString(), anyInt(), anyInt());
    }

    private void createScope() {
        variantPrinter.setOut(new PrintStream(outContent));
        Reference ref = referenceResource.getReference(region);
        initialScope = new Scope<>(bam, region, ref, referenceResource, 0, new HashSet<>(),
                variantPrinter, new InitialData());
    }
}
