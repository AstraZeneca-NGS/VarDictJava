package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.integrationtests.IntegrationTest;
import com.astrazeneca.vardict.integrationtests.utils.CSVReferenceManager;
import com.astrazeneca.vardict.printers.SystemOutVariantPrinter;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;
import org.mockito.*;
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
import static org.testng.Assert.assertTrue;

public class PipelineModesTest {
    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private Region region = new Region("chr7", 55259400, 55259600, "gene_name");
    private String sample = "test_bam";
    private String bam;
    private VariantPrinter variantPrinter = new SystemOutVariantPrinter();
    private Map<String, Integer> chrLengths;

    @Spy
    private ReferenceResource referenceResource;

    @InjectMocks
    SimpleMode simpleMode;

    @InjectMocks
    SomaticMode somaticMode;

    @InjectMocks
    AmpliconMode ampliconMode;

    @BeforeClass
    public void initResources() {
        CSVReferenceManager.init();
    }

    @BeforeMethod
    public void setUpStreams() throws IOException {
        bam = IntegrationTest.class.getResource("L861Q.bam").getPath();
        chrLengths = readChr(bam);
        GlobalReadOnlyScope.init(new Configuration(), chrLengths, sample, null, null);
        instance().conf.bam = new Configuration.BamNames(bam);
        instance().conf.freq = 0.001;
        MockitoAnnotations.initMocks(this);
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
        mockReferenceResource("hg19.fa");
        List<List<Region>> segments = new ArrayList<>();
        List<Region> regions = new ArrayList<>();
        regions.add(region);
        segments.add(regions);

        simpleMode = new SimpleMode(segments, referenceResource);
        Reference ref = referenceResource.getReference(region);
        variantPrinter.setOut(new PrintStream(outContent));

        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = simpleMode.pipeline(bam, region, ref, referenceResource, 0,
                new HashSet<>(), variantPrinter, new DirectThreadExecutor());

        Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> tuple2Scope = pipeline.get();
        Set<Map.Entry<Integer, Vars>> entryVars = tuple2Scope.data._2.entrySet();
        List<Vars> actualListOfVars = entryVars.stream().map(Map.Entry::getValue).collect(Collectors.toList());

        assertEquals(actualListOfVars.size(), 1);
        Vars actualVars = actualListOfVars.get(0);
        assertEquals(actualVars.variants.size(), 1);
        assertNotNull(actualVars.referenceVariant);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testSimpleModePipelineException()  {
        mockReferenceResource("hg19.fa");
        List<List<Region>> segments = new ArrayList<>();
        List<Region> regions = new ArrayList<>();
        regions.add(region);
        segments.add(regions);

        simpleMode = new SimpleMode(segments, referenceResource);
        Reference ref = referenceResource.getReference(region);
        variantPrinter.setOut(new PrintStream(outContent));
        System.setErr(new PrintStream(outContent));
        instance().conf.bam = null;
        Configuration.MAX_EXCEPTION_COUNT = 5;

        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = simpleMode.pipeline(bam, region, ref, referenceResource, 0,
                new HashSet<>(), variantPrinter, new DirectThreadExecutor());
        pipeline.join();
        assertTrue(outContent.toString().contains("NullPointerException"));
    }

    @Test
    public void testSomaticModePipelineNormal() throws ExecutionException, InterruptedException {
        mockReferenceResource("hg19.fa");
        List<List<Region>> segments = new ArrayList<>();
        List<Region> regions = new ArrayList<>();
        regions.add(region);
        segments.add(regions);

        somaticMode = new SomaticMode(segments, referenceResource);
        Reference ref = referenceResource.getReference(region);
        variantPrinter.setOut(new PrintStream(outContent));

        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = somaticMode.pipeline(bam, region, ref, referenceResource, 0,
                new HashSet<>(), variantPrinter, new DirectThreadExecutor());

        Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> tuple2Scope = pipeline.get();
        Set<Map.Entry<Integer, Vars>> entryVars = tuple2Scope.data._2.entrySet();
        List<Vars> actualListOfVars = entryVars.stream().map(Map.Entry::getValue).collect(Collectors.toList());

        assertEquals(actualListOfVars.size(), 1);
        Vars actualVars = actualListOfVars.get(0);
        assertEquals(actualVars.variants.size(), 1);
        assertNotNull(actualVars.referenceVariant);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testSomaticModePipelineException() {
        mockReferenceResource("hg19.fa");
        List<List<Region>> segments = new ArrayList<>();
        List<Region> regions = new ArrayList<>();
        regions.add(region);
        segments.add(regions);

        somaticMode = new SomaticMode(segments, referenceResource);
        Reference ref = referenceResource.getReference(region);
        variantPrinter.setOut(new PrintStream(outContent));
        System.setErr(new PrintStream(outContent));
        instance().conf.bam = null;
        Configuration.MAX_EXCEPTION_COUNT = 5;

        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = somaticMode.pipeline(bam, region, ref, referenceResource, 0,
                new HashSet<>(), variantPrinter, new DirectThreadExecutor());
        pipeline.join();
        assertTrue(outContent.toString().contains("NullPointerException"));
    }

    @Test
    public void testAmpliconModePipelineNormal() throws ExecutionException, InterruptedException {
        GlobalReadOnlyScope.clear();
        GlobalReadOnlyScope.init(new Configuration(), chrLengths, sample, null, "100:0.95");
        instance().conf.bam = new Configuration.BamNames(bam);
        instance().conf.freq = 0.001;

        mockReferenceResource("hg19.fa");
        List<List<Region>> segments = new ArrayList<>();
        List<Region> regions = new ArrayList<>();
        regions.add(region);
        segments.add(regions);

        ampliconMode = new AmpliconMode(segments, referenceResource);
        Reference ref = referenceResource.getReference(region);
        variantPrinter.setOut(new PrintStream(outContent));

        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = ampliconMode.pipeline(bam, region, ref, referenceResource, 0,
                new HashSet<>(), variantPrinter, new DirectThreadExecutor());

        Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> tuple2Scope = pipeline.get();
        Set<Map.Entry<Integer, Vars>> entryVars = tuple2Scope.data._2.entrySet();
        List<Vars> actualListOfVars = entryVars.stream().map(Map.Entry::getValue).collect(Collectors.toList());

        assertEquals(actualListOfVars.size(), 178);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testAmpliconModePipelineException()  {
        GlobalReadOnlyScope.clear();
        GlobalReadOnlyScope.init(new Configuration(), chrLengths, sample, null, "100:0.95");
        Configuration.MAX_EXCEPTION_COUNT = 5;

        mockReferenceResource("hg19.fa");
        List<List<Region>> segments = new ArrayList<>();
        List<Region> regions = new ArrayList<>();
        regions.add(region);
        segments.add(regions);

        ampliconMode = new AmpliconMode(segments, referenceResource);
        Reference ref = referenceResource.getReference(region);
        variantPrinter.setOut(new PrintStream(outContent));
        System.setErr(new PrintStream(outContent));

        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = ampliconMode.pipeline(bam, region, ref, referenceResource, 0,
                new HashSet<>(), variantPrinter, new DirectThreadExecutor());

        pipeline.join();
        assertTrue(outContent.toString().contains("NullPointerException"));
    }

    private void mockReferenceResource(String fastaFileName) {
        Mockito.doAnswer(invocation -> {
            Object[] args = invocation.getArguments();
            return CSVReferenceManager.getReader(fastaFileName).queryFasta((String) args[1], (Integer) args[2], (Integer) args[3]);
        }).when(referenceResource).retrieveSubSeq(any(), anyString(), anyInt(), anyInt());
    }
}
