package com.astrazeneca.vardict.integrationtests;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.integrationtests.utils.CSVReferenceManager;
import com.astrazeneca.vardict.integrationtests.utils.OutputFileCreator;
import com.astrazeneca.vardict.integrationtests.utils.VarDictInput;
import com.astrazeneca.vardict.*;
import org.apache.commons.cli.ParseException;
import org.mockito.InjectMocks;
import org.mockito.Mockito;
import org.mockito.MockitoAnnotations;
import org.mockito.Spy;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.*;
import java.util.*;

import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.ArgumentMatchers.anyString;

public class IntegrationTest {

    private final static String DEFAULT_ARGS = "-z -c 1 -S 2 -E 3 -g 4 -G ";
    private final static String PATH_TO_TESTCASES = "testdata/integrationtestcases/";

    private ByteArrayOutputStream outContent;

    @Spy
    private ReferenceResource referenceResource;

    @InjectMocks
    private VarDictLauncher vardictLauncher;

    private static class TestCase {
        VarDictInput input;
        String output;

        @Override
        public String toString() {
            return input.toString();
        }
    }

    @BeforeClass
    public void initResources() {
        CSVReferenceManager.init();
    }

    private void mockReferenceResource(String fastaFileName) throws Exception {
        Mockito.doAnswer(invocation -> {
            Object[] args = invocation.getArguments();
            return CSVReferenceManager.getReader(fastaFileName).queryFasta((String) args[1], (Integer) args[2], (Integer) args[3]);
        }).when(referenceResource).retrieveSubSeq(anyString(), anyString(), anyInt(), anyInt());
    }

    @BeforeMethod
    public void setUpStreams() {
        MockitoAnnotations.initMocks(this);
        outContent = new ByteArrayOutputStream();
        System.setOut(new PrintStream(outContent));
    }

    @AfterMethod
    public void cleanUpStreams() throws IOException {
        GlobalReadOnlyScope.clear();
        System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
        outContent.close();
    }

    @DataProvider(name = "testCases")
    public Object[][] testCases() {
        try {
            return Arrays.stream(new File(PATH_TO_TESTCASES).listFiles())
                    .map(this::parseFile)
                    .map(testCase -> new Object[]{testCase})
                    .toArray(Object[][]::new);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    @Test(dataProvider = "testCases")
    public void integrationTest(Object testCaseObject) throws Exception {
        TestCase testCase = (TestCase) testCaseObject;
        mockReferenceResource(testCase.input.fasta);

        runVarDict(getArgs(testCase.input));

        assertResults(testCase);
    }

    private void assertResults(TestCase testCase) {
        String actualOutput = outContent.toString();
        if(!testCase.output.equals(actualOutput)) {
            File actualJson = OutputFileCreator.fileCreator(testCase.input.toString(), "IntegrationTestActual", actualOutput);
            File expectedJson = OutputFileCreator.fileCreator(testCase.input.toString(), "IntegrationTestExpected", testCase.output);
            Assert.fail("Data aren't identical. Data type: " + testCase.input.toString()
                    + ", IntegrationTest"
                    + "\n Expected: " + expectedJson.getAbsolutePath()
                    + "\n Actual: " + actualJson.getAbsolutePath());
        }
    }

    private String[] getArgs(VarDictInput varDictInput) throws Exception {
        File bed = makeBedFile(varDictInput);
        String[] args;
        String[] bams = varDictInput.bam.split("\\|");
        if (bams.length > 1) {
            args = buildArgs(varDictInput, bed,
                    IntegrationTest.class.getResource(bams[0]).getPath() + "|" +
                            IntegrationTest.class.getResource(bams[1]).getPath());
        } else {
            args = buildArgs(varDictInput, bed, IntegrationTest.class.getResource(varDictInput.bam).getPath());
        }
        return args;
    }

    private TestCase parseFile(File testCase) {
        TestCase testData = new TestCase();
        try(Scanner in = new Scanner(testCase)) {
            testData.input = VarDictInput.fromCSVLine(in.nextLine());
            StringBuilder output = new StringBuilder();
            while (in.hasNext()) {
                output.append(in.nextLine() + "\n");
            }
            testData.output = output.toString();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            throw new UncheckedIOException(e);
        }
        return testData;
    }

    private static String[] buildArgs(VarDictInput varDictInput, File bed, String bam) {
        return (bed
                + " -b " + bam
                + " "
                + (varDictInput.args.equals("") ? "" : (varDictInput.args + " "))
                + DEFAULT_ARGS
                + varDictInput.fasta + ".csv")
                .split("\\s");
    }


    private void runVarDict(String[] args) throws ParseException {
        Configuration conf = new CmdParser().parseParams(args);
        vardictLauncher = new VarDictLauncher(referenceResource);
        vardictLauncher.start(conf);
    }

    public static File makeBedFile(VarDictInput varDictInput) throws IOException {
        File bed = File.createTempFile("tmpbed", ".bed");
        bed.deleteOnExit();

        try (FileWriter bedWriter = new FileWriter(bed)) {
            String bedline = varDictInput.startAmp == null ?
                    varDictInput.chr + "\t" + varDictInput.start + "\t" + varDictInput.end + "\ttestbed" :
                    varDictInput.chr + "\t" + varDictInput.start + "\t" + varDictInput.end + "\ttestbed" + "\t0" + "\t." + "\t" + varDictInput.startAmp + "\t" + varDictInput.endAmp;
            bedWriter.write(bedline);
        }
        return bed;
    }

    @DataProvider(name = "incorrectBedFiles")
    public Object[][] incorrectBedFiles() {
        return new Object[][]{
                {"beds/incorrect_4columns.bed"},
                {"beds/incorrect_8columns.bed"}
        };
    }

    @Test(dataProvider = "incorrectBedFiles", expectedExceptions = NumberFormatException.class)
    public void testIncorrectBedFileRead(String bedFilePath) throws Exception {
        File bam = new File(IntegrationTest.class.getResource("L861Q.bam").getPath());
        File bed = new File(IntegrationTest.class.getResource(bedFilePath).getPath());

        String fastaPath ="hg19.fa";
        mockReferenceResource(fastaPath);

        String[] args = (bed + " -b " + bam + " " + DEFAULT_ARGS + fastaPath + "csv").split("\\s");

        runVarDict(args);
    }

    @DataProvider(name = "correctBedFiles")
    public Object[][] correctBedFiles() {
        return new Object[][]{
                {"beds/correct_5columns.bed"},
                {"beds/correct_9columns.bed"},
                {"beds/correct_8columns_not_amplicon.bed"}
        };
    }

    @Test(dataProvider = "correctBedFiles")
    public void testCorrectBedFileRead(String bedFilePath) throws Exception {
        File bam = new File(IntegrationTest.class.getResource("L861Q.bam").getPath());
        File bed = new File(IntegrationTest.class.getResource(bedFilePath).getPath());
        StringBuilder output = new StringBuilder();
        try(Scanner in = new Scanner(new File(PATH_TO_TESTCASES + "Simple;hg19.fa;L861Q.bam;chr7;55259400-55259600;.txt"))) {
            while (in.hasNext()) {
                output.append(in.nextLine() + "\n");
            }
        }
        String fastaPath ="hg19.fa";
        mockReferenceResource(fastaPath);
        String[] args = (bed + " -b " + bam + " " + DEFAULT_ARGS + fastaPath + ".csv").split("\\s");

        runVarDict(args);

        String actualOutput = "Simple,hg19.fa,L861Q.bam,chr7,55259400,55259600\n" + outContent.toString();
        Assert.assertEquals(actualOutput, output.toString());
    }

    @Test
    public void testCorrectAmpliconBedFileRead() throws Exception {
        File bam = new File(IntegrationTest.class.getResource("L861Q.bam").getPath());
        File bed = new File(IntegrationTest.class.getResource("beds/correct_8columns_amplicon.bed").getPath());
        String fastaPath ="hg19.fa";
        mockReferenceResource(fastaPath);
        String[] args = (bed + " -b " + bam + " " + DEFAULT_ARGS + fastaPath + ".csv").split("\\s");
        runVarDict(args);
    }
}

