package com.astrazeneca.vardict.integrationtests;

import com.astrazeneca.vardict.integrationtests.utils.CSVReferenceManager;
import com.astrazeneca.vardict.integrationtests.utils.OutputFileCreator;
import com.astrazeneca.vardict.integrationtests.utils.VarDictInput;
import com.astrazeneca.vardict.*;
import org.mockito.invocation.InvocationOnMock;
import org.powermock.api.mockito.PowerMockito;
import org.powermock.core.classloader.annotations.PrepareForTest;
import org.powermock.modules.testng.PowerMockObjectFactory;
import org.testng.Assert;
import org.testng.IObjectFactory;
import org.testng.annotations.*;

import java.io.*;
import java.util.*;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.anyInt;
import static org.mockito.Matchers.anyString;

@PrepareForTest(ReferenceResource.class)
public class IntegrationTest {

    private final static String DEFAULT_ARGS = "-z -c 1 -S 2 -E 3 -g 4 -G ";
    private final static String PATH_TO_TESTCASES = "testdata/integrationtestcases/";

    private ByteArrayOutputStream outContent;

    private static class TestCase {
        VarDictInput input;
        String output;

        @Override
        public String toString() {
            return input.toString();
        }
    }

    @ObjectFactory
    public IObjectFactory getObjectFactory() {
        return new PowerMockObjectFactory();
    }

    @BeforeClass
    public void initResources() {
        CSVReferenceManager.init();
    }

    static void mockReferenceResource(String fastaFileName) throws Exception {
        PowerMockito.mockStatic(ReferenceResource.class);

        PowerMockito.when(ReferenceResource.retriveSubSeq(anyString(), anyString(), anyInt(), anyInt()))
                .thenAnswer(invocation -> {
                    Object[] args = invocation.getArguments();
                    return CSVReferenceManager.getReader(fastaFileName).queryFasta((String) args[1], (Integer) args[2], (Integer) args[3]);
                });

        PowerMockito.when(ReferenceResource.getREF(any(), any(), anyString(), anyInt())).thenAnswer(InvocationOnMock::callRealMethod);
    }

    @BeforeMethod
    public void setUpStreams() {
        outContent = new ByteArrayOutputStream();
        System.setOut(new PrintStream(outContent));
    }

    @AfterMethod
    public void cleanUpStreams() throws IOException {
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

        Main.main(getArgs(testCase.input));
        assertResults(testCase);
    }

    private void assertResults(TestCase testCase) throws IOException {
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
}

