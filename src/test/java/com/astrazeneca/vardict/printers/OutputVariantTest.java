package com.astrazeneca.vardict.printers;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.variations.Variant;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class OutputVariantTest {
    private VariantPrinter printer = new SystemOutVariantPrinter();
    private ByteArrayOutputStream outContent;

    @BeforeMethod
    public void setUpStreams() {
        outContent = new ByteArrayOutputStream();
        GlobalReadOnlyScope.init(new Configuration(), null, null, null, "", new HashMap<>(), new HashMap<>());
    }

    @AfterMethod
    public void cleanUpStreams() throws IOException {
        outContent.close();
        GlobalReadOnlyScope.clear();
    }

    @Test
    public void testColumnsSimpleVariant(){
        int position = 1234567;
        Region region = new Region("chr1", 1234560, 1234570, ".");
        OutputVariant simpleOutputVariant = new SimpleOutputVariant(new Variant(), region, "", position);

        printer.setOut(new PrintStream(outContent));
        printer.print(simpleOutputVariant);

        assertEquals(outContent.toString().split("\t").length, 36);
    }

    @Test
    public void testColumnsNullVariant(){
        int position = 1234567;
        Region region = new Region("chr1", 1234560, 1234570, ".");
        OutputVariant simpleNullOutputVariant = new SimpleOutputVariant(null, region, "", position);

        printer.setOut(new PrintStream(outContent));
        printer.print(simpleNullOutputVariant);

        assertEquals(outContent.toString().split("\t").length, 36);
    }

    @Test
    public void testColumnsSomaticVariant(){
        Variant variant = new Variant();
        Region region = new Region("chr1", 1234560, 1234570, ".");
        OutputVariant somaticOutputVariant = new SomaticOutputVariant(variant, variant, variant, variant, region, "", "", "");

        printer.setOut(new PrintStream(outContent));
        printer.print(somaticOutputVariant);

        assertEquals(outContent.toString().split("\t").length, 55);
    }

    @Test
    public void testColumnsNullSomaticVariant(){
        Region region = new Region("chr1", 1234560, 1234570, ".");
        OutputVariant somaticOutputVariant = new SomaticOutputVariant(null, null, null,null, region, "", "", "");

        printer.setOut(new PrintStream(outContent));
        printer.print(somaticOutputVariant);

        assertEquals(outContent.toString().split("\t").length, 55);
    }

    @Test
    public void testColumnsAmpliconVariant() {
        int position = 1234567;
        Region region = new Region("chr1", 1234560, 1234570, ".");
        List<Tuple.Tuple2<Variant, String>> goodVariants = new ArrayList<>();
        List<Tuple.Tuple2<Variant, String>> badVariants = new ArrayList<>();
        OutputVariant ampliconOutputVariant = new AmpliconOutputVariant(new Variant(), region, goodVariants, badVariants, position, 1, 0, false);

        printer.setOut(new PrintStream(outContent));
        printer.print(ampliconOutputVariant);

        assertEquals(outContent.toString().split("\t").length, 38);
    }

    @Test
    public void testColumnsNullAmpliconVariant(){
        int position = 1234567;
        Region region = new Region("chr1", 1234560, 1234570, ".");
        OutputVariant ampliconNullOutputVariant = new AmpliconOutputVariant(null, region, null, null, position, 0,0, false);

        printer.setOut(new PrintStream(outContent));
        printer.print(ampliconNullOutputVariant);

        assertEquals(outContent.toString().split("\t").length, 38);
    }
}
