package com.astrazeneca.vardict;

import com.astrazeneca.vardict.data.Region;
import org.mockito.Mockito;
import org.mockito.MockitoAnnotations;
import org.mockito.Spy;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import static java.util.Collections.singletonList;
import static org.testng.Assert.assertEquals;

public class RegionBuilderTest {
    private RegionBuilder regionBuilder;

    private Configuration defaultTestConfig =  getTestConfig(new RegionBuilder.BedRowFormat(0, 1, 2, 1, 2, 3), 0);

    private Configuration getTestConfig(RegionBuilder.BedRowFormat bedRowFormat, int chrColumn) {
        Configuration conf = new Configuration();
        conf.bedRowFormat = bedRowFormat;
        conf.zeroBased = true;
        conf.columnForChromosome = chrColumn;
        conf.numberNucleotideToExtend = 0;
        conf.delimiter = "\t";
        return conf;
    }

    @BeforeMethod
    public void setUp() {
        MockitoAnnotations.initMocks(this);
    }

    private void mockCorrectChromosome(RegionBuilder spy) {
        Mockito.doReturn("Y").when(spy).correctChromosome(Mockito.any(), Mockito.anyString());
    }

    @DataProvider(name = "nonAmp")
    public Object[][] nonAmp() {
        List<String> raws_1 = new LinkedList<String>(){{
            add("Y\t59023460\t59028460\tY_59023460_1");
            add("Y\t59028310\t59033310\tY_59028310_2");
            add("Y\t59033160\t59038160\tY_59033160_3");
            add("Y\t59038010\t59043010\tY_59038010_4");
            add("Y\t59042860\t59047860\tY_59042860_5");
        }};
        List<List<Region>> result_1 = new LinkedList<List<Region>>() {{
            add(singletonList(new Region("Y", 59023461, 59028460, "Y_59023460_1")));
            add(singletonList(new Region("Y", 59028311, 59033310, "Y_59028310_2")));
            add(singletonList(new Region("Y", 59033161, 59038160, "Y_59033160_3")));
            add(singletonList(new Region("Y", 59038011, 59043010, "Y_59038010_4")));
            add(singletonList(new Region("Y", 59042861, 59047860, "Y_59042860_5")));
        }};

        List<String> raws_2 = new LinkedList<String>(){{
            add("Y\t59023460\t59028460\t59028460");
            add("Y\t59028310\t59033310\t59033310");
            add("Y\t59033160\t59038160\t59038160");
            add("Y\t59038010\t59043010\t59043010");
            add("Y\t59042860\t59047860\t59047860");
        }};
        List<List<Region>> result_2 = new LinkedList<List<Region>>() {{
            add(new LinkedList<>());
            add(new LinkedList<>());
            add(new LinkedList<>());
            add(new LinkedList<>());
            add(new LinkedList<>());
        }};
        Configuration conf = getTestConfig(new RegionBuilder.BedRowFormat(0, 1, 2, 1, 2, 3), -1);
        return new Object[][] {
                {raws_1, defaultTestConfig, result_1},
                {raws_2, conf, result_2}
        };
    }

    @Test(dataProvider = "nonAmp")
    public void testNonAmpToRegions(List<String> raws, Object config, List<List<Region>> expected){
        regionBuilder = new RegionBuilder(new HashMap<>(), (Configuration) config);
        RegionBuilder spy = Mockito.spy(regionBuilder);
        mockCorrectChromosome(spy);
        assertEquals(spy.buildRegions(raws, true), expected);
    }

    @DataProvider(name = "singleRegion")
    public Object[][] singleRegion() {
        return new Object[][] {
                {
                        "Y:59023460-59028460:Y_59023460_1",
                        singletonList(singletonList(new Region("Y", 59023461, 59028460, "Y_59023460_1")))
                },
                {
                        "Y:59028310-59033310:Y_59028310_2",
                        singletonList(singletonList(new Region("Y", 59028311, 59033310, "Y_59028310_2")))
                },
                {
                        "Y:59033160-59038160:Y_59033160_3",
                        singletonList(singletonList(new Region("Y", 59033161, 59038160, "Y_59033160_3")))
                },
                {
                        "Y:59038010-59043010:Y_59038010_4",
                        singletonList(singletonList(new Region("Y", 59038011, 59043010, "Y_59038010_4")))
                },
                {
                        "Y:59042860-59047860:Y_59042860_5",
                        singletonList(singletonList(new Region("Y", 59042861, 59047860, "Y_59042860_5")))
                }
        };
    }

    @Test(dataProvider = "singleRegion")
    public void testSingleRegion(String region, Object expected) {
        defaultTestConfig.regionOfInterest = region;
        regionBuilder = new RegionBuilder(new HashMap<>(), defaultTestConfig);
        RegionBuilder spy = Mockito.spy(regionBuilder);
        mockCorrectChromosome(spy);
        assertEquals(spy.buildRegionFromConfiguration(), expected);
    }

    @DataProvider(name = "amp")
    public Object[][] amp() {
        List<String> raws = new LinkedList<String>(){{
            add("Y\t59023460\t59028460\tY_59023460_1\t0\t.\t59023460\t59028460");
            add("Y\t59028310\t59033310\tY_59028310_2\t0\t.\t59028310\t59033310");
            add("Y\t59033160\t59038160\tY_59033160_3\t0\t.\t59033160\t59038160");
            add("Y\t59038010\t59043010\tY_59038010_4\t0\t.\t59038010\t59043010");
            add("Y\t59042860\t59047860\tY_59042860_5\t0\t.\t59042860\t59047860");
        }};
        List<List<Region>> result = new LinkedList<List<Region>>() {{
            add(new LinkedList<Region>() {{
                add(new Region("Y", 59023461, 59028460, "Y_59023460_1", 59023461, 59028460));
                add(new Region("Y", 59028311, 59033310, "Y_59028310_2", 59028311, 59033310));
                add(new Region("Y", 59033161, 59038160, "Y_59033160_3", 59033161, 59038160));
                add(new Region("Y", 59038011, 59043010, "Y_59038010_4", 59038011, 59043010));
                add(new Region("Y", 59042861, 59047860, "Y_59042860_5", 59042861, 59047860));
            }});
        }};
        return new Object[][] {
                {raws, result}
        };
    }

    @Test(dataProvider = "amp")
    public void testAmp(List<String> raws, List<List<Region>> expected) {
        regionBuilder = new RegionBuilder(new HashMap<>(), defaultTestConfig);
        RegionBuilder spy = Mockito.spy(regionBuilder);
        mockCorrectChromosome(spy);
        assertEquals(spy.buildAmpRegions(raws, true), expected);
    }
}
