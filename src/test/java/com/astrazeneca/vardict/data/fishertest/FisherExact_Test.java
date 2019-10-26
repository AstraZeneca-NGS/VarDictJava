package com.astrazeneca.vardict.data.fishertest;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class FisherExact_Test {
    @Test
    public void test_logdc() {
        FisherExact fisher_test = new FisherExact(11, 12, 1, 2);
        List<Double> logdc = fisher_test.getLogdc();
        List<Double> expected = new ArrayList<Double>() {{
            add(-2.469639177657212);
            add(-1.0345546523678895);
            add(-0.867500567704723);
            add(-1.9661128563728325);
        }};
        Assert.assertTrue(logdc.equals(expected));
        // -2.4696392 -1.0345547 -0.8675006 -1.966112
    }

    @DataProvider(name = "counts")
    public Object[][] counts() {
        return new Object[][]{
                { new FisherExact(121, 55, 18, 23), 0.00378, 0.99908, 0.00287, 2.79657 },
                { new FisherExact(121, 5, 18, 23), 0.0, 1.0, 0.0, 29.86184 },
                { new FisherExact(37, 76, 1, 1), 1.0, 0.55362, 0.89275, 0.49015 },
                { new FisherExact(0, 0, 1, 0), 1.0, 1.0, 1.0, 0.0 },
                { new FisherExact(1, 0, 1, 0), 1.0, 1.0, 1.0, 0.0 },
                { new FisherExact(0, 0, 0, 0), 1.0, 1.0, 1.0, 0.0 },
                { new FisherExact(0, 1, 1, 0), 1.0, 0.5, 1.0, 0.0 },
                { new FisherExact(1, 1, 1, 0), 1.0, 0.66667, 1.0, 0.0 },
                { new FisherExact(1, 1, 1, 1), 1.0, 0.83333, 0.83333, 1.0 },
                { new FisherExact(10, 10, 10, 1), 0.04722, 0.02599, 0.99802, 0.10703 },
                { new FisherExact(10, 10, 10, 0), 0.01099, 0.00615, 1.0, 0.0 },
                { new FisherExact(69,1,	74,	95), 0.0, 1.0, 0.0, 87.68597 },
        };
    }

    @Test(dataProvider = "counts")
    public void test_fexact(FisherExact fisherExact,
                            double pvalue, double pvalueLess, double pValueGreater, double oddRatio) {
        double p = fisherExact.getPValue();
        System.err.println(p);
        assert p == pvalue;

        p = fisherExact.getPValueLess();
        System.err.println(p);
        assert p == pvalueLess;

        p = fisherExact.getPValueGreater();
        System.err.println(p);
        assert p == pValueGreater;

        String odd = fisherExact.getOddRatio();
        System.err.println(odd);
        assert Double.valueOf(odd) == oddRatio;
    }
}