package com.astrazeneca.vardict.data.fishertest;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static com.astrazeneca.vardict.Utils.roundHalfEven;

public class FisherExact_Test {
    @Test
    public void test_logdc() {
        FisherExact fisher_test = new FisherExact(11, 12, 1, 2);
        List<Double> logdc = fisher_test.getLogdc();
        List<Double> expected = new ArrayList<Double>() {{
            add(-2.4696392);
            add(-1.0345547);
            add(-0.8675006);
            add(-1.9661129);
        }};
        Assert.assertTrue(logdc.equals(expected));
        // -2.4696392 -1.0345547 -0.8675006 -1.9661129
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
                { new FisherExact(41, 86,	1,	1), 0.54688, 0.54687, 0.89571, 0.47973 },
                { new FisherExact(130, 189,	1,	0), 0.40937, 0.40937, 1.0, 0.0 }, // R: 0.409375 round to 0.40938 in both
                { new FisherExact(83,  40,	1,	2), 0.25746, 0.96473, 0.25746, 4.09908 },
                { new FisherExact(74, 117,	1,	0), 0.39062, 0.39063, 1.0, 0.0 },
                { new FisherExact(60, 62,	2, 0), 0.49593, 0.24797, 1.0, 0.0 },
                { new FisherExact(43, 83,	1, 1), 1.0, 0.57111, 0.88361, 0.52091 },
                { new FisherExact(78,	40,	1,	1), 1.0, 0.88515, 0.56849, 1.93844 },
        };
    }

    @Test(dataProvider = "counts")
    public void test_fexact(FisherExact fisherExact,
                            double pvalue, double pvalueLess, double pValueGreater, double oddRatio) {
        double p = fisherExact.getPValue();
        assert roundHalfEven("0.00000", p) == pvalue;

        p = fisherExact.getPValueLess();
        assert roundHalfEven("0.00000", p)  == pvalueLess;

        p = fisherExact.getPValueGreater();
        assert roundHalfEven("0.00000", p)  == pValueGreater;

        String odd = fisherExact.getOddRatio();
        assert roundHalfEven("0.00000", Double.valueOf(odd)) == oddRatio;
    }
}