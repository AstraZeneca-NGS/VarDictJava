package com.astrazeneca.vardict.data.fishertest;

import org.apache.commons.math3.distribution.HypergeometricDistribution;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import static com.astrazeneca.vardict.Utils.roundHalfEven;

/**
 * EXPERIMENTAL FEATURE.
 * <p>
 * Implementation of FisherExact Test as it is implemented in R.
 * <p>
 * R implementation of Fisher Test for oddratio uses conditional MLE (maximum likelihood estimation)
 * that wasn't found in standard libraries for Java.
 * <p>
 * Reason to replace R fisher test with this implementation is a slow R `textConnection` function.
 * In other case we have to use temp files to process VarDict result in R faster, and this is not a good option.
 */

public class FisherExact {
    private List<Double> logdc;
    private int m;
    private int n;
    private int k;
    private int x;
    private int lo;
    private int hi;
    private double PvalueLess;
    private double PvalueGreater;
    private double PvalueTwoSided;
    private List<Integer> support;

    // Seems that Java and R have differences with round half even (JDK-8227248 example, it will round value in memory)
    public static double RESULT_ROUND_R = 1E5;

    public FisherExact(int refFwd, int refRev, int altFwd, int altRev) {
        m = refFwd + refRev;
        n = altFwd + altRev;
        k = refFwd + altFwd;
        x = refFwd;
        lo = Math.max(0, k - n);
        hi = Math.min(k, m);
        support = new ArrayList<>();
        for (int j = lo; j <= hi; j++) {
            support.add(j);
        }
        logdc = logdcDhyper(m, n, k);

        calculatePValue();
    }

    // Density of the central hypergeometric distribution on its support: store for once as this is needed quite a bit.
    private List<Double> logdcDhyper(int m, int n, int k) {
        List<Double> logdc = new ArrayList<>();

        for (int element : support) {
            if (m + n == 0) {
                logdc.add(0.0);
                continue;
            }
            // m + n - total number of successes, m - number of successes (reference) k - sample size (forward)
            HypergeometricDistribution dhyper = new HypergeometricDistribution(m + n, m, k);
            Double value = dhyper.logProbability(element);
            if (value.isNaN()) {
                value = 0.0;
            }
            logdc.add(roundHalfEven("0.0000000", value));
        }
        return logdc;
    }

    // Determine the MLE for ncp by solving E(X) = x, where the expectation is with respect to H.
    // Note that in general the conditional distribution of x given the marginals is a non-central hypergeometric
    // distribution H with non-centrality parameter ncp, the odds ratio.
    // The null conditional independence is equivalent to the hypothesis that the odds ratio equals one. `Exact`
    // inference can be based on observing that in general, given all marginal totals fixed, the first element of the
    // contingency table has a non-central hypergeometric distribution with non-centrality parameter given by odds
    // ratio (Fisher, 1935). The alternative for a one-sided test is based on the odds ratio, so alternative =
    // 'greater' is a test of the odds ratio being bigger than or = 1.
    private Double mle(double x) {
        double eps = Math.ulp(1.0);
        if (x == lo) return 0.0;
        if (x == hi) return Double.POSITIVE_INFINITY;
        double mu = mnhyper(1.0);
        double root;
        if (mu > x) {
            Function<Double, Double> f = t -> mnhyper(t) - x;
            root = UnirootZeroIn.zeroinC(0, 1, f, Math.pow(eps, 0.25));
        } else if (mu < x) {
            Function<Double, Double> f = t -> mnhyper(1.0 / t) - x;
            root = 1.0 / UnirootZeroIn.zeroinC(eps, 1, f, Math.pow(eps, 0.25));
        } else {
            root = 1.0;
        }
        return root;
    }

    private Double mnhyper(Double ncp) {
        if (ncp == 0) return (double) lo;
        if (ncp.isInfinite()) return (double) hi;
        else {
            List<Double> dnhyperResult = dnhyper(ncp);
            List<Double> multiply = new ArrayList<>();
            for (int i = 0; i < support.size(); i++) {
                multiply.add(support.get(i) * dnhyperResult.get(i));
            }
            double b = multiply.stream().mapToDouble(a -> a).sum();
            return b;
        }
    }

    private List<Double> dnhyper(Double ncp) {
        List<Double> result = new ArrayList<>();
        for (int i = 0; i < support.size(); i++) {
            result.add(logdc.get(i) + Math.log(ncp) * support.get(i));
        }
        double maxResult = Collections.max(result);
        List<Double> exponentResult = new ArrayList<>();

        for (double el : result) {
            exponentResult.add(Math.exp(el - maxResult));
        }
        result = new ArrayList<>();
        double sum = exponentResult.stream().mapToDouble(a -> a).sum();
        for (double element : exponentResult) {
            result.add(element / sum);
        }
        return result;
    }

    public String getOddRatio() {
        Double oddRatio = mle(x);
        if (oddRatio.isInfinite()) {
            return "Inf";
        } else if (oddRatio == Math.round(oddRatio)) {
            return new DecimalFormat("0").format(oddRatio);
        } else {
            return String.valueOf(round_as_r(oddRatio));
        }
    }

    public double getPValue() {
        return round_as_r(PvalueTwoSided);
    }

    public List<Double> getLogdc() {
        logdc = logdcDhyper(m, n, k);
        return logdc;
    }

    public double getPValueGreater() {
        return round_as_r(PvalueGreater);
    }

    public double getPValueLess() {
        return round_as_r(PvalueLess);
    }

    private double round_as_r(double value) {
        value = roundHalfEven("0", value * RESULT_ROUND_R);
        value = value/RESULT_ROUND_R;
        value = value == 0.0 ? 0 : (value == 1.0 ? 1 : value);
        return value;
    }

    private void calculatePValue() {
        PvalueLess = pnhyper(x, false);
        PvalueGreater = pnhyper(x, true);

        double relErr = 1 + 1E-7;
        List<Double> d = dnhyper(1.0);
        double sum = 0.0;
        for (Double el : d) {
            if (el <= d.get(x - lo) * relErr) {
                sum += el;
            }
        }
        PvalueTwoSided = sum;
    }

    private double pnhyper(int q, boolean upper_tail) {
        if (m + n == 0) {
            return 1.0;
        }
        if (upper_tail) {
            HypergeometricDistribution dhyper = new HypergeometricDistribution(m + n, m, k);
            return dhyper.upperCumulativeProbability(q);
        } else {
            HypergeometricDistribution dhyper = new HypergeometricDistribution(m + n, m, k);
            return dhyper.cumulativeProbability(q);
        }
    }
}


