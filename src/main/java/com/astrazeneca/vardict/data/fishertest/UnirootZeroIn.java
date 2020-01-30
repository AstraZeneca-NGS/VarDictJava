package com.astrazeneca.vardict.data.fishertest;


import java.util.function.Function;

/**
 * EXPERIMENTAL FEATURE.
 * <p>
 * Implementation of zeroin function from Fortran Netlib and C libraries
 */

public class UnirootZeroIn {
    /**
     * An estimate to the root from zeroin C implementation
     *
     * @param ax  Left border of the range
     * @param bx  Right border the root is seeked
     * @param f   Function under investigation
     * @param tol Acceptable tolerance
     * @return root
     */
    public static double zeroinC(double ax, double bx, Function<Double, Double> f, double tol) {
        double a, b, c;           /* Abscissae, descr. see above	*/
        double fa;                /* f(a)				*/
        double fb;                /* f(b)				*/
        double fc;                /* f(c)				*/
        double EPSILON = Math.ulp(1.0);
        a = ax;
        b = bx;
        fa = f.apply(a);
        fb = f.apply(b);
        c = a;
        fc = fa;

        /* Main iteration loop	*/
        for (; ; ) {
            double prev_step = b - a;   /* Distance from the last but one to the last approximation	*/
            double tol_act;             /* Actual tolerance		*/
            double p;                   /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment */
            double q;
            double new_step;            /* Step at this iteration */

            /* Swap data for b to be the best approximation	*/
            if (Math.abs(fc) < Math.abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol_act = 2 * EPSILON * Math.abs(b) + tol / 2.0;
            new_step = (c - b) / 2.0;

            /* Acceptable approx. is found	*/
            if (Math.abs(new_step) <= tol_act || fb == 0.0) return b;

            /* Decide if the interpolation can be tried. If prev_step was large enough and was in true direction.
            Interpolatiom may be tried */
            if (Math.abs(prev_step) >= tol_act && Math.abs(fa) > Math.abs(fb)) {
                double t1, cb, t2;
                cb = c - b;
                if (a == c) {        /* If we have only two distinct points linear interpolation can only be applied */
                    t1 = fb / fa;
                    p = cb * t1;
                    q = 1.0 - t1;
                } else {             /* Quadric inverse interpolation*/
                    q = fa / fc;
                    t1 = fb / fc;
                    t2 = fb / fa;
                    p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                    q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
                }
                if (p > 0.0) q = -q;        /* p was calculated with the opposite sign; make p positive*/
                else p = -p;                /* and assign possible minus to	q*/

                /* If b+p/q falls in [b,c] and isn't too large it is accepted	*/
                /* If p/q is too large then the	bissection procedure can reduce [b,c] range to more extent*/
                if (p < (0.75 * cb * q - Math.abs(tol_act * q) / 2.0) && p < Math.abs(prev_step * q / 2.0)) {
                    new_step = p / q;
                }
            }

            /* Adjust the step to be not less than tolerance*/
            if (Math.abs(new_step) < tol_act) {
                if (new_step > 0.0) new_step = tol_act;
                else new_step = -tol_act;
            }

            a = b;
            fa = fb;            /* Save the previous approx.	*/
            b += new_step;
            fb = f.apply(b);    /* Do step to a new approxim.	*/

            /* Adjust c for it to have a sign*/
            if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
                c = a;
                fc = fa;
            }
        }
    }
}


