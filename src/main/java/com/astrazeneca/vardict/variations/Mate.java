package com.astrazeneca.vardict.variations;

public class Mate {
    public int mateStart_ms;
    public int mateEnd_me;
    public int mateLength_mlen;
    public int start_s;
    public int end_e;
    public int pmean_rp;
    public double qmean_q;
    public int Qmean_Q;
    public int nm;

    public Mate() {
    }

    public Mate(int mateStart_ms, int mateEnd_me, int mateLength_mlen, int start_s, int end_e, int pmean_rp,
                double qmean_q, int qmean_Q, int nm) {
        this.mateStart_ms = mateStart_ms;
        this.mateEnd_me = mateEnd_me;
        this.mateLength_mlen = mateLength_mlen;
        this.start_s = start_s;
        this.end_e = end_e;
        this.pmean_rp = pmean_rp;
        this.qmean_q = qmean_q;
        this.Qmean_Q = qmean_Q;
        this.nm = nm;
    }

    @Override
    public String toString() {
        return "Mate{" +
                "mateStart_ms=" + mateStart_ms +
                ", mateEnd_me=" + mateEnd_me +
                ", mateLength_mlen=" + mateLength_mlen +
                ", start_s=" + start_s +
                ", end_e=" + end_e +
                ", pmean_rp=" + pmean_rp +
                ", qmean_q=" + qmean_q +
                ", Qmean_Q=" + Qmean_Q +
                ", nm=" + nm +
                '}';
    }
}
