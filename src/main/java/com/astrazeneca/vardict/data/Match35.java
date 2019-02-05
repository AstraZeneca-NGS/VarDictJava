package com.astrazeneca.vardict.data;

/**
 * Contains information about matching position on 3 and 5 ends of sequence
 */
public class Match35 {
    /**
     * Start matched position of 5'
     */
    public int matched5end;
    /**
     * Start matched position vof 3'
     */
    public int matched3End;
    /**
     * maximum matched length
     */
    public int maxMatchedLength;

    public Match35(int matched5end, int matched3End, int maxMatchedLength) {
        this.matched5end = matched5end;
        this.matched3End = matched3End;
        this.maxMatchedLength = maxMatchedLength;
    }
}