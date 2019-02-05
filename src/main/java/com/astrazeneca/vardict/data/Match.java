package com.astrazeneca.vardict.data;

/**
 * Matched information between sequence and reference/ Contains base position and matched sequence
 */
public class Match {
    public int basePosition; //$bp
    public String matchedSequence;

    public Match(int basePosition, String matchedSequence) {
        this.basePosition = basePosition;
        this.matchedSequence = matchedSequence;
    }
}
