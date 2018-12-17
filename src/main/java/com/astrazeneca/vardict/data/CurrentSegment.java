package com.astrazeneca.vardict.data;

/**
 * Contains information about current processed segment (chromosome, start and end of the region)
 */
public class CurrentSegment {
    public String chr;
    public int start;
    public int end;

    public CurrentSegment(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }
}
