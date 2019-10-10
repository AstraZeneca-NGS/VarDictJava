package com.astrazeneca.vardict.data;

public class ModifiedCigar {
    public int position;
    public String cigar;
    public String querySequence;
    public String queryQuality;

    public ModifiedCigar(int position, String cigar,String querySequence, String queryQuality) {
        this.position = position;
        this.cigar = cigar;
        this.querySequence = querySequence;
        this.queryQuality = queryQuality;
    }
}
