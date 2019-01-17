package com.astrazeneca.vardict.variations;

import java.util.*;

/**
 * Class to store data about softclips.
 */
public class Sclip extends Variation {
    /**
     * Map of position of high quality base in the base sequence (from SAM record)
     * to base on this position and it's count
     */
    public TreeMap<Integer, TreeMap<Character, Integer>> nt = new TreeMap<>();
    /**
     * Map of position of high quality base in the base sequence (from SAM record)
     * to base on this position and it's variation
     */
    public TreeMap<Integer, Map<Character, Variation>> seq = new TreeMap<>();

    /**
     * The consensus sequence in soft-clipped reads.
     */
    public String sequence;
    public boolean used;

    /**
     * Additional fields for structural variation
     * //TODO: can be extracted to separate class but collections must be changed
     */
    public int start;
    public int end;
    public int mstart;
    public int mend;
    public int mlen;
    public int disc;
    public int softp;

    /**
     * Map of softclip positions to their counts for SV
     */
    public Map<Integer, Integer> soft = new LinkedHashMap<>();
    public List<Mate> mates = new ArrayList<>();
}

