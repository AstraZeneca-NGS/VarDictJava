package com.astrazeneca.vardict.variations;

import java.util.*;

public class Sclip extends Variation {
    public TreeMap<Integer, Map<Character, Integer>> nt = new TreeMap<>();
    public TreeMap<Integer, Map<Character, Variation>> seq = new TreeMap<>();
    public String sequence;
    public boolean used;

    //structural variation fields
    public int start;
    public int end;
    public int mstart;
    public int mend;
    public int mlen;
    public int disc;
    public int softp;
    public Map<Integer, Integer> soft = new HashMap<>();
    public List<Mate> mates = new ArrayList<>();
}

