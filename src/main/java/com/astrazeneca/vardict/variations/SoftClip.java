package com.astrazeneca.vardict.variations;

import java.util.Map;
import java.util.TreeMap;

public class SoftClip extends Variation {
    public TreeMap<Integer, Map<Character, Integer>> nt = new TreeMap<>();
    public TreeMap<Integer, Map<Character, Variation>> seq = new TreeMap<>();
    public String sequence;
    public boolean used;
}
