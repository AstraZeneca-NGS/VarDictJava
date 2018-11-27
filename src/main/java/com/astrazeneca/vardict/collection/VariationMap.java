package com.astrazeneca.vardict.collection;

import com.astrazeneca.vardict.variations.Variation;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * LinkedHashMap contains additional field to store information about found Structural variants
 * @param <K> commonly String - description string of variant
 * @param <V> commonly Variation
 */
public class VariationMap<K,V> extends LinkedHashMap<K,V> {
    public SV sv;

    /**
     * Structural variant content
     */
    public static class SV {
        public String type;
        public int pairs;
        public int splits;
        public int clusters;
    }

    /**
     * If map already contains SV, return is. Else it creates SV and put new Variation to description String "SV"
     * @param hash map of description strings and variations on position
     * @param start start position to search for structural variants
     * @return structural variant data (pairs, clusters, splits).
     */
    public static SV getSV(Map<Integer, VariationMap<String, Variation>> hash,
                           int start) {
        VariationMap<String, Variation> map = hash.get(start);
        if (map == null) {
            map = new VariationMap<>();
            hash.put(start, map);
        }
        SV sv = map.sv;
        if (sv == null) {
            sv = new VariationMap.SV();
            map.sv = sv;
            map.put("SV", new Variation());
        }
        if (!map.containsKey("SV")) {
            map.put("SV", new Variation());
        }
        return sv;
    }

    /**
     * Removes SV description string from map and delete information about SV.
     * @param hash map of description strings and variations on position
     * @param start start position to search for structural variants
     */
    public static void removeSV (Map<Integer, VariationMap<String, Variation>> hash, int start) {
        hash.get(start).sv = null;
        if (hash.get(start).containsKey("SV")) {
            hash.get(start).remove("SV");
        }
    }
}
