package com.astrazeneca.vardict.collection;

import com.astrazeneca.vardict.variations.Variation;

import java.util.LinkedHashMap;
import java.util.Map;

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

    public static SV getSV(Map<Integer, VariationMap<String, Variation>> hash,
                    int start,
                    Variation vref) {
        VariationMap<String, Variation> map = hash.get(start);
        if (map == null) {
            map = new VariationMap<>();
            hash.put(start, map);
        }
        SV sv = map.sv;
        if (sv == null) {
            sv = new VariationMap.SV();
            map.sv = sv;
            map.put("SV", vref);
        }
        return sv;
    }
}
