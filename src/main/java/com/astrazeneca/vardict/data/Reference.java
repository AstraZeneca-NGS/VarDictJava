package com.astrazeneca.vardict.data;

import java.util.*;

/**
 * Class storage reference sequences for the region and list of loaded region.
 * Used for comparing reference sequences and query sequences.
 */
public class Reference {
    public List<LoadedRegion> loadedRegions;
    public Map<Integer, Character> referenceSequences;
    public Map<String, List<Integer>> seed;

    public Reference() {
        this.loadedRegions = new ArrayList<>();
        this.referenceSequences = new HashMap<>();
        this.seed = new HashMap<>();
    }

    /**
     * Class to store already loaded region in reference to avoid excess operations on it
     */
    public static class LoadedRegion {
        public String chr;
        public int sequenceStart;
        public int sequenceEnd;

        public LoadedRegion(String chr, int sequenceStart, int sequenceEnd) {
            this.chr = chr;
            this.sequenceStart = sequenceStart;
            this.sequenceEnd = sequenceEnd;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            LoadedRegion that = (LoadedRegion) o;
            return sequenceStart == that.sequenceStart &&
                    sequenceEnd == that.sequenceEnd &&
                    Objects.equals(chr, that.chr);
        }

        @Override
        public int hashCode() {
            return Objects.hash(chr, sequenceStart, sequenceEnd);
        }
    }
}
