package com.astrazeneca.vardict;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.LocalDateTime;
import java.util.*;

public class ReferenceResource {
    private static final int DEFAULT_REFERENCE_EXTENSION = 1200;
    private static ThreadLocal<Map<String, IndexedFastaSequenceFile>> threadLocalFastaFiles = ThreadLocal.withInitial(HashMap::new);

    synchronized private static IndexedFastaSequenceFile fetchFasta(String file) {
        return threadLocalFastaFiles.get().computeIfAbsent(
                file,
                (f) -> {
                    try {
                        return new IndexedFastaSequenceFile(new File((f)));
                    } catch (FileNotFoundException e) {
                        throw new IllegalArgumentException("Couldn't open reference file: " + f, e);
                    }
                }
        );
    }

    /**
     * Method convert fasta data to array contains header and nucleotide bases.
     * @param fasta path to fasta file
     * @param chr chromosome name of region
     * @param start start position of region
     * @param end end position of region
     * @return array of nucleotide bases in the region of fasta
     */
    public static String[] retrieveSubSeq(String fasta, String chr, int start, int end) {
        IndexedFastaSequenceFile idx = fetchFasta(fasta);
        ReferenceSequence seq = idx.getSubsequenceAt(chr, start, end);
        byte[] bases = seq.getBases();
        return new String[] { ">" + chr + ":" + start + "-" + end, bases != null ? new String(bases) : "" };
    }

    /**
     * Get part of reference sequence in
     * @param region region of interest
     * @param chrs map of chromosome lengths
     * @param conf Vardict Configuration
     * @return reference object contains sequence map: key - position, value - base and seed map
     */
    public static Reference getREF(Region region, Map<String, Integer> chrs, Configuration conf) {
        Reference ref = new Reference();
        return getREF(region, chrs, conf, DEFAULT_REFERENCE_EXTENSION, ref);
    }

    /**
     * Get part of reference sequence
     * @param region region of interest
     * @param chrs map of chromosome lengths
     * @param conf Vardict Configuration contains number of nucleotides to extend and seed1, seed2 lengths
     * @param extension extension for reference
     * @param ref reference
     * @return reference object contains sequence map: key - position, value - base and seed map
     */
    public static Reference getREF(Region region, Map<String, Integer> chrs, Configuration conf, int extension, Reference ref) {
        int sequenceStart = region.start - conf.numberNucleotideToExtend - extension < 1 ? 1
                : region.start - conf.numberNucleotideToExtend - extension;
        int len = chrs.containsKey(region.chr) ? chrs.get(region.chr) : 0;
        int sequenceEnd = region.end + conf.numberNucleotideToExtend + extension > len ?
                len : region.end + conf.numberNucleotideToExtend + extension;
        if (conf.y) {
            System.err.println("TIME: Getting REF: " + LocalDateTime.now());
        }

        String[] subSeq = retrieveSubSeq(conf.fasta, region.chr, sequenceStart, sequenceEnd);
        //Header doesn't used
        //String header = subSeq[0];
        String exon = subSeq[1];

        if (isLoaded(region.chr, sequenceStart, sequenceEnd, ref)) {
            return ref;
        }
        Reference.LoadedRegion loadedRegion = new Reference.LoadedRegion(region.chr, sequenceStart, sequenceEnd);
        ref.loadedRegions.add(loadedRegion);

        for (int i = 0; i < exon.length(); i++) { // TODO why '<=' in Perl?
            // don't process it more than once
            if (ref.referenceSequences.containsKey(i + sequenceStart)) {
                continue;
            }
            ref.referenceSequences.put(i + sequenceStart, Character.toUpperCase(exon.charAt(i)));

            // Fill the seed map by sequences of seed1 and seed2 length
            if (exon.length() - i > conf.seed1) {
                String keySequence = exon.substring(i, i + conf.seed1).toUpperCase();
                ref = addPositionsToSeedSequence(ref, sequenceStart, i, keySequence);
            }
            if (exon.length() - i > conf.seed2) {
                String keySequence = exon.substring(i, i + conf.seed2).toUpperCase();
                ref = addPositionsToSeedSequence(ref, sequenceStart, i, keySequence);
            }
        }

        if (conf.y) {
            System.err.println("TIME: Got REF: " + LocalDateTime.now());
        }

        return ref;
    }

    /**
     * Method adds the key sequence for current position in the seed map
     * @param ref Reference object
     * @param sequenceStart position of sequence in chromosome
     * @param i current position in exon to add in seed Map
     * @param keySequence sequence length of seed1 parameter from chromosome
     * @return updated Reference
     */
    public static Reference addPositionsToSeedSequence(Reference ref, int sequenceStart, int i, String keySequence) {
        List<Integer> seedPositions = ref.seed.getOrDefault(keySequence, new ArrayList<>());
        seedPositions.add(i + sequenceStart);
        ref.seed.put(keySequence, seedPositions);
        return ref;
    }

    /**
     * Check whether a region is already loaded in reference
     * @return true if region is already loaded and false if not
     */
    public static boolean isLoaded(String chr, int sequenceStart, int sequenceEnd, Reference reference) {
        if (reference.loadedRegions.isEmpty()) {
            return false;
        }
        for (Reference.LoadedRegion region : reference.loadedRegions) {
            if (chr.equals(region.chr) && sequenceStart >= region.sequenceStart && sequenceEnd <= region.sequenceEnd) {
                return true;
            }
        }
        return false;
    }
}
