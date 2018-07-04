package com.astrazeneca.vardict;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class ReferenceResource {

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
     *
     * @param fasta path to fasta file
     * @param chr chromosome name of region
     * @param start start position of region
     * @param end end position of region
     * @return array of nucleotide bases in the region of fasta
     */
    public static String[] retriveSubSeq(String fasta, String chr, int start, int end) {
        IndexedFastaSequenceFile idx = fetchFasta(fasta);
        ReferenceSequence seq = idx.getSubsequenceAt(chr, start, end);
        byte[] bases = seq.getBases();
        return new String[] { ">" + chr + ":" + start + "-" + end, bases != null ? new String(bases) : "" };
    }

    /**
     * Get part of reference sequence
     * @param region region of interest
     * @return reference sequence map: key - position, value - base
     * @throws IOException
     */
    public static Map<Integer, Character> getREF(Region region, Map<String, Integer> chrs, String fasta, int numberNucleotideToExtend) {
        Map<Integer, Character> ref = new HashMap<>();

        int s_start = region.start - numberNucleotideToExtend - 700 < 1 ? 1 : region.start - numberNucleotideToExtend - 700;
        int len = chrs.containsKey(region.chr) ? chrs.get(region.chr) : 0;
        int s_end = region.end + numberNucleotideToExtend + 700 > len ?
                len : region.end + numberNucleotideToExtend + 700;

        String[] subSeq = retriveSubSeq(fasta, region.chr, s_start, s_end);
//        String header = subSeq[0];
        String exon = subSeq[1];
        for (int i = s_start; i < s_start + exon.length(); i++) { // TODO why '<=' ?
            ref.put(i, Character.toUpperCase(exon.charAt(i - s_start)));
        }

        return ref;
    }
}
