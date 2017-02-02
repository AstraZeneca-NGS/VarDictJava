package com.astrazeneca.vardict;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static com.astrazeneca.GlobalReadOnlyScope.instance;

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

    private static String[] retriveSubSeq(String fasta, String chr, int start, int end) throws IOException {
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
    //perl version: 568
    public static Map<Integer, Character> getReference(Region region) {
        Map<Integer, Character> ref = new HashMap<>();

        final int numberNucleotideToExtend = instance().conf.numberNucleotideToExtend;

        int s_start = region.start - numberNucleotideToExtend - 700 < 1 ?
                1 :
                region.start - numberNucleotideToExtend - 700;

        int len = instance().chrLengths.containsKey(region.chr) ?
                instance().chrLengths.get(region.chr) :
                0;

        int s_end = region.end + numberNucleotideToExtend + 700 > len ?
                len :
                region.end + numberNucleotideToExtend + 700;

        String[] subSeq = new String[0];
        try {
            subSeq = retriveSubSeq(instance().conf.fasta, region.chr, s_start, s_end);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
//        String header = subSeq[0];
        String exon = subSeq[1];
        for (int i = s_start; i < s_start + exon.length(); i++) { // TODO why '<=' ?
            ref.put(i, Character.toUpperCase(exon.charAt(i - s_start)));
        }

        return ref;
    }
}
