package com.astrazeneca.vardict.data;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.exception.RegionBoundariesException;
import com.astrazeneca.vardict.exception.WrongFastaOrBamException;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.LocalDateTime;
import java.util.*;

import static com.astrazeneca.vardict.data.Patterns.UNABLE_FIND_CONTIG;
import static com.astrazeneca.vardict.data.Patterns.WRONG_START_OR_END;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.substr;

/**
 * Utility class for access to reference sequences
 */
public class ReferenceResource {
    /**
     * Fasta files store in thread local variables to avoid multithreading issues.
     */
    private ThreadLocal<Map<String, IndexedFastaSequenceFile>> threadLocalFastaFiles = ThreadLocal.withInitial(HashMap::new);

    synchronized private IndexedFastaSequenceFile fetchFasta(String file) {
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
    public String[] retrieveSubSeq(String fasta, String chr, int start, int end) {
        try {
            IndexedFastaSequenceFile idx = fetchFasta(fasta);
            ReferenceSequence seq = idx.getSubsequenceAt(chr, start, end);
            byte[] bases = seq.getBases();
            return new String[]{">" + chr + ":" + start + "-" + end, bases != null ? new String(bases) : ""};
        } catch (SAMException e){
            if (UNABLE_FIND_CONTIG.matcher(e.getMessage()).find()){
                throw new WrongFastaOrBamException(chr, e);
            } else if (WRONG_START_OR_END.matcher(e.getMessage()).find()){
                throw new RegionBoundariesException(chr, start, end, e);
            } else {
                throw e;
            }
        }
    }

    /**
     * Get part of reference sequence with default extension.
     * @param region region of interest
     * @return reference object contains sequence map: key - position, value - base and seed map
     */
    public Reference getReference(Region region) {
        Reference ref = new Reference();
        return getReference(region, instance().conf.referenceExtension, ref);
    }

    /**
     * Get part of reference sequence with customizable extension.
     * @param region region of interest
     * @param extension extension for reference
     * @param ref reference
     * @return reference object contains sequence map: key - position, value - base and seed map
     */
    public Reference getReference(Region region, int extension, Reference ref) {
        int sequenceStart = region.start - instance().conf.numberNucleotideToExtend - extension < 1 ? 1
                : region.start - instance().conf.numberNucleotideToExtend - extension;
        int len = instance().chrLengths.containsKey(region.chr) ? instance().chrLengths.get(region.chr) : 0;
        int sequenceEnd = region.end + instance().conf.numberNucleotideToExtend + extension > len ?
                len : region.end + instance().conf.numberNucleotideToExtend + extension;
        if (instance().conf.y) {
            System.err.println("TIME: Getting REF: " + LocalDateTime.now());
        }

        String[] subSeq = retrieveSubSeq(instance().conf.fasta, region.chr, sequenceStart, sequenceEnd);
        //Header doesn't used
        //String header = subSeq[0];
        String exon = subSeq[1].toUpperCase();

        if (isLoaded(region.chr, sequenceStart, sequenceEnd, ref)) {
            return ref;
        }
        Reference.LoadedRegion loadedRegion = new Reference.LoadedRegion(region.chr, sequenceStart, sequenceEnd);
        ref.loadedRegions.add(loadedRegion);

        // To process ends of chromosomes without decreasing by SEED1
        int siteEnd = len == sequenceEnd ? exon.length() : exon.length() - Configuration.SEED_1;
        for (int i = 0; i < siteEnd; i++) {
            // don't process it more than once
            if (ref.referenceSequences.containsKey(i + sequenceStart)) {
                continue;
            }
            ref.referenceSequences.put(i + sequenceStart, exon.charAt(i));

            // Do not create adaptor sequences for the very end of chromosome
            if (len == sequenceEnd && i > exon.length() - Configuration.SEED_1) {
                continue;
            }
            // Fill the seed map by sequences of SEED_1 and SEED_2 length
            String keySequence = substr(exon, i, Configuration.SEED_1);
            ref = addPositionsToSeedSequence(ref, sequenceStart, i, keySequence);

            keySequence = substr(exon, i, Configuration.SEED_2);
            ref = addPositionsToSeedSequence(ref, sequenceStart, i, keySequence);

        }

        if (instance().conf.y) {
            System.err.println("TIME: Got REF: " + LocalDateTime.now());
        }

        return ref;
    }

    /**
     * Method adds the key sequence for current position in the seed map
     * @param ref Reference object
     * @param sequenceStart position of sequence in chromosome
     * @param i current position in exon to add in seed Map
     * @param keySequence sequence length of SEED_1 parameter from chromosome
     * @return updated Reference
     */
    private Reference addPositionsToSeedSequence(Reference ref, int sequenceStart, int i, String keySequence) {
        List<Integer> seedPositions = ref.seed.getOrDefault(keySequence, new ArrayList<>());
        seedPositions.add(i + sequenceStart);
        ref.seed.put(keySequence, seedPositions);
        return ref;
    }

    /**
     * Check whether a region is already loaded in reference
     * @param chr string name of the chromosome
     * @param sequenceStart start of the region of interest
     * @param sequenceEnd end of the region of interest
     * @param reference Reference object
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
