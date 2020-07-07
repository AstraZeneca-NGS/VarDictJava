package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.SAMRecord;

import java.util.*;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.modules.CigarParser.getMateReferenceName;

/**
 * Class for primary processing of SAM record and filtering reads that don't fit the requirements (low quality,
 * not primary alignment, etc.).
 */
public class RecordPreprocessor {
    private final static Random RND = new Random(System.currentTimeMillis());

    final Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    final Map<Integer, VariationMap<String, Variation>> insertionVariants;
    final Map<Integer, Integer> refCoverage;
    final Map<Integer, Sclip> softClips5End;
    final Map<Integer, Sclip> softClips3End;

    private final Deque<String> bams;
    private final Region region;
    public int totalReads;
    public int duplicateReads;

    private SamView currentReader;
    private HashSet<String> duplicates; //$dup
    private int firstMatchingPosition;

    public RecordPreprocessor(String[] bams, Region region, InitialData data) {
        this.nonInsertionVariants = data.nonInsertionVariants;
        this.insertionVariants = data.insertionVariants;
        this.refCoverage = data.refCoverage;
        this.softClips5End = data.softClips5End;
        this.softClips3End = data.softClips3End;

        this.bams = new ArrayDeque<>();
        for (String bam : bams) {
            this.bams.push(bam);
        }
        this.region = region;

        nextReader();
    }

    /**
     * Updates current reader of SAM/BAM file (SamView) and gets the record from it.
     * @return current SAM Record to process in CigarParser
     */
    public SAMRecord nextRecord() {
        SAMRecord record = nextOnCurrentReader();
        if (record == null) {
            nextReader();
            record = nextOnCurrentReader();
        }
        return record;
    }

    /**
     * Creates SamView which will contains iterator for reads of this region in SAM/BAM file.
     */
    private void nextReader() {
        if (bams.isEmpty())
            return;

        if (currentReader != null) {
            currentReader.close();
        }
        currentReader = new SamView(bams.pollLast(), instance().conf.samfilter, region, instance().conf.validationStringency);

        //$dup contains already seen reads. For each seen read dup contains either POS-RNEXT-PNEXT or POS-CIGAR
        // (if next segment in template is unmapped).
        duplicates = new HashSet<>();
        //position of first matching base (POS in SAM)
        firstMatchingPosition = -1;
    }

    /**
     * While there are records in SAM/BAM file, processed them.
     * @return next SAM record from the SAM/BAM file
     */
    private SAMRecord nextOnCurrentReader() {
        SAMRecord record;
        do {
            if ((record = currentReader.read()) == null) {
                return null;
            }
        } while (!preprocessRecord(record));

        return record;
    }

    public void close() {
        currentReader.close();
    }

    /**
     * Filtered each record if it doesn't fit the requirements and calculates total reads and duplicate reads numbers.
     * @param record SAM Record which is processed
     * @return true if record must be processed further and false if it must be skipped.
     */
    private boolean preprocessRecord(SAMRecord record) {
        if (instance().conf.isDownsampling() && RND.nextDouble() <= instance().conf.downsampling) {
            return false;
        }
        final String querySequence = record.getReadString();
        final int mappingQuality = record.getMappingQuality();

        // Ignore low mapping quality reads
        if (instance().conf.hasMappingQuality() && mappingQuality < instance().conf.mappingQuality) {
            return false;
        }

        //Skip not primary alignment reads
        if (record.isSecondaryAlignment() && !instance().conf.samfilter.equals("0")) {
            return false;
        }
        // Skip reads where sequence is not stored in read
        if (querySequence.length() == 1 && querySequence.charAt(0) == '*') {
            return false;
        }
        totalReads++;

        final String mateReferenceName = getMateReferenceName(record);

        // filter duplicated reads if option -t is set
        if (instance().conf.removeDuplicatedReads) {
            if (record.getAlignmentStart() != firstMatchingPosition) {
                duplicates.clear();
            }
            if (record.getMateAlignmentStart() < 10) {
                //POS-RNEXT-PNEXT
                String dupKey = record.getAlignmentStart() + "-" + mateReferenceName + "-" + record.getMateAlignmentStart();
                if (duplicates.contains(dupKey)) {
                    duplicateReads++;
                    return false;
                }
                duplicates.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            } else if (record.getReadPairedFlag() && record.getMateUnmappedFlag()) {
                //POS-CIGAR
                String dupKey = record.getAlignmentStart() + "-" + record.getCigarString();
                if (duplicates.contains(dupKey)) {
                    duplicateReads++;
                    return false;
                }
                duplicates.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            }
        }
        return true;
    }

    /**
     * Fixes chromosome name. If option -C set, it will remove prefix "chr".
     * @param region region with chromosome
     * @return fixed chromosome name
     */
    public static String getChrName(Region region) {
        String chrName;
        if (instance().conf.chromosomeNameIsNumber && region.chr.startsWith("chr")) { //remove prefix 'chr' if option -C is set
            chrName = region.chr.substring("chr".length());
        } else {
            chrName = region.chr;
        }
        return chrName;
    }

}
