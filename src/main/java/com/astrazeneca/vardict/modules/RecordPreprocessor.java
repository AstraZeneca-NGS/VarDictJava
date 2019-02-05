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
    public int dupReads;

    private SamView currentReader;
    private HashSet<String> dup;
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

    public SAMRecord next() {

        SAMRecord record = nextOnCurrentReader();
        if (record == null) {
            nextReader();
            record = nextOnCurrentReader();
        }

        return record;
    }

    public static String getChrName(Region region) {
        String chrName;
        if (instance().conf.chromosomeNameIsNumber && region.chr.startsWith("chr")) { //remove prefix 'chr' if option -C is set
            chrName = region.chr.substring("chr".length());
        } else {
            chrName = region.chr;
        }
        return chrName;
    }

    private void nextReader() {
        if (bams.isEmpty())
            return;

        if (currentReader != null) {
            currentReader.close();
        }
        String samfilter = instance().conf.samfilter == null || instance().conf.samfilter.isEmpty() ? "" : instance().conf.samfilter;
        currentReader = new SamView(bams.pollLast(), samfilter, region, instance().conf.validationStringency);

        //dup contains already seen reads. For each seen read dup contains either POS-RNEXT-PNEXT or POS-CIGAR (if next segment in template is unmapped).
        dup = new HashSet<>();
        //position of first matching base (POS in SAM)
        firstMatchingPosition = -1;
    }

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

    private boolean preprocessRecord(SAMRecord record) {
        if (instance().conf.isDownsampling() && RND.nextDouble() <= instance().conf.downsampling) {
            return false;
        }

        final String querySequence = record.getReadString();
        // TODO: htsjdk already has this functionality
        final Flags flag = new Flags(record.getFlags());
        final int mappingQuality = record.getMappingQuality();

        if (instance().conf.hasMappingQuality() && mappingQuality < instance().conf.mappingQuality) { // ignore low mapping quality reads
            return false;
        }

        if (flag.isNotPrimaryAlignment() && instance().conf.samfilter != null) {
            return false;
        }

        if (querySequence.length() == 1 && querySequence.charAt(0) == '*') {
            return false;
        }
        totalReads++;

        final String mateReferenceName = getMateReferenceName(record);

        // filter duplicated reads if option -t is set
        if (instance().conf.removeDuplicatedReads) {
            if (record.getAlignmentStart() != firstMatchingPosition) {
                dup.clear();
            }
            if (record.getMateAlignmentStart() < 10) {
                String dupKey = record.getAlignmentStart() + "-" + mateReferenceName + "-" + record.getMateAlignmentStart();
                if (dup.contains(dupKey)) {
                    dupReads++;
                    return false;
                }
                dup.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            } else if (flag.isUnmappedMate()) {
                String dupKey = record.getAlignmentStart() + "-" + record.getCigarString();
                if (dup.contains(dupKey)) {
                    dupReads++;
                    return false;
                }
                dup.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            }
        }

        return true;
    }
}
