package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.vardict.pipeline.data.Flags;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.SamView;
import htsjdk.samtools.SAMRecord;

import java.io.IOException;
import java.util.*;

import static com.astrazeneca.GlobalReadOnlyScope.instance;

public class RecordPreprocessor {

    private final static Random RND = new Random(System.currentTimeMillis());

    private final Deque<String> bams;
    private final Region region;

    private SamView currentReader;
    private HashSet<String> dup;
    private int firstMatchingPosition;

    public RecordPreprocessor(String[] bams, Region region) {
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
            try {
                currentReader.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        String samfilter = instance().conf.samfilter == null || instance().conf.samfilter.isEmpty() ? "" : instance().conf.samfilter;
        currentReader = new SamView(bams.pollLast(), samfilter, region, instance().conf.validationStringency);

        dup = new HashSet<>();
        firstMatchingPosition = -1;
    }

    private SAMRecord nextOnCurrentReader() {
        SAMRecord record;
        do {
            try {
                if ((record = currentReader.read()) == null) {
                    return null;
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        } while (!preprocessRecord(record));

        return record;
    }

    public void close() {
        try {
            currentReader.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private boolean preprocessRecord(SAMRecord record) {

        // perl version: 609
        //dup contains already seen reads. For each seen read dup contains either POS-RNEXT-PNEXT or POS-CIGAR (if next segment in template is unmapped).
        //position of first matching base (POS in SAM)
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

        final String mateReferenceName = SAMFileParser.getMateReferenceName(record);

        // filter duplicated reads if option -t is set
        if (instance().conf.removeDuplicatedReads) {
            if (record.getAlignmentStart() != firstMatchingPosition) {
                dup.clear();
            }
            if (record.getMateAlignmentStart() < 10) {
                String dupKey = record.getAlignmentStart() + "-" + mateReferenceName + "-" + record.getMateAlignmentStart();
                if (dup.contains(dupKey)) {
                    return false;
                }
                dup.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            } else if (flag.isUnmappedMate()) {
                String dupKey = record.getAlignmentStart() + "-" + record.getCigarString();
                if (dup.contains(dupKey)) {
                    return false;
                }
                dup.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            }
        }

        return true;
    }

}
