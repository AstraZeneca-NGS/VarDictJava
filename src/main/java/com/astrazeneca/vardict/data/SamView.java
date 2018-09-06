package com.astrazeneca.vardict.data;

import htsjdk.samtools.*;

import java.util.HashMap;
import java.util.Map;

/**
 * To avoid performance issues only one SamView object on the same file can be opened per thread
 */
public class SamView implements AutoCloseable {
    private static ThreadLocal<Map<String, SamReader>> threadLocalSAMReaders = ThreadLocal.withInitial(HashMap::new);

    private SAMRecordIterator iterator;
    private int filter = 0;

    public SamView(String file, String samfilter, Region region, ValidationStringency stringency) {

        iterator = fetchReader(file, stringency)
                .queryOverlapping(region.chr, region.start, region.end);
        if (!"".equals(samfilter)) {
            filter = Integer.decode(samfilter);
        }
    }

    public SAMRecord read() {
        while(iterator.hasNext()) {
            SAMRecord record = iterator.next();
            if (filter != 0 && (record.getFlags() & filter) != 0) {
                continue;
            }
            return record;
        }
        return null;
    }

    @Override
    public void close() {
        iterator.close();
    }

    synchronized private static SamReader fetchReader(String file, ValidationStringency stringency) {
        return threadLocalSAMReaders.get().computeIfAbsent(
                file,
                (f) -> SamReaderFactory.makeDefault().validationStringency(stringency).open(SamInputResource.of(f))
        );
    }

}
