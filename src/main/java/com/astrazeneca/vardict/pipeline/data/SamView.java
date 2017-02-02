package com.astrazeneca.vardict.pipeline.data;

import com.astrazeneca.vardict.Region;
import htsjdk.samtools.*;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class SamView implements AutoCloseable {

    private SAMRecordIterator iterator;
    private int filter = 0;

    private static ThreadLocal<Map<String, SamReader>> threadLocalSAMReaders = ThreadLocal.withInitial(HashMap::new);

    synchronized private static SamReader fetchReader(String file, ValidationStringency stringency) {
        return threadLocalSAMReaders.get().computeIfAbsent(
                file,
                (f) -> SamReaderFactory.makeDefault().validationStringency(stringency).open(SamInputResource.of(f))
        );
    }

    public SamView(String file, String samfilter, Region region, ValidationStringency stringency) {
        iterator = fetchReader(file, stringency)
                .queryOverlapping(region.chr, region.start, region.end);
        if (!"".equals(samfilter)) {
            filter = Integer.decode(samfilter);
        }
    }

    public SAMRecord read() throws IOException {
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
    public void close() throws IOException {
        iterator.close();
    }
}
