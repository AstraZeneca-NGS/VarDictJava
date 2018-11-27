package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;

/**
 * Class to start the process of SAM/BAM file with initializing of record processing for current region.
 */
public class SAMFileParser implements Module<InitialData, RecordPreprocessor> {
    @Override
    public Scope<RecordPreprocessor> process(Scope<InitialData> scope) {
        return new Scope<>(
                scope,
                new RecordPreprocessor(scope.bam.split(":"), scope.region, scope.data)
        );
    }
}
