package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.variations.SoftClip;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.VariationUtils;
import htsjdk.samtools.SAMRecord;

import java.util.Map;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Utils.charAt;
import static com.astrazeneca.utils.Utils.substr;
import static java.lang.String.format;

/**
 * Construct a variant structure.
 */
public class SAMFileParser implements Module<Object, RecordPreprocessor> {

    @Override
    public Scope<RecordPreprocessor> process(Scope<Object> scope) {
        return new Scope<>(
                scope,
                new RecordPreprocessor(scope.bam.split(":"), scope.region)
        );
    }

    /**
     * Adjust MNP when there're breakpoints within MNP (multi-nucleotide polymorphism)
     * @param nonInsertionVariants noninsertion variant structure
     * @param mnp multi-nucleotide polymorphism structure
     * @param refCoverage reference coverage
     * @param reference reference bases by position
     * @param softClip5End soft clipped at 5'
     * @param softClip3End soft clipped at 3'
     */
    //perl version: 1659
    static void adjustMNP(Map<Integer, Map<String, Variation>> nonInsertionVariants,
                          Map<Integer, Map<String, Integer>> mnp,
                          Map<Integer, Integer> refCoverage,
                          Map<Integer, Character> reference,
                          Map<Integer, SoftClip> softClip3End,
                          Map<Integer, SoftClip> softClip5End) {
        Configuration configuration = instance().conf;
        for (Map.Entry<Integer, Map<String, Integer>> entry : mnp.entrySet()) {
            final Integer p = entry.getKey();
            Map<String, Integer> v = entry.getValue();

            for (Map.Entry<String, Integer> en : v.entrySet()) {
                final String vn = en.getKey();
                final Map<String, Variation> hashP = nonInsertionVariants.get(p);
                if (hashP == null) {
                    continue;
                }
                final Variation vref = hashP.get(vn);
                if (vref == null ) { // The variant is likely already been used by indel realignment
                    continue;
                }
                final String mnt = vn.replaceFirst("&", "");
                for (int i = 0; i < mnt.length() - 1; i++) {
                    String left = substr(mnt, 0, i + 1);
                    String right = substr(mnt, -(mnt.length() - i - 1));
                    {
                        Variation tref = hashP.get(left);
                        if (tref != null) {
                            if (tref.varsCount < vref.varsCount && tref.meanPosition / tref.varsCount <= i + 1) {
                                if (configuration.y) {
                                    System.err.printf(" AdjMnt Left: %s %s %s\n", p, vn, tref.varsCount);
                                }
                                VariationUtils.adjCnt(vref, tref);
                                hashP.remove(left);
                            }
                        }
                    }
                    if (nonInsertionVariants.containsKey(p + i + 1)) {
                        Variation tref = nonInsertionVariants.get(p + i + 1).get(right);
                        if (tref != null) {
                            if (tref.varsCount < vref.varsCount && tref.meanPosition / tref.varsCount <= mnt.length() - i - 1) {
                                if (configuration.y) {
                                    System.err.printf(" AdjMnt Right: %s %s %s\n", p, vn, tref.varsCount);
                                }
                                VariationUtils.adjCnt(vref, tref);
                                VariationUtils.incCnt(refCoverage, p, tref.varsCount);
                                nonInsertionVariants.get(p + i + 1).remove(right);
                            }
                        }
                    }

                }
                if (softClip3End.containsKey(p)) {
                    final SoftClip sc3v = softClip3End.get(p);
                    if (!sc3v.used) {
                        final String seq = VariationUtils.findconseq(sc3v);
                        if (seq.startsWith(mnt)) {
                            if(seq.length() == mnt.length() || isMatchRef(seq.substring(mnt.length()), reference, p + mnt.length(), 1)) {
                                VariationUtils.adjCnt(nonInsertionVariants.get(p).get(vn), sc3v);
                                VariationUtils.incCnt(refCoverage, p, sc3v.varsCount);
                                sc3v.used = true;
                            }
                        }
                    }
                }
                if (softClip5End.containsKey(p + mnt.length())) {
                    final SoftClip sc5v = softClip5End.get(p + mnt.length());
                    if (!sc5v.used) {
                        String seq =  VariationUtils.findconseq(sc5v);
                        if (!seq.isEmpty() && seq.length() >= mnt.length()) {
                            seq =  new StringBuffer(seq).reverse().toString();
                            if (seq.endsWith(mnt)) {
                                if (seq.length() == mnt.length() || isMatchRef(seq.substring(0, seq.length() - mnt.length()), reference, p - 1, -1)) {
                                    VariationUtils.adjCnt(nonInsertionVariants.get(p).get(vn), sc5v);
                                    VariationUtils.incCnt(refCoverage, p, sc5v.varsCount);
                                    sc5v.used = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //perl version: 2765
    //utility method for adjustMNP method
    static boolean isMatchRef(String sequence, // subsequence consensus sequence in soft-clipped reads
                                      Map<Integer, Character> reference, // reference bases by position
                                      int position, // key for MNP map
                                      int direction) {
        if (instance().conf.y) {
            System.err.println(format("      Matching REF %s %s %s", sequence, position, direction));
        }
        int mm = 0;
        for (int n = 0; n < sequence.length(); n++) {
            final Character refCh = reference.get(position + direction * n);
            if (refCh == null || charAt(sequence, direction == 1 ? n : direction * n - 1) != refCh) {
                mm++;
            }
        }
        return mm <= 3 && mm / (double)sequence.length() < 0.15;
    }

    static String getMateReferenceName(SAMRecord record) {
        // TODO: htsjdk has "*" by default
        if (record.getMateReferenceName() == null) {
            return "*";
        }

        if (record.getReferenceName().equals(record.getMateReferenceName())) {
            return "=";
        }
        return record.getMateReferenceName();
    }
}
