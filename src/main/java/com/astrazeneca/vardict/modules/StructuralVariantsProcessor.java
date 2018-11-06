package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.variations.Cluster;
import com.astrazeneca.vardict.variations.Mate;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.SequenceUtil;

import java.io.IOException;
import java.util.*;

import static com.astrazeneca.vardict.collection.VariationMap.getSV;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.ReferenceResource.isLoaded;
import static com.astrazeneca.vardict.collection.VariationMap.SV;
import static com.astrazeneca.vardict.modules.SAMFileParser.getMrnm;
import static com.astrazeneca.vardict.modules.SAMFileParser.ismatchref;
import static com.astrazeneca.vardict.modules.SAMFileParser.parseSAM;
import static com.astrazeneca.vardict.modules.ToVarsBuilder.CURSEG;
import static com.astrazeneca.vardict.modules.VariationRealigner.*;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.Utils.*;

import static java.lang.Math.abs;
import static java.util.Collections.reverseOrder;
import static java.util.Comparator.comparing;
import static java.util.stream.Collectors.toList;

public class StructuralVariantsProcessor {
    private static ThreadLocal<Map<Integer, List<Sclip>>> SOFTP2SV = ThreadLocal.withInitial(HashMap::new);

    /**
     * Prepare and fill SV structures for deletions, duplications, inversions.
     * @param rlen max read length
     * @param conf Configuration
     * @param svStructures Contains all structural variants ends for deletions, duplications, inversions, insertions
     and structures for structural variant analysis
     * @param record parsed SAMRecord
     * @param queryQuality string of query qualities from record
     * @param nm number of mismatches from record
     * @param dir direction of the strand
     * @param cigar cigar from record (modified by realignment)
     * @param mdir direction of mate strand
     * @param start adjusted start position of read
     * @param rlen2 the total length, including soft-clipped bases
     */
    public static void prepareSVStructuresForAnalysis(int rlen,
                                                      Configuration conf,
                                                      SVStructures svStructures,
                                                      SAMRecord record,
                                                      String queryQuality,
                                                      int nm,
                                                      boolean dir,
                                                      Cigar cigar,
                                                      boolean mdir,
                                                      int start,
                                                      int rlen2) {
        int mstart = record.getMateAlignmentStart();
        int mend = mstart + rlen2;
        int end = start;
        List<String> msegs = globalFind(ALIGNED_LENGTH_MND, cigar.toString());
        end += sum(msegs);

        int soft5 = 0;

        jregex.Matcher matcher = BEGIN_NUM_S_OR_BEGIN_NUM_H.matcher(cigar.toString());
        if (matcher.find()) {
            int tt = toInt(matcher.group(1));
            if (tt != 0 && queryQuality.charAt(tt - 1) - 33 > conf.goodq) {
                soft5 = start;
            }
        }
        int soft3 = 0;

        matcher = END_NUM_S_OR_NUM_H.matcher(cigar.toString());
        if (matcher.find()) {
            int tt = toInt(matcher.group(1));
            if (tt != 0 && queryQuality.charAt(record.getReadLength() - tt) - 33 > conf.goodq) {
                soft3 = end;
            }
        }
        int dirNum = dir ? -1 : 1;
        int mdirNum = mdir ? 1 : -1;
        int MIN_D = 75;
        if (getMrnm(record).equals("=")) {
            int mlen = record.getInferredInsertSize();
            if (record.getStringAttribute(SAMTag.MC.name()) != null
                    && MC_Z_NUM_S_ANY_NUM_S.matcher(record.getStringAttribute(SAMTag.MC.name())).find()) {
                // Ignore those with mates mapped with softcliping at both ends
            } else if (record.getIntegerAttribute(SAMTag.MQ.name()) != null
                    && record.getIntegerAttribute(SAMTag.MQ.name()) < 15) {
                // Ignore those with mate mapping quality less than 15
            } else if (dirNum * mdirNum == -1 && (mlen * dirNum) > 0 ) {
                // deletion candidate
                mlen = mstart > start ? mend - start : end - mstart;
                if(abs(mlen) > conf.INSSIZE + conf.INSSTDAMT * conf.INSSTD ) {
                    if (dirNum == 1) {
                        if (svStructures.svfdel.size() == 0
                                || start - svStructures.svdelfend > Configuration.MINSVCDIST * rlen) {
                            Sclip sclip = new Sclip();
                            sclip.cnt = 0;
                            svStructures.svfdel.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svfdel), start, end, mstart,
                                mend, dirNum, rlen2, mlen, soft3, rlen/2.0,
                                queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                        svStructures.svdelfend = end;
                    } else {
                        if (svStructures.svrdel.size() == 0
                                || start - svStructures.svdelrend > Configuration.MINSVCDIST * rlen) {
                            Sclip sclip = new Sclip();
                            sclip.cnt = 0;
                            svStructures.svrdel.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svrdel), start, end, mstart,
                                mend, dirNum, rlen2, mlen, soft5, rlen/2.0,
                                queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                        svStructures.svdelrend = end;
                    }

                    if (!svStructures.svfdel.isEmpty()
                            && abs(start - svStructures.svdelfend) <= Configuration.MINSVCDIST * rlen ){
                        adddisccnt(getLastSVStructure(svStructures.svfdel));
                    }
                    if (!svStructures.svrdel.isEmpty()
                            && abs(start - svStructures.svdelrend) <= Configuration.MINSVCDIST * rlen ){
                        adddisccnt(getLastSVStructure(svStructures.svrdel));
                    }
                    if (!svStructures.svfdup.isEmpty() && abs(start - svStructures.svdupfend) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svfdup));
                    }
                    if (!svStructures.svrdup.isEmpty() && abs(start - svStructures.svduprend) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svrdup));
                    }
                    if (!svStructures.svfinv5.isEmpty() && abs(start - svStructures.svinvfend5) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svfinv5));
                    }
                    if (!svStructures.svrinv5.isEmpty() && abs(start - svStructures.svinvrend5) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svrinv5));
                    }
                    if (!svStructures.svfinv3.isEmpty() && abs(start - svStructures.svinvfend3) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svfinv3));
                    }
                    if (!svStructures.svrinv3.isEmpty() && abs(start - svStructures.svinvrend3) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svrinv3));
                    }
                }
            } else if (dirNum * mdirNum == -1 && dirNum * mlen < 0) {
                //duplication
                if (dirNum == 1) {
                    if (svStructures.svfdup.size() == 0
                            || start - svStructures.svdupfend > Configuration.MINSVCDIST * rlen) {
                        Sclip sclip = new Sclip();
                        sclip.cnt = 0;
                        svStructures.svfdup.add(sclip);
                    }
                    addSV(getLastSVStructure(svStructures.svfdup), start, end, mstart,
                            mend, dirNum, rlen2, mlen, soft3, rlen/2.0,
                            queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                    svStructures.svdupfend = end;
                } else {
                    if (svStructures.svrdup.size() == 0
                            || start - svStructures.svduprend > Configuration.MINSVCDIST * rlen) {
                        Sclip sclip = new Sclip();
                        sclip.cnt = 0;
                        svStructures.svrdup.add(sclip);
                    }
                    addSV(getLastSVStructure(svStructures.svrdup), start, end, mstart,
                            mend, dirNum, rlen2, mlen, soft5, rlen/2.0,
                            queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                    svStructures.svduprend = end;
                }
                if (!svStructures.svfdup.isEmpty()
                        && abs(start - svStructures.svdupfend) <= Configuration.MINSVCDIST * rlen ) {
                    getLastSVStructure(svStructures.svfdup).disc++;
                }
                if (!svStructures.svrdup.isEmpty()
                        && abs(start - svStructures.svduprend) <= Configuration.MINSVCDIST * rlen ) {
                    getLastSVStructure(svStructures.svrdup).disc++;
                }
                if (!svStructures.svfdel.isEmpty() && abs(start - svStructures.svdelfend) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svfdel));
                }
                if (!svStructures.svrdel.isEmpty() && abs(start - svStructures.svdelrend) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svrdel));
                }
                if (!svStructures.svfinv5.isEmpty() && abs(start - svStructures.svinvfend5) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svfinv5));
                }
                if (!svStructures.svrinv5.isEmpty() && abs(start - svStructures.svinvrend5) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svrinv5));
                }
                if (!svStructures.svfinv3.isEmpty() && abs(start - svStructures.svinvfend3) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svfinv3));
                }
                if (!svStructures.svrinv3.isEmpty() && abs(start - svStructures.svinvrend3) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svrinv3));
                }
            } else if (dirNum * mdirNum == 1) { // Inversion
                if (dirNum == 1 && mlen != 0 ) {
                    if (mlen < -3 * rlen) {
                        if (svStructures.svfinv3.size() == 0
                                || start - svStructures.svinvfend3 > Configuration.MINSVCDIST * rlen) {
                            Sclip sclip = new Sclip();
                            sclip.cnt = 0;
                            svStructures.svfinv3.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svfinv3), start, end, mstart,
                                mend, dirNum, rlen2, mlen, soft3, rlen/2.0,
                                queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                        svStructures.svinvfend3 = end;
                        getLastSVStructure(svStructures.svfinv3).disc++;
                    } else if (mlen > 3 * rlen) {
                        if (svStructures.svfinv5.size() == 0
                                || start - svStructures.svinvfend5 > Configuration.MINSVCDIST * rlen) {
                            Sclip sclip = new Sclip();
                            sclip.cnt = 0;
                            svStructures.svfinv5.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svfinv5), start, end, mstart,
                                mend, dirNum, rlen2, mlen, soft3, rlen/2.0,
                                queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                        svStructures.svinvfend5 = end;
                        getLastSVStructure(svStructures.svfinv5).disc++;
                    }
                } else if (mlen != 0) {
                    if (mlen < -3 * rlen) {
                        if (svStructures.svrinv3.size() == 0
                                || start - svStructures.svinvrend3 > Configuration.MINSVCDIST * rlen) {
                            Sclip sclip = new Sclip();
                            sclip.cnt = 0;
                            svStructures.svrinv3.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svrinv3), start, end, mstart,
                                mend, dirNum, rlen2, mlen, soft5, rlen/2.0,
                                queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                        svStructures.svinvrend3 = end;
                        getLastSVStructure(svStructures.svrinv3).disc++;

                    } else if (mlen > 3 * rlen) {
                        if (svStructures.svrinv5.size() == 0
                                || start - svStructures.svinvrend5 > Configuration.MINSVCDIST * rlen) {
                            Sclip sclip = new Sclip();
                            sclip.cnt = 0;
                            svStructures.svrinv5.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svrinv5), start, end, mstart,
                                mend, dirNum, rlen2, mlen, soft5, rlen/2.0,
                                queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                        svStructures.svinvrend5 = end;
                        getLastSVStructure(svStructures.svrinv5).disc++;
                    }
                }
                if (mlen != 0) {
                    if (!svStructures.svfdel.isEmpty() && (start - svStructures.svdelfend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svfdel));
                    }
                    if (!svStructures.svrdel.isEmpty() && (start - svStructures.svdelrend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svrdel));
                    }
                    if (!svStructures.svfdup.isEmpty() && (start - svStructures.svdupfend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svfdup));
                    }
                    if (!svStructures.svrdup.isEmpty() && (start - svStructures.svduprend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svrdup));
                    }
                }
            }
        } else { // Inter-chr translocation
            // to be implemented
            String mchr = getMrnm(record);
            if (record.getStringAttribute(SAMTag.MC.name()) != null
                    && MC_Z_NUM_S_ANY_NUM_S.matcher(record.getStringAttribute(SAMTag.MC.name())).find()) {
                // Ignore those with mates mapped with softcliping at both ends
            } else if (record.getIntegerAttribute(SAMTag.MQ.name()) != null
                    && record.getIntegerAttribute(SAMTag.MQ.name()) < 15) {
                // Ignore those with mate mapping quality less than 15
            } else if (dirNum == 1) {
                if (svStructures.svffus.get(mchr) == null
                        || start - svStructures.svfusfend.get(mchr) > Configuration.MINSVCDIST * rlen) {
                    Sclip sclip = new Sclip();
                    sclip.cnt = 0;
                    List<Sclip> sclips = svStructures.svffus.getOrDefault(mchr, new ArrayList<>());
                    sclips.add(sclip);
                    svStructures.svffus.put(mchr, sclips);
                }
                int svn = svStructures.svffus.get(mchr).size() - 1;
                addSV(svStructures.svffus.get(mchr).get(svn), start, end, mstart,
                        mend, dirNum, rlen2, 0, soft3, rlen/2.0,
                        queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                svStructures.svfusfend.put(mchr, end);
                svStructures.svffus.get(mchr).get(svn).disc++;
            } else {
                if (svStructures.svrfus.get(mchr) == null
                        || start - svStructures.svfusrend.get(mchr) > Configuration.MINSVCDIST * rlen) {
                    Sclip sclip = new Sclip();
                    sclip.cnt = 0;
                    List<Sclip> sclips = svStructures.svrfus.getOrDefault(mchr, new ArrayList<>());
                    sclips.add(sclip);
                    svStructures.svrfus.put(mchr, sclips);
                }
                int svn = svStructures.svrfus.get(mchr).size() - 1;
                addSV(svStructures.svrfus.get(mchr).get(svn), start, end, mstart,
                        mend, dirNum, rlen2, 0, soft5, rlen/2.0,
                        queryQuality.charAt(15) - 33, record.getMappingQuality(), nm, conf);
                svStructures.svfusrend.put(mchr, end);
                svStructures.svrfus.get(mchr).get(svn).disc++;
            }
            //TODO: There are no abs before (start - int) in perl, but otherwise here will be negative numbers
            //adddisccnt( $svrdel[$#svrdel] ) if ( @svrdel && ($start - $svdelrend) <= 25 );
            if (!svStructures.svfdel.isEmpty() && start - svStructures.svdelfend <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svfdel));
            }
            if (!svStructures.svrdel.isEmpty() && start - svStructures.svdelrend <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svrdel));
            }
            if (!svStructures.svfdup.isEmpty() && start - svStructures.svdupfend <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svfdup));
            }
            if (!svStructures.svrdup.isEmpty() && start - svStructures.svduprend <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svrdup));
            }
            if (!svStructures.svfinv5.isEmpty() && start - svStructures.svinvfend5 <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svfinv5));
            }
            if (!svStructures.svrinv5.isEmpty() && start - svStructures.svinvrend5 <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svrinv5));
            }
            if (!svStructures.svfinv3.isEmpty() && start - svStructures.svinvfend3 <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svfinv3));
            }
            if (!svStructures.svrinv3.isEmpty() && start - svStructures.svinvrend3 <= 25) {
                adddisccnt(getLastSVStructure(svStructures.svrinv3));
            }
        }
     }

    private static Sclip getLastSVStructure(List<Sclip> svStructure) {
        return svStructure.get(svStructure.size() - 1);
    }

    /**
     * Find SVs for each structural variant structure
     * @param region region of interest
     * @param bam BAM file name
     * @param chrs map of chromosome lengths
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param rlen max read length
     * @param reference reference in a given region
     * @param conf Configuration
     * @param hash map of variants on positions
     * @param iHash map of insertions on positions
     * @param cov map of coverage on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param svStructures Contains all structural variants ends for deletions, duplications, inversions, insertions
    and structures for structural variant analysis
     * @param referenceResource object for access to reference map
     * @throws IOException if BAM file can't be read
     */
    static void findAllSVs(Region region, String bam, Map<String, Integer> chrs, String sample,
                                  Set<String> splice, String ampliconBasedCalling, int rlen, Reference reference,
                                  Configuration conf, Map<Integer, VariationMap<String, Variation>> hash,
                                  Map<Integer, VariationMap<String, Variation>> iHash, Map<Integer, Integer> cov,
                                  Map<Integer, Sclip> sclip3, Map<Integer, Sclip> sclip5,
                                  SVStructures svStructures, ReferenceResource referenceResource) throws IOException {
        if (conf.y) {
            System.err.println("Start Structural Variants: DEL\n");
        }
        findDEL(conf, hash, iHash, cov, sclip5, sclip3, reference, referenceResource, region, bam, svStructures.svfdel,
                svStructures.svrdel, rlen, chrs, sample, splice, ampliconBasedCalling);
        if (conf.y) {
            System.err.println("Start Structural Variants: INV\n");
        }
        findINV(conf, hash, iHash, cov, sclip5, sclip3, reference, referenceResource, region, bam, svStructures.svfinv3,
                svStructures.svrinv3, svStructures.svfinv5, svStructures.svrinv5,
                rlen, chrs, sample, splice, ampliconBasedCalling);
        if (conf.y) {
            System.err.println("Start Structural Variants\n");
        }
        findsv(hash, cov, sclip5, sclip3, reference, region, svStructures.svfdel, svStructures.svrdel, conf, rlen);
        if (conf.y) {
            System.err.println("Start Structural Variants: DEL discordant pairs only\n");
        }
        findDELdisc(conf, hash, cov, sclip5, sclip3, reference, referenceResource, region, svStructures.svfdel,
                svStructures.svrdel, splice, rlen, chrs) ;
        if (conf.y) {
            System.err.println("Start Structural Variants: INV discordant pairs only\n");
        }
        findINVdisc(hash, cov, sclip5, sclip3, reference, referenceResource, region, svStructures.svfinv3, svStructures.svrinv3,
                svStructures.svfinv5, svStructures.svrinv5, conf, rlen, chrs);
        if (conf.y) {
            System.err.println("Start Structural Variants: DUP discordant pairs only\n");
        }
        findDUPdisc(hash, iHash, cov, sclip5, sclip3, reference, referenceResource, region, bam, sample, splice,
                ampliconBasedCalling, svStructures.svfdup, svStructures.svrdup, conf, rlen, chrs);
    }

    /**
     * Add structural variant to current structural variant structure and update
     * @param sdref current structural variant structure
     */
    private static void addSV (Sclip sdref,
                               int start_s,
                               int end_e,
                               int mateStart_ms,
                               int mateEnd_me,
                               int dir,
                               int rlen,
                               int mlen,
                               int softp,
                               double pmean_rp,
                               double qmean,
                               double Qmean,
                               double nm,
                               Configuration conf) {
        sdref.cnt++;
        sdref.incDir(dir == 1 ? false : true);

        if (qmean >= conf.goodq) {
            sdref.hicnt++;
        } else {
            sdref.locnt++;
        }

        if (sdref.start == 0 || sdref.start >= start_s) {
            sdref.start = start_s;
        }
        if (sdref.end == 0 || sdref.end <= end_e) {
            sdref.end = end_e;
        }
        sdref.mates.add(new Mate(mateStart_ms, mateEnd_me, mlen,
                start_s, end_e, pmean_rp, qmean, Qmean, nm));

        if (sdref.mstart == 0 || sdref.mstart >= mateStart_ms) {
            sdref.mstart = mateStart_ms;
        }
        if (sdref.mend == 0 || sdref.mend <= mateEnd_me ) {
            sdref.mend = mateStart_ms + rlen;
        }
        if (softp != 0) {
            if (dir == 1) {
                if (abs(softp - sdref.end) < 10 ) {
                    int softCount = sdref.soft.getOrDefault(softp, 0);
                    softCount++;
                    sdref.soft.put(softp, softCount);
                }
            } else {
                if (abs(softp - sdref.start) < 10 ) {
                    int softCount = sdref.soft.getOrDefault(softp, 0);
                    softCount++;
                    sdref.soft.put(softp, softCount);
                }
            }
        }
    }

    /**
     * Filter all structural variants structures
     * @param rlen max read length
     * @param conf Configuration
     * @param svStructures Contains all structural variants ends for deletions, duplications, inversions, insertions
    and structures for structural variant analysis
     */
     static void filterAllSVStructures(int rlen,
                                             Configuration conf,
                                             SVStructures svStructures) {
        filterSV(svStructures.svfinv3, rlen, conf);
        filterSV(svStructures.svrinv3, rlen, conf);
        filterSV(svStructures.svfinv5, rlen, conf);
        filterSV(svStructures.svrinv5, rlen, conf);
        filterSV(svStructures.svfdel, rlen, conf);
        filterSV(svStructures.svrdel, rlen, conf);
        filterSV(svStructures.svfdup, rlen, conf);
        filterSV(svStructures.svrdup, rlen, conf);

        for (Map.Entry<String, List<Sclip>> svv : svStructures.svffus.entrySet()) {
            filterSV(svv.getValue(), rlen, conf);
        }
        for (Map.Entry<String, List<Sclip>> svv : svStructures.svrfus.entrySet()) {
            filterSV(svv.getValue(), rlen, conf);
        }
        for (Map.Entry<Integer, List<Sclip>> entry : SOFTP2SV.get().entrySet()) {
            Integer key = entry.getKey();
            List<Sclip> sclips = entry.getValue();
            sclips.sort(comparing((Sclip sclip) -> sclip.cnt).reversed());
            SOFTP2SV.get().put(key, sclips);
        }
    }

    /**
     * Filter SV  by checking created clusters
     * @param svList_sva contains list of SVs to filter
     * @param rlen read length
     * @param conf configuration of Vardict
     */
    static void filterSV(List<Sclip> svList_sva,
                         int rlen,
                         Configuration conf) {
        for (Sclip sv: svList_sva) {
            Cluster cluster = checkCluster(sv.mates, rlen, conf);

            if (cluster.mateStart_ms != 0) {
                sv.mstart = cluster.mateStart_ms;
                sv.mend = cluster.mateEnd_me;
                sv.cnt = cluster.cnt;
                sv.mlen = cluster.mateLength_mlen;
                sv.start = cluster.start_s;
                sv.end = cluster.end_e;
                sv.pmean = cluster.pmean_rp;
                sv.qmean = cluster.qmean_q;
                sv.Qmean = cluster.Qmean_Q;
                sv.nm = cluster.nm;
            } else {
                sv.used = true;
            }
            // Too many unhappy mates are false positive
            if (sv.disc != 0 && sv.cnt/(double) sv.disc < 0.5) {
                if (!(sv.cnt/(double) sv.disc >= 0.35 && sv.cnt >= 5)) {
                    sv.used = true;
                }
            }

            List<Map.Entry<Integer, Integer>> soft = new ArrayList<>(sv.soft.entrySet());
            soft.sort(comparing((Map.Entry<Integer, Integer> entry) -> entry.getValue(), Integer::compareTo).reversed());

            sv.softp = soft.size() > 0 ? soft.get(0).getKey() : 0;
            if (sv.softp != 0) {
                List<Sclip> sclips = SOFTP2SV.get().getOrDefault(sv.softp, new ArrayList<>());
                sclips.add(sv);
                SOFTP2SV.get().put(sv.softp, sclips);
            }

            if (conf.y) {
                System.err.printf("SV cluster: %s %s %s %s Cnt: %s Discordant Cnt: %s Softp: %s Used: %s\n",
                        cluster.start_s, cluster.end_e,
                        cluster.mateStart_ms, cluster.mateEnd_me, sv.cnt, sv.disc, sv.softp, sv.used);
            }
        }
    }

    /**
     * Check cluster for start and end of SV mates
     * @param mates list of mates for current SV
     * @param rlen read length
     * @param conf configuration of Vardict
     * @return Cluster with updated starts and ends
     */
    static Cluster checkCluster(List<Mate> mates,
                                int rlen,
                                Configuration conf) {
        mates.sort(comparing(mate -> mate.mateStart_ms, Integer::compareTo));

        List<Cluster> clusters = new ArrayList<>();
        Mate firstMate = mates.get(0);
        clusters.add(new Cluster(0, firstMate.mateStart_ms, firstMate.mateEnd_me, firstMate.start_s, firstMate.end_e));

        int cur = 0;
        for (Mate mate_m : mates) {
            Cluster currentCluster = clusters.get(cur);
            if (mate_m.mateStart_ms - currentCluster.mateEnd_me > Configuration.MINSVCDIST * rlen) {
                cur++;
                clusters.add(cur, new Cluster(0, mate_m.mateStart_ms, mate_m.mateEnd_me, mate_m.start_s, mate_m.end_e));
                currentCluster = clusters.get(cur);
            }

            currentCluster.cnt++;
            currentCluster.mateLength_mlen += mate_m.mateLength_mlen;

            if (mate_m.mateEnd_me > currentCluster.mateEnd_me) {
                currentCluster.mateEnd_me = mate_m.mateEnd_me;
            }
            if (mate_m.start_s < currentCluster.start_s) {
                currentCluster.start_s = mate_m.start_s;
            }
            if (mate_m.end_e > currentCluster.end_e) {
                currentCluster.end_e = mate_m.end_e;
            }
            currentCluster.pmean_rp += mate_m.pmean_rp;
            currentCluster.qmean_q += mate_m.qmean_q;
            currentCluster.Qmean_Q += mate_m.Qmean_Q;
            currentCluster.nm += mate_m.nm;
        }
        clusters.sort(comparing((Cluster cluster) -> cluster.cnt).reversed());

        if (conf.y) {
            System.err.print("Clusters; ");
            clusters.forEach(cluster -> System.err.print(join("; ", cluster.cnt, cluster.start_s, cluster.end_e,
                    cluster.mateStart_ms, cluster.mateEnd_me)));
            System.err.println(join("; ", "; out of", mates.size()));
        }
        Cluster firstCluster = clusters.get(0);
        return firstCluster.cnt / (double) mates.size() >= 0.60
                ? new Cluster(firstCluster.mateStart_ms, firstCluster.mateEnd_me, firstCluster.cnt,
                firstCluster.mateLength_mlen/firstCluster.cnt, firstCluster.start_s, firstCluster.end_e,
                firstCluster.pmean_rp, firstCluster.qmean_q, firstCluster.Qmean_Q, firstCluster.nm)
                : new Cluster(0,0,0,0,0,0,0,0.0,0,0);
    }

    /**
     * Mark SV clusters as used
     * @param start start of the region $s
     * @param end end of the region $e
     * @param structuralVariants_sv contains list of list of Structural variants $sv
     * @param rlen read length
     * @param conf configuration of Vardict
     * @return tuple of coverage, count of SVs overlaping with mates, number of pairs
     */
    public static Tuple.Tuple3<Integer, Integer, Integer> markSV(int start,
                                                          int end,
                                                          List<List<Sclip>> structuralVariants_sv,
                                                          int rlen,
                                                          Configuration conf) {
        int cov = 0;
        int pairs = 0;
        int cnt = 0;

        for (List<Sclip> currentSclips_sr : structuralVariants_sv) {
            for (Sclip sv_r : currentSclips_sr) {
                Tuple.Tuple2<Integer, Integer> tuple = sv_r.start < sv_r.mstart
                        ? new Tuple.Tuple2<>(sv_r.end, sv_r.mstart)
                        : new Tuple.Tuple2<>(sv_r.mend, sv_r.start);
                if (conf.y) {
                    System.err.printf("   Marking SV %s %s %s %s cnt: %s\n", start, end, tuple._1, tuple._2, sv_r.cnt);
                }
                if (isOverlap(start, end, tuple, rlen) ) {
                    if (conf.y) {
                        System.err.printf("       SV %s %s %s %s cnt: %s marked\n", start, end, tuple._1, tuple._2, sv_r.cnt);
                    }
                    sv_r.used = true;
                    cnt++;
                    pairs += sv_r.cnt;
                    cov += (int) ((sv_r.cnt * rlen)/(sv_r.end - sv_r.start)) + 1;
                }
            }
        }
        return new Tuple.Tuple3<>(cov, cnt, pairs);
    }

    /**
     * Mark DUP clusters as used
     * @param start start of the region $s
     * @param end end of the region $e
     * @param structuralVariants_sv contains list of list of Structural variants $sv
     * @param rlen read length
     * @param conf  configuration of Vardict
     * @return tuple of count of SVs overlaping with mates, number of pairs
     */
    static Tuple.Tuple2<Integer, Integer> markDUPSV(int start,
                                                    int end,
                                                    List<List<Sclip>> structuralVariants_sv,
                                                    int rlen,
                                                    Configuration conf) {
        int cov = 0;
        int pairs = 0;
        int cnt = 0;

        for (List<Sclip> currentSclips_sr : structuralVariants_sv) {
            for (Sclip sv_r : currentSclips_sr) {
                Tuple.Tuple2<Integer, Integer> tuple = sv_r.start < sv_r.mstart
                        ? new Tuple.Tuple2<>(sv_r.start, sv_r.mend)
                        : new Tuple.Tuple2<>(sv_r.mstart, sv_r.end);
                if (conf.y) {
                    System.err.printf("   Marking DUP SV %s %s %s %s cnt: %s\n", start, end, tuple._1, tuple._2, sv_r.cnt);
                }
                if (isOverlap(start, end, tuple, rlen)) {
                    if (conf.y) {
                        System.err.printf("       DUP SV %s %s %s %s cnt: %s marked\n", start, end, tuple._1, tuple._2, sv_r.cnt);
                    }
                    sv_r.used = true;
                    cnt++;
                    pairs += sv_r.cnt;
                    cov += (int) ((sv_r.cnt * rlen)/(sv_r.end - sv_r.start)) + 1;
                }
            }
        }
        return new Tuple.Tuple2<>(cnt, pairs);
    }

    /**
     * Determine overlapping
     * @param start1 start of the first SV $s1
     * @param end1 end of the first SV $e1
     * @param tuple2 contains start and end of second SV
     * @param rlen read length
     * @return true if SV overlaps with mate
     */
    static boolean isOverlap (int start1,
                              int end1,
                              Tuple.Tuple2<Integer, Integer> tuple2,
                              int rlen){
        int start2 = tuple2._1;
        int end2 = tuple2._2;
        if (start1 >= end2 || start2 >= end1) {
            return false;
        }
        List<Integer> positions = Arrays.asList(start1, end1, start2, end2);
        positions.sort(Integer::compareTo);

        int ins = positions.get(2) - positions.get(1);
        if (ins/(double)(end1 - start1) > 0.75 && ins/(double)(end2 - start2) > 0.75 ) {
            return true;
        }
        if (positions.get(1) - positions.get(0) + positions.get(3) - positions.get(2) < 3 * rlen) {
            return true;
        }
        return false;
    }

    /**
     * Given a candidate SV identified by clipped reads, check whether there're mates to support it
     * @param chr chromosome name
     * @param start start of the region
     * @param end end of the region
     * @param sv list of SVs in list of clusters
     * @param RLEN max read length
     * @param conf Configuration
     * @return tuple if pairs, mean position, mean base quality, mean mapping quality and number of mismatches
     */
    static Tuple.Tuple5<Integer, Double, Double, Double, Double> checkPairs(String chr,
                                                                               int start,
                                                                               int end,
                                                                               List<List<Sclip>> sv,
                                                                               int RLEN,
                                                                               Configuration conf) {
        int pairs = 0;
        double pmean = 0;
        double qmean = 0;
        double Qmean = 0;
        double nm = 0;

        Tuple.Tuple5<Integer, Double, Double, Double, Double> tuple5 = tuple(pairs, pmean, qmean, Qmean, nm);

        for (List<Sclip> svcluster: sv) {
            for(Sclip svr : svcluster) {
                if (svr.used) {
                    continue;
                }
                int s = (svr.start + svr.end) / 2;
                int e = (svr.mstart + svr.mend) / 2;
                if (s > e) {
                    int temp = s;
                    s = e;
                    e = temp;
                }
                if (!isOverlap(start, end, new Tuple.Tuple2<>(s, e), RLEN)) {
                    continue;
                }
                if (svr.cnt > pairs) {
                    tuple5 = tuple(svr.cnt, svr.pmean, svr.qmean, svr.Qmean, svr.nm);
                    pairs = svr.cnt;
                }
                svr.used = true;
                if (conf.y) {
                    System.err.printf("      Pair [%s:%s-%s] overlapping [%s:%s-%s] found and marked.\n", chr, s, e, chr, start, end);
                }
            }
        }
        return tuple5;
    }

    /**
     * Add discordant count
     * @param svref
     */
    private static void adddisccnt(Sclip svref) {
        svref.disc++;
    }

    /**
     * Find candidate SVs on 3' and 5' ends
     * @param hash map of variants on positions
     * @param cov map of coverage on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param reference reference in a given region
     * @param region region of interest
     * @param svfdel list of DEL SVs in forward direction
     * @param svrdel list of DEL SVs in reverse direction
     * @param conf Configuration
     * @param RLEN max read length
     */
    static void findsv(Map<Integer, VariationMap<String, Variation>> hash,
                       Map<Integer, Integer> cov,
                       Map<Integer, Sclip> sclip5,
                       Map<Integer, Sclip> sclip3,
                       Reference reference,
                       Region region,
                       List<Sclip> svfdel,
                       List<Sclip> svrdel,
                       Configuration conf,
                       int RLEN) {
        Map<Integer, Character> ref = reference.referenceSequences;
        List<Tuple.Tuple3<Integer, Sclip, Integer>> tmp5 = new ArrayList<>();
        for(Map.Entry<Integer, Sclip> entry : sclip5.entrySet()) {
            int p = entry.getKey();
            Sclip sc5v = entry.getValue();
            if (sc5v.used) {
                continue;
            }
            if (p < CURSEG.get()._2 || p > CURSEG.get()._3) {
                continue;
            }
            tmp5.add(new Tuple.Tuple3<>(p, sc5v, sc5v.cnt));
        }
        tmp5.sort(comparing((Tuple.Tuple3<Integer, Sclip, Integer> tuple) -> tuple._3).reversed());

        List<Tuple.Tuple3<Integer, Sclip, Integer>> tmp3 = new ArrayList<>();
        for(Map.Entry<Integer, Sclip> entry : sclip3.entrySet()) {
            int p = entry.getKey();
            Sclip sc3v = entry.getValue();
            if (sc3v.used) {
                continue;
            }
            if (p < CURSEG.get()._2 || p > CURSEG.get()._3) {
                continue;
            }
            tmp3.add(new Tuple.Tuple3<>(p, sc3v, sc3v.cnt));
        }
        tmp3.sort(comparing((Tuple.Tuple3<Integer, Sclip, Integer> tuple) -> tuple._3).reversed());

        for (Tuple.Tuple3<Integer, Sclip, Integer> tuple5 : tmp5) {
            int p5 = tuple5._1;
            Sclip sc5v = tuple5._2;
            int cnt5 = tuple5._3;
            if (cnt5 < conf.minr) {
                break;
            }
            if (sc5v.used) {
                continue;
            }
            if (SOFTP2SV.get().containsKey(p5) && SOFTP2SV.get().get(p5).get(0).used) {
                continue;
            }
            String seq = findconseq(sc5v, conf, 0);
            if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                continue;
            }
            if (conf.y) {
                System.err.printf("  Finding SV 5': %s %s cnt: %s\n", seq, p5, cnt5);
            }
            Tuple.Tuple2<Integer, String> tuple2 = findMatch(conf, seq, reference, p5, -1);
            int bp = tuple2._1;
            String EXTRA = tuple2._2;

            if (bp != 0) {
                // candidate deletion
                if (bp < p5) {
                    Tuple.Tuple5<Integer, Double, Double, Double, Double> tuplePairs = checkPairs(region.chr, bp, p5, Arrays.asList(svfdel, svrdel), RLEN, conf);
                    int pairs = tuplePairs._1;
                    double pmean = tuplePairs._2;
                    double qmean = tuplePairs._3;
                    double Qmean = tuplePairs._4;
                    double nm = tuplePairs._5;
                    if (pairs == 0) {
                        continue;
                    }
                    p5--;
                    bp++;
                    int dellen = p5 - bp + 1;

                    final Variation vref = getVariation(hash, bp, "-" + dellen);
                    vref.cnt = 0;

                    SV sv = getSV(hash, bp);
                    sv.type = "DEL";
                    sv.pairs += pairs;
                    sv.splits += cnt5;
                    sv.clusters += pairs != 0 ? 1 : 0;

                    if (!cov.containsKey(bp)) {
                        cov.put(bp, pairs + sc5v.cnt);
                    }
                    if (cov.containsKey(p5 + 1) && cov.get(bp) < cov.get(p5 + 1)) {
                        cov.put(bp, cov.get(p5 + 1));
                    }
                    adjCnt(vref, sc5v, conf);
                    Variation tmp = new Variation();
                    tmp.cnt = pairs;
                    tmp.hicnt = pairs;
                    tmp.dirPlus = (pairs / 2);
                    tmp.dirMinus = pairs - (pairs / 2);
                    tmp.pmean = pmean;
                    tmp.qmean = qmean;
                    tmp.Qmean = Qmean;
                    tmp.nm = nm;
                    adjCnt(vref, tmp, conf);
                    if (conf.y) {
                        System.err.println("    Finding candidate deletion 5'");
                    }
                } else { // candidate duplication
                }
            } else { // candidate inversion
                Tuple.Tuple2<Integer, String> tuple2Rev = findMatchRev(conf, seq, reference, p5, -1);
                bp = tuple2Rev._1;
                EXTRA = tuple2Rev._2;
                if (bp == 0) {
                    continue;
                }
                if (!(abs(bp - p5) > Configuration.SVFLANK)) {
                    continue;
                }
                if (bp > p5) { // bp at 3' side
                } else { // bp at 5' side
                    int temp = bp;
                    bp = p5;
                    p5 = temp;
                }
                bp--;
                while (ref.containsKey(bp + 1)
                        && isHasAndEquals(complement(ref.get(bp + 1)), ref, p5 - 1)) {
                    p5--;
                    if (p5 != 0) {
                        bp++;
                    }
                }

                String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp - Configuration.SVFLANK + 1, bp));
                String ins3 = SequenceUtil.reverseComplement(joinRef(ref, p5, p5 + Configuration.SVFLANK - 1));
                int mid = bp - p5 - ins5.length() - ins3.length() + 1;

                String vn = "-" + (bp - p5 + 1) + "^" + ins5 + "<inv" + mid + ">" + ins3 + EXTRA;
                if (mid <= 0) {
                    String tins = SequenceUtil.reverseComplement(joinRef(ref, p5, bp));
                    vn = "-" + (bp - p5 + 1) + "^" + tins + EXTRA;
                }

                final Variation vref = getVariation(hash, p5, vn);
                // TODO: 26/07/18 Seem that next statement is redundant since it assigns cnt to 0 if it is already 0
                //	    $hash->{ $p5 }->{ $vn }->{ cnt } = 0 unless( $hash->{ $p5 }->{ $vn }->{ cnt } );

                SV sv = getSV(hash, p5);
                sv.type = "INV";
                sv.pairs += 0;
                sv.splits += cnt5;
                sv.clusters += 0;

                adjCnt(vref, sc5v, conf);
                incCnt(cov, p5, cnt5);
                if (cov.containsKey(bp) && cov.get(p5) < cov.get(bp)) {
                    cov.put(p5, cov.get(bp));
                }
                if (conf.y) {
                    System.err.printf("    Found INV: %d %s Cnt:%d\n", p5, vn, cnt5);
                }
            }
        }
        for (Tuple.Tuple3<Integer, Sclip, Integer> tuple3 : tmp3) {
            int p3 = tuple3._1;
            Sclip sc3v = tuple3._2;
            int cnt3 = tuple3._3;
            if (cnt3 < conf.minr) {
                break;
            }
            if (sc3v.used) {
                continue;
            }
            if (SOFTP2SV.get().containsKey(p3) && SOFTP2SV.get().get(p3).get(0).used) {
                continue;
            }
            String seq = findconseq(sc3v, conf, 0);
            if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                continue;
            }
            if (conf.y) {
                System.err.printf("  Finding SV 3': %s %s cnt: %s\n", seq, p3, cnt3);
            }

            Tuple.Tuple2<Integer, String> tuple2 = findMatch(conf, seq, reference, p3, 1);
            int bp = tuple2._1;
            String EXTRA = tuple2._2;
            if (bp != 0) {
                if (bp > p3) { // candidate deletion
                    Tuple.Tuple5<Integer, Double, Double, Double, Double> tuplePairs = checkPairs(region.chr, p3, bp, Arrays.asList(svfdel, svrdel), RLEN, conf);
                    int pairs = tuplePairs._1;
                    double pmean = tuplePairs._2;
                    double qmean = tuplePairs._3;
                    double Qmean = tuplePairs._4;
                    double nm = tuplePairs._5;
                    if (pairs == 0) {
                        continue;
                    }
                    int dellen = bp - p3;
                    bp--;

                    while (isHasAndEquals(bp, ref, p3 - 1)) {
                        bp--;
                        if (bp != 0) {
                            p3--;
                        }
                    }

                    final Variation vref = getVariation(hash, p3, "-" + dellen);
                    vref.cnt = 0;

                    SV sv = getSV(hash, p3);
                    sv.type = "DEL";
                    sv.pairs += pairs;
                    sv.splits += cnt3;
                    sv.clusters += pairs != 0 ? 1 : 0;

                    if (!cov.containsKey(p3)) {
                        cov.put(p3, pairs + sc3v.cnt);
                    }
                    if (cov.containsKey(bp) && cov.get(bp) < cov.get(p3)) {
                        cov.put(bp, cov.get(p3));
                    }

                    adjCnt(vref, sc3v, conf);
                    Variation tmp = new Variation();
                    tmp.cnt = pairs;
                    tmp.hicnt = pairs;
                    tmp.dirPlus = (int) (pairs / 2);
                    tmp.dirMinus = pairs - (int) (pairs / 2);
                    tmp.pmean = pmean;
                    tmp.qmean = qmean;
                    tmp.Qmean = Qmean;
                    tmp.nm = nm;
                    adjCnt(vref, tmp, conf);
                    if (conf.y) {
                        System.err.println("    Finding candidate deletion 3'");
                    }
                } else { // candidate duplication
                }
            } else { // candidate inversion
                Tuple.Tuple2<Integer, String> tuple2Rev = findMatchRev(conf, seq, reference, p3, 1);
                bp = tuple2Rev._1;
                EXTRA = tuple2Rev._2;

                if (bp == 0) {
                    continue;
                }
                if (abs(bp - p3) <= Configuration.SVFLANK) {
                    continue;
                }
                if (bp < p3) { // bp at 5' side
                    int tmp = bp;
                    bp = p3;
                    p3 = tmp;
                    p3++;
                    bp--;
                } else { // bp at 3' side
                }

                while (ref.containsKey(bp + 1)
                        && isHasAndEquals(complement(ref.get(bp + 1)), ref, p3 - 1)) {
                    p3--;
                    if (p3 != 0) {
                        bp++;
                    }
                }
                String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp - Configuration.SVFLANK + 1, bp));
                String ins3 = SequenceUtil.reverseComplement(joinRef(ref, p3, p3 + Configuration.SVFLANK - 1));
                int mid = bp - p3 - 2 * Configuration.SVFLANK + 1;

                String vn = "-" + (bp - p3 + 1) + "^" + EXTRA + ins5 + "<inv" + mid + ">" + ins3;
                if (mid <= 0) {
                    String tins = SequenceUtil.reverseComplement(joinRef(ref, p3, bp));
                    vn = "-" + (bp - p3 + 1) + "^" + EXTRA + tins;
                }

                final Variation vref = getVariation(hash, p3, vn);
                // TODO: 26/07/18 Seem that next statement is redundant since it assigns cnt to 0 if it is already 0
                //$hash->{ $p3 }->{ $vn }->{ cnt } = 0 unless( $hash->{ $p3 }->{ $vn }->{ cnt } );

                SV sv = getSV(hash, p3);
                sv.type = "INV";
                sv.pairs += 0;
                sv.splits += cnt3;
                sv.clusters += 0;

                adjCnt(vref, sc3v, conf);
                incCnt(cov, p3, cnt3);

                if (cov.containsKey(bp) && cov.get(p3) < cov.get(bp)) {
                    cov.put(p3, cov.get(bp));
                }
                if (conf.y) {
                    System.err.printf("    Found INV: %s BP: %s Cov: %s %s %s EXTRA: %s Cnt: %s\n",
                            p3, bp, cov.get(p3), cov.get(bp), vn, EXTRA, cnt3);
                }
            }
        }
    }

    /**
     * Find INV SV with discordant pairs only
     * Would only consider those with supports from both orientations
     * (svfinv5) --&gt; | &lt;-- (svrinv5) ..... (svfinv3) --&gt; | &lt;-- (svrinv3)
     * @param hash map of variants on positions
     * @param cov map of coverage on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param reference reference in a given region
     * @param region region of interest
     * @param svfinv3 list of INV SVs in forward direction of 3'
     * @param svfinv5 list of INV SVs in forward direction of 5'
     * @param svrinv3 list of INV SVs in reverse direction of 3'
     * @param svrinv5 list of INV SVs in reverse direction of 5'
     * @param conf Configuration
     * @param referenceResource object for access to reference map
     * @param RLEN max read length
     * @param chrs map of chromosome lengths
     */
    static void findINVdisc (Map<Integer, VariationMap<String, Variation>> hash,
                             Map<Integer, Integer> cov,
                             Map<Integer, Sclip> sclip5,
                             Map<Integer, Sclip> sclip3,
                             Reference reference, ReferenceResource referenceResource, Region region,
                             List<Sclip> svfinv3,
                             List<Sclip> svrinv3,
                             List<Sclip> svfinv5,
                             List<Sclip> svrinv5,
                             Configuration conf,
                             int RLEN,
                             Map<String, Integer> chrs) {
        Map<Integer, Character> ref = reference.referenceSequences;
        for (Sclip invf5 : svfinv5) {
            if (invf5.used) {
                continue;
            }
            int cnt = invf5.cnt;
            int me = invf5.mend;
            int ms = invf5.mstart;
            int end = invf5.end;
            int start = invf5.start;
            double nm = invf5.nm;
            double pmean = invf5.pmean;
            double qmean = invf5.qmean;
            double Qmean = invf5.Qmean;
            if (!(Qmean/ (double) cnt > Configuration.DISCPAIRQUAL)) {
                continue;
            }
            for (Sclip invr5 : svrinv5) {
                if (invr5.used) {
                    continue;
                }
                int rcnt = invr5.cnt;
                int rstart = invr5.start;
                int rms = invr5.mstart;
                double rnm = invr5.nm;
                double rpmean = invr5.pmean;
                double rqmean = invr5.qmean;
                double rQmean = invr5.Qmean;
                if (!(rQmean/(double)rcnt > Configuration.DISCPAIRQUAL)) {
                    continue;
                }
                if (!(cnt + rcnt > conf.minr + 5)) {
                    continue;
                }
                if (isOverlap(end, me, new Tuple.Tuple2<>(rstart, rms), RLEN)) {
                    int bp = abs((end+rstart)/2);
                    int pe = abs((me+rms)/2);
                    if (!ref.containsKey(pe)) {
                        Region modifiedRegion = Region.newModifiedRegion(region, pe - 150, pe + 150);
                        referenceResource.getREF(modifiedRegion, chrs, conf, 300, reference);
                    }
                    int len = pe - bp + 1;

                    String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp, bp + Configuration.SVFLANK - 1));
                    String ins3 = SequenceUtil.reverseComplement(joinRef(ref, pe - Configuration.SVFLANK + 1 ,pe));
                    String ins = ins3 + "<inv" + (len - 2 * Configuration.SVFLANK) + ">" + ins5;
                    if (len - 2 * Configuration.SVFLANK <= 0) {
                        ins = SequenceUtil.reverseComplement(joinRef(ref, bp, pe));
                    }
                    if (conf.y) {
                        System.err.printf("  Found INV with discordant pairs only 5': cnt: %d Len: %d %d-%d<->%d-%d %s\n",
                                cnt, len, end, rstart, me, rms, ins);
                    }
                    final Variation vref = getVariation(hash, bp,"-" + len + "^" + ins);

                    // TODO: 26/07/18 Seem that next statement is redundant since it assigns cnt to 0 if it is already 0
                    //$hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } = 0 unless( $hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } );

                    invf5.used = true;
                    invr5.used = true;
                    vref.pstd = true;
                    vref.qstd = true;

                    Variation tmp = new Variation();
                    tmp.cnt = cnt + rcnt;
                    tmp.hicnt = cnt + rcnt;
                    tmp.dirPlus = cnt;
                    tmp.dirMinus = rcnt;
                    tmp.qmean = qmean + rqmean;
                    tmp.pmean = pmean + rpmean;
                    tmp.Qmean = Qmean + rQmean;
                    tmp.nm = nm + rnm;
                    adjCnt(vref, tmp, conf);

                    SV sv = getSV(hash, bp);
                    sv.type = "INV";
                    sv.pairs += cnt;
                    sv.splits += sclip5.containsKey(start) ? sclip5.get(start).cnt : 0;
                    sv.splits += sclip5.containsKey(ms) ? sclip5.get(ms).cnt : 0;
                    sv.clusters ++;

                    if (!cov.containsKey(bp)) {
                        cov.put(bp, 2 * cnt);
                    }
                    markSV(bp, pe, Arrays.asList(svfinv3, svrinv3), RLEN, conf);
                }
            }
        }

        for (Sclip invf3 : svfinv3) {
            if (invf3.used) {
                continue;
            }
            int cnt = invf3.cnt;
            int me = invf3.mend;
            int end = invf3.end;
            double nm = invf3.nm;
            double pmean = invf3.pmean;
            double qmean = invf3.qmean;
            double Qmean = invf3.Qmean;
            for (Sclip invr3 : svrinv3) {
                if (invr3.used) {
                    continue;
                }
                int rcnt = invr3.cnt;
                int rstart = invr3.start;
                int rms = invr3.mstart;
                double rnm = invr3.nm;
                double rpmean = invr3.pmean;
                double rqmean = invr3.qmean;
                double rQmean = invr3.Qmean;
                if (!(rQmean/ (double)rcnt > Configuration.DISCPAIRQUAL)) {
                    continue;
                }
                if (!(cnt + rcnt > conf.minr + 5)) {
                    continue;
                }
                if (isOverlap(me, end, new Tuple.Tuple2<>(rms, rstart), RLEN)) {
                    int pe = abs((end + rstart)/2);
                    int bp = abs((me + rms)/2);
                    if (!ref.containsKey(bp)) {
                        Region modifiedRegion = Region.newModifiedRegion(region, bp - 150, bp + 150);
                        referenceResource.getREF(modifiedRegion, chrs, conf, 300, reference);
                    }
                    int len = pe - bp + 1;
                    String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp, bp + Configuration.SVFLANK - 1));
                    String ins3 = SequenceUtil.reverseComplement(joinRef(ref, pe - Configuration.SVFLANK + 1, pe));

                    String ins = ins3+"<inv" + (len - 2* Configuration.SVFLANK) + ">"+ ins5;
                    if (len - 2 * Configuration.SVFLANK <= 0 ) {
                        ins = SequenceUtil.reverseComplement(joinRef(ref, bp, pe));
                    }
                    if (conf.y) {
                        System.err.printf("  Found INV with discordant pairs only 3': cnt: %d Len: %d %d-%d<->%d-%d %s\n",
                                cnt, len, me, rms, end, rstart, ins);
                    }
                    final Variation vref = getVariation(hash, bp,"-" + len + "^" + ins);
                    // TODO: 26/07/18 Seem that next statement is redundant since it assigns cnt to 0 if it is already 0
                    //$hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } = 0 unless( $hash->{ $bp }->{ "-${len}^$ins" }->{ cnt } );

                    invf3.used = true;
                    invr3.used = true;
                    vref.pstd = true;
                    vref.qstd = true;

                    Variation tmp = new Variation();
                    tmp.cnt = cnt + rcnt;
                    tmp.hicnt = cnt + rcnt;
                    tmp.dirPlus = cnt;
                    tmp.dirMinus = rcnt;
                    tmp.qmean = qmean + rqmean;
                    tmp.pmean = pmean + rpmean;
                    tmp.Qmean = Qmean + rQmean;
                    tmp.nm = nm + rnm;

                    adjCnt(vref, tmp, conf);

                    SV sv = getSV(hash, bp);
                    sv.type = "INV";
                    sv.pairs += cnt;
                    sv.splits += sclip3.containsKey(end + 1) ? sclip3.get(end + 1).cnt : 0;
                    sv.splits += sclip3.containsKey(me + 1) ? sclip3.get(me + 1).cnt : 0;
                    sv.clusters++;

                    if(!cov.containsKey(bp)) {
                        cov.put(bp, 2 * cnt);
                    }
                    markSV(bp, pe, Arrays.asList(svfinv5, svrinv5), RLEN, conf);
                }
            }
        }
    }

    /**
     * Find DUP SVs with discordant pairs only
     * (svrdup) |&lt;--.........--&gt;| (svfdup)
     * @param hash map of variants on positions
     * @param iHash map of insertions on positions
     * @param cov map of coverage on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param reference reference in a given region
     * @param region region of interest
     * @param bam BAM file name
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param svfdup list of DUP SVs in discordant pairs in forward direction
     * @param svrdup list of DUP SVs in discordant pairs in reverse direction
     * @param conf Configuration
     * @param RLEN max read length
     * @param referenceResource object for access to reference map
     * @param chrs map of chromosome lengths
     * @throws IOException if BAM file can't be read
     */
    static void findDUPdisc (Map<Integer, VariationMap<String, Variation>> hash,
                             Map<Integer, VariationMap<String, Variation>> iHash,
                             Map<Integer, Integer> cov,
                             Map<Integer, Sclip> sclip5,
                             Map<Integer, Sclip> sclip3,
                             Reference reference, ReferenceResource referenceResource,
                             Region region,
                             String bam,
                             String sample,
                             Set<String> splice,
                             String ampliconBasedCalling,
                             List<Sclip> svfdup,
                             List<Sclip> svrdup,
                             Configuration conf,
                             int RLEN,
                             Map<String, Integer> chrs) throws IOException {
        Map<Integer, Character> ref = reference.referenceSequences;
        for (Sclip dup : svfdup) {
            if (dup.used) {
                continue;
            }
            int ms = dup.mstart;
            int me = dup.mend;
            int cnt = dup.cnt;
            int end = dup.end;
            int start = dup.start;
            double pmean = dup.pmean;
            double qmean = dup.qmean;
            double Qmean = dup.Qmean;
            double nm = dup.nm;
            if (!(cnt >= conf.minr + 5)) {
                continue;
            }
            if (!(Qmean/cnt > Configuration.DISCPAIRQUAL)) {
                continue;
            }
            int mlen = end - ms + RLEN/cnt;
            int bp = ms - (RLEN/cnt)/2;
            int pe = end;

            if(!ReferenceResource.isLoaded(region.chr, ms, me, reference)) {
                if (!ref.containsKey(bp)) {
                    referenceResource.getREF(Region.newModifiedRegion(region, bp - 150, bp + 150), chrs, conf, 300, reference);
                }
                parseSAM(Region.newModifiedRegion(region, ms - 200, me + 200), bam, chrs, sample, splice,
                        ampliconBasedCalling, RLEN, reference, referenceResource, conf, hash,
                        iHash, cov, sclip3, sclip5, true);
            }

            int cntf = cnt;
            int cntr = cnt;
            double qmeanf = qmean;
            double qmeanr = qmean;
            double Qmeanf = Qmean;
            double Qmeanr = Qmean;
            double pmeanf = pmean;
            double pmeanr = pmean;
            double nmf = nm;
            double nmr = nm;

            if (!dup.soft.isEmpty()) {
                List<Map.Entry<Integer, Integer>> soft = new ArrayList<>(dup.soft.entrySet());
                soft.sort(comparing((Map.Entry<Integer, Integer> entry) -> entry.getValue(),
                        Integer::compareTo).reversed());

                if (soft.size() > 0) {
                    pe = soft.get(0).getKey();
                }
                if (!sclip3.containsKey(pe)) {
                    continue;
                }
                if (sclip3.get(pe).used) {
                    continue;
                }
                Sclip currentSclip3 =  sclip3.get(pe);
                cntf = currentSclip3.cnt;
                qmeanf = currentSclip3.qmean;
                Qmeanf = currentSclip3.Qmean;
                pmeanf = currentSclip3.pmean;
                nmf = currentSclip3.nm;

                String seq = findconseq(currentSclip3, conf, 0);
                Tuple.Tuple2<Integer, String> tuple = findMatch(conf, seq, reference, bp, 1);
                int tbp = tuple._1;
                String EXTRA = tuple._2;

                if (tbp != 0 && tbp < pe) {
                    currentSclip3.used = true;
                    while(isHasAndEquals(pe - 1, ref, tbp - 1)) {
                        tbp--;
                        if (tbp != 0) {
                            pe--;
                        }
                    }
                    mlen = pe - tbp;
                    bp = tbp;
                    pe--;
                    end = pe;
                    if (sclip5.containsKey(bp)) {
                        Sclip currentSclip5 = sclip5.get(bp);
                        cntr = currentSclip5.cnt;
                        qmeanr = currentSclip5.qmean;
                        Qmeanr = currentSclip5.Qmean;
                        pmeanr = currentSclip5.pmean;
                        nmr = currentSclip5.nm;
                    }
                }
            }

            String ins5 = joinRef(ref, bp, bp + Configuration.SVFLANK - 1);
            String ins3 = joinRef(ref, pe - Configuration.SVFLANK + 1, pe);
            String ins = ins5 + "<dup" + (mlen - 2 * Configuration.SVFLANK) + ">" + ins3;

            final Variation vref = getVariation(hash, bp, "+" + ins);
            vref.cnt = 0;

            SV sv = getSV(hash, bp);
            sv.type = "DUP";
            sv.pairs += cnt;
            sv.splits += dup.softp != 0 && sclip3.get(dup.softp) != null
                    ? sclip3.get(dup.softp).cnt : 0;
            sv.clusters++;

            if (conf.y) {
                System.err.printf("  Found DUP with discordant pairs only (forward): cnt: %d BP: %d END: %d %s Len: %d %d-%d<->%d-%d\n",
                        cnt, bp, pe, ins, mlen, start, end, ms, me);
            }
            int tcnt = cntr + cntf;

            Variation tmp = new Variation();
            tmp.cnt = tcnt;
            tmp.extracnt = tcnt;
            tmp.hicnt = tcnt;
            tmp.dirPlus = cntf;
            tmp.dirMinus = cntr;
            tmp.qmean = qmeanf + qmeanr;
            tmp.pmean = pmeanf + pmeanr;
            tmp.Qmean = Qmeanf + Qmeanr;
            tmp.nm = nmf + nmr;

            adjCnt(vref, tmp, conf);
            dup.used = true;
            if(!cov.containsKey(bp)) {
                cov.put(bp, tcnt);
            }
            if (cov.containsKey(end) && cov.get(bp) < cov.get(end)) {
                cov.put(bp, cov.get(end));
            }

            Tuple.Tuple2<Integer, Integer> tuple = markDUPSV(bp, pe, Collections.singletonList(svrdup), RLEN, conf);
            int clusters = tuple._1;
            sv.clusters += clusters;
        }

        for (Sclip dup : svrdup) {
            if (dup.used) {
                continue;
            }
            int ms = dup.mstart;
            int me = dup.mend;
            int cnt = dup.cnt;
            int end = dup.end;
            int start = dup.start;
            double pmean = dup.pmean;
            double qmean = dup.qmean;
            double Qmean = dup.Qmean;
            double nm = dup.nm;
            if (cnt < conf.minr + 5) {
                continue;
            }
            if (!(Qmean/cnt > Configuration.DISCPAIRQUAL)) {
                continue;
            }
            int mlen = me - start + RLEN/cnt;
            int bp = start - (RLEN/cnt)/2;
            int pe = mlen + bp - 1;
            int tpe = pe;
            if (!ReferenceResource.isLoaded(region.chr, ms, me, reference) ) {
                if (!ref.containsKey(pe)) {
                    referenceResource.getREF(Region.newModifiedRegion(region, pe - 150, pe + 150), chrs, conf, 300, reference);
                }
                parseSAM(Region.newModifiedRegion(region, ms - 200, me + 200), bam, chrs, sample, splice,
                        ampliconBasedCalling, RLEN, reference, referenceResource, conf, hash,
                        iHash, cov, sclip3, sclip5, true);
            }
            int cntf = cnt;
            int cntr = cnt;
            double qmeanf = qmean;
            double qmeanr = qmean;
            double Qmeanf = Qmean;
            double Qmeanr = Qmean;
            double pmeanf = pmean;
            double pmeanr = pmean;
            double nmf = nm;
            double nmr = nm;
            if (!dup.soft.isEmpty()) {
                List<Map.Entry<Integer, Integer>> soft = new ArrayList<>(dup.soft.entrySet());
                soft.sort(comparing((Map.Entry<Integer, Integer> entry) -> entry.getValue(),
                        Integer::compareTo).reversed());

                if (soft.size() > 0) {
                    bp = soft.get(0).getKey();
                }
                if (!sclip5.containsKey(bp)) {
                    continue;
                }
                Sclip currentSclip5 = sclip5.get(bp);
                if (currentSclip5.used) {
                    continue;
                }
                cntr = currentSclip5.cnt;
                qmeanr = currentSclip5.qmean;
                Qmeanr = currentSclip5.Qmean;
                pmeanr = currentSclip5.pmean;
                nmr = currentSclip5.nm;
                String seq = findconseq(currentSclip5, conf, 0);
                Tuple.Tuple2<Integer, String> tuple = findMatch(conf, seq, reference, pe, -1);
                int tbp = tuple._1;
                String EXTRA = tuple._2;
                if (tbp != 0 && tbp > bp) {
                    currentSclip5.used = true;
                    pe = tbp;
                    mlen = pe - bp + 1;
                    tpe = pe + 1;
                    while(isHasAndEquals(tpe, ref, bp + (tpe-pe-1))) {
                        tpe++;
                    }
                    if (sclip3.containsKey(tpe)) {
                        Sclip currentSclip3 = sclip3.get(tpe);
                        cntf = currentSclip3.cnt;
                        qmeanf = currentSclip3.qmean;
                        Qmeanf = currentSclip3.Qmean;
                        pmeanf = currentSclip3.pmean;
                        nmf = currentSclip3.nm;
                    }
                }
            }

            String ins5 = joinRef(ref, bp, bp + Configuration.SVFLANK - 1);
            String ins3 = joinRef(ref, pe - Configuration.SVFLANK + 1, pe);
            String ins = ins5 + "<dup" + (mlen - 2 * Configuration.SVFLANK) + ">" + ins3;

            final Variation vref = getVariation(hash, bp, "+" + ins);
            vref.cnt = 0;

            SV sv = getSV(hash, bp);
            sv.type = "DUP";
            sv.pairs += cnt;
            sv.splits += sclip5.containsKey(bp) ? sclip5.get(bp).cnt : 0;
            sv.splits += sclip3.containsKey(tpe) ? sclip3.get(tpe).cnt : 0;
            sv.clusters++;

            if (conf.y) {
                System.err.printf("  Found DUP with discordant pairs only (reverse): cnt: %d BP: %d Len: %d %d-%d<->%d-%d\n",
                        cnt, bp, mlen, start, end, ms, me);
            }
            int tcnt = cntr + cntf;
            Variation tmp = new Variation();
            tmp.cnt = tcnt;
            tmp.extracnt = tcnt;
            tmp.hicnt = tcnt;
            tmp.dirPlus = cntf;
            tmp.dirMinus = cntr;
            tmp.qmean = qmeanf + qmeanr;
            tmp.pmean = pmeanf + pmeanr;
            tmp.Qmean = Qmeanf + Qmeanr;
            tmp.nm = nmf + nmr;
            adjCnt(vref, tmp, conf);

            dup.used = true;
            if (!cov.containsKey(bp)) {
                cov.put(bp, tcnt);
            }
            if (cov.containsKey(me) && cov.get(bp) < cov.get(me)) {
                cov.put(bp, cov.get(me));
            }
            Tuple.Tuple2<Integer, Integer> tuple = markDUPSV(bp, pe, Collections.singletonList(svfdup), RLEN, conf);
            int clusters = tuple._1;
            sv.clusters += clusters;
        }
    }

    /**
     * Find matching on reverse strand (without MM defined, by default MM = 3)
     * @param conf Configuration
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    static Tuple.Tuple2<Integer, String> findMatchRev(Configuration conf,
                                                      String seq,
                                                      Reference REF,
                                                      int position,
                                                      int dir) {
        int MM = 3;
        return findMatchRev(conf, seq, REF, position, dir, Configuration.SEED_1, MM);
    }

    /**
     * Find matching on reverse strand (with MM defined)
     * @param conf Configuration
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @param SEED seed length from configuration
     * @param MM the minimum matches for a read to be considered
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    static Tuple.Tuple2<Integer, String> findMatchRev(Configuration conf,
                                                      String seq,
                                                      Reference REF,
                                                      int position,
                                                      int dir,
                                                      int SEED,
                                                      int MM) {
        // dir = 1 means from 3' soft clip
        if (dir == 1) {
            seq = reverse(seq);
        }
        seq = complement(seq);
        if (conf.y) {
            System.err.printf("    Working MatchRev %d %s %d%n", position, seq, dir);
        }
        String extra = "";
        for (int i = seq.length() - SEED; i >= 0; i--) {
            String seed = substr(seq, i, SEED);
            if (REF.seed.containsKey(seed)) {
                List<Integer> seeds = REF.seed.get(seed);
                if (seeds.size() == 1) {
                    Integer firstSeed = seeds.get(0);
                    int bp = dir == 1 ? firstSeed + seq.length() - i - 1 : firstSeed - i;
                    if (ismatchref(seq, REF.referenceSequences, bp, -1 * dir, conf.y, MM)) {
                        if (conf.y) {
                            System.err.printf(
                                    "      Found SV BP (reverse): %d BP: %d SEEDpos: %d %d %s %d %s \n",
                                    dir, bp, firstSeed, position, seed, i, seq
                            );
                        }
                        return Tuple.tuple(bp, extra);
                    } else {
                        // for complex indels, allowing some mismatches at the end up to 15bp or 20% length
                        String sseq = seq;
                        int eqcnt = 0;
                        for (int j = 1; j <= 15; j++) {
                            bp -= dir;
                            sseq = dir == -1 ? substr(sseq, 1) : substr(sseq, 0, -1);
                            if (dir == -1) {
                                if (isHasAndNotEquals(charAt(sseq, 0), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, 1), REF.referenceSequences, bp+1)) {
                                    continue;
                                }
                                extra = substr(seq, 0, j);
                            } else {
                                if (isHasAndNotEquals(charAt(sseq, -1), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, -2), REF.referenceSequences, bp - 1 )) {
                                    continue;
                                }
                                extra = substr(seq, -j);
                            }
                            if (eqcnt >= 3 && eqcnt / (double)j > 0.5) {
                                break;
                            }
                            if (conf.y) {
                                System.err.printf(
                                        "      FoundSEED SV BP (reverse): %d BP: %d SEEDpos%d %d %s %d %s EXTRA: %s%n",
                                        dir, bp, firstSeed, position, seed, i, seq, extra);
                            }
                            if (ismatchref(sseq, REF.referenceSequences, bp, -1 * dir, conf.y,1)) {
                                return Tuple.tuple(bp, extra);
                            }
                        }
                    }
                }
            }
        }
        return Tuple.tuple(0, "");
    }

    /**
     * Find matching on forward strand (without MM defined, by default MM = 3)
     * @param conf Configuration
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    static Tuple.Tuple2<Integer, String> findMatch(Configuration conf,
                                                   String seq,
                                                   Reference REF,
                                                   int position,
                                                   int dir) {
        int MM = 3;
        return findMatch(conf, seq, REF, position, dir, Configuration.SEED_1, MM);
    }

    /**
     * Find matching on forward strand (with MM defined)
     * @param conf Configuration
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @param SEED seed length from configuration
     * @param MM the minimum matches for a read to be considered
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    static Tuple.Tuple2<Integer, String> findMatch(Configuration conf,
                                                   String seq,
                                                   Reference REF,
                                                   int position,
                                                   int dir,
                                                   int SEED,
                                                   int MM) {
        if (dir == -1) {
            seq = reverse(seq); // dir==-1 means 5' clip
        }
        if (conf.y) {
            System.err.printf("    Working Match %d %s %d SEED: %d\n", position, seq, dir, SEED);
        }
        String extra = "";
        for(int i = seq.length() - SEED; i >= 0; i--) {
            String seed = substr(seq, i, SEED);
            if (REF.seed.containsKey(seed)) {
                List<Integer> seeds = REF.seed.get(seed);
                if (seeds.size() == 1 ) {
                    Integer firstSeed = seeds.get(0);
                    int bp = dir == 1 ? firstSeed - i : firstSeed + seq.length() - i - 1;

                    if (ismatchref(seq, REF.referenceSequences, bp, dir, conf.y, MM)) {
                        int mm = dir == -1 ? -1 : 0;
                        while(isHasAndNotEquals(charAt(seq, mm), REF.referenceSequences, bp)) {
                            extra += substr(seq, mm, 1);
                            bp += dir;
                            mm += dir;
                        }
                        if (!extra.isEmpty() && dir == -1 ) {
                            extra = reverse(extra);
                        }
                        if (conf.y) {
                            System.err.printf("      Found SV BP: %d BP: %d SEEDpos %s %d %s %d %s extra: %s\n",
                                    dir, bp, firstSeed, position, seed, i, seq, extra);
                        }
                        return Tuple.tuple(bp, extra);
                    } else { // for complex indels, allowing some mismatches at the end up to 15bp or 20% length
                        String sseq = seq;
                        int eqcnt = 0;
                        for(int ii = 1; ii <= 15; ii++) {
                            bp += dir;
                            sseq = dir == 1 ? substr(sseq, 1) : substr(sseq, 0, -1);
                            if(dir == 1) {
                                if (isHasAndNotEquals(charAt(sseq, 0), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, 1), REF.referenceSequences, bp + 1)) {
                                    continue;
                                }
                                extra = substr(seq, 0, ii);
                            } else {
                                if (isHasAndNotEquals(charAt(sseq, -1), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, -2), REF.referenceSequences, bp - 1)) {
                                    continue;
                                }
                                extra = substr(seq, -ii);
                            }
                            if  (eqcnt >= 3 && eqcnt/(double)ii > 0.5) {
                                break;
                            }
                            if (conf.y) {
                                System.err.printf("      FoundSEED SV BP: %d BP: %d SEEDpos%s %d %s %d %s EXTRA: %s\n",
                                        dir, bp, firstSeed, position, seed, i, seq, extra);
                            }
                            if (ismatchref(sseq, REF.referenceSequences, bp, dir, conf.y, 1) ) {
                                return Tuple.tuple(bp, extra);
                            }
                        }
                    }
                }
            }
        }
        return Tuple.tuple(0, "");
    }

    /**
     * Find DEL SV with discordant pairs only
     * (svfdel) --&gt; | ............ | &lt;-- (svrdel)
     * @param hash map of variants on positions
     * @param cov map of coverage on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param REF reference in a given region
     * @param region region of interest
     * @param splice set of strings representing spliced regions
     * @param svfdel list of DEL SVs in discordant pairs in forward direction
     * @param svrdel list of DEL SVs in discordant pairs in reverse direction
     * @param conf Configuration
     * @param referenceResource object for access to reference map
     * @param rlen max read length
     * @param chrs map of chromosome lengths
     */
    static void findDELdisc(Configuration conf,
                            Map<Integer, VariationMap<String, Variation>> hash,
                            Map<Integer, Integer> cov,
                            Map<Integer, Sclip> sclip5,
                            Map<Integer, Sclip> sclip3,
                            Reference REF, ReferenceResource referenceResource,
                            Region region,
                            List<Sclip> svfdel,
                            List<Sclip> svrdel,
                            Set<String> splice,
                            int rlen,
                            Map<String, Integer> chrs) {
        int MINDIST = 8 * rlen; // the minimum distance between two clusters
        for (Sclip del : svfdel) {
            if (del.used) {
                continue;
            }
            if (!splice.isEmpty() && abs(del.mlen) < 250000) {
                continue; // more stringent for RNA-Seq
            }

            if (del.cnt < conf.minr + 5) {
                continue;
            }
            if (del.mstart <= del.end + MINDIST) {
                continue;
            }
            if (del.Qmean / del.cnt <= Configuration.DISCPAIRQUAL) {
                continue;
            }
            int mlen = del.mstart - del.end - rlen / (del.cnt + 1);
            if (!(mlen > 0 && mlen > MINDIST)) {
                continue;
            }
            int bp = del.end + (rlen / (del.cnt + 1)) / 2;
            if (del.softp != 0) {
                bp = del.softp;
            }
            final Variation vref = getVariation(hash, bp, "-" + mlen);
            vref.cnt = 0;
            SV sv = getSV(hash, bp);
            sv.type = "DEL";
            sv.splits += sclip3.containsKey(del.end + 1) ? sclip3.get(del.end + 1).cnt : 0;
            sv.splits += sclip5.containsKey(del.mstart) ? sclip5.get(del.mstart).cnt : 0;
            sv.pairs += del.cnt;
            sv.clusters++;

            if (conf.y) {
                System.err.printf(
                        "  Found DEL with discordant pairs only: cnt: %d BP: %d Len: %d %d-%d<->%d-%d%n",
                        del.cnt, bp, mlen, del.start, del.end, del.mstart, del.mend
                );
            }
            Variation tv = new Variation();
            tv.cnt = 2 * del.cnt;
            tv.hicnt = 2 * del.cnt;
            tv.dirPlus = del.cnt;
            tv.dirMinus = del.cnt;
            tv.qmean = 2 * del.qmean;
            tv.pmean = 2 * del.pmean;
            tv.Qmean = 2 * del.Qmean;
            tv.nm = 2 * del.nm;
            adjCnt(vref, tv, conf);
            if (!cov.containsKey(bp)) {
                cov.put(bp, 2 * del.cnt);
            }
            del.used = true;
            markSV(del.end, del.mstart, Arrays.asList(svrdel), rlen, conf);
        }
        for (Sclip del : svrdel) {
            if (del.used) {
                continue;
            }
            if (!splice.isEmpty() && abs(del.mlen) < 250000) {
                continue; // more stringent for RNA-Seq
            }
            if (del.cnt < conf.minr + 5) {
                continue;
            }
            if (del.start <= del.mend + MINDIST) {
                continue;
            }
            if (del.Qmean / del.cnt <= Configuration.DISCPAIRQUAL) {
                continue;
            }
            int mlen = del.start - del.mend - rlen / (del.cnt + 1);
            if (!(mlen > 0 && mlen > MINDIST)) {
                continue;
            }
            int bp = del.mend + ((rlen / (del.cnt + 1)) / 2);

            final Variation ref = getVariation(hash, bp, "-" + mlen);
            ref.cnt = 0;

            SV sv = getSV(hash, bp);
            sv.type = "DEL";
            sv.splits += sclip3.containsKey(del.mend + 1) ? sclip3.get(del.mend + 1).cnt : 0;
            sv.splits += sclip5.containsKey(del.start) ? sclip5.get(del.start).cnt : 0;
            sv.pairs += del.cnt;
            sv.clusters += 1;

            if (conf.y) {
                System.err.printf(
                        "  Found DEL with discordant pairs only (reverse): cnt: %d BP: %d Len: %d %d-%d<->%d-%d%n",
                        del.cnt, bp, mlen, del.start, del.end, del.mstart, del.mend
                );
            }
            if (del.softp != 0 && sclip5.containsKey(del.softp)) {
                sclip5.get(del.softp).used = true;
            }
            Variation tv = new Variation();
            tv.cnt = 2 * del.cnt;
            tv.hicnt = 2 * del.cnt;
            tv.dirPlus = del.cnt;
            tv.dirMinus = del.cnt;
            tv.qmean = 2 * del.qmean;
            tv.pmean = 2 * del.pmean;
            tv.Qmean = 2 * del.Qmean;
            tv.nm = 2 * del.nm;
            adjCnt(ref, tv, conf);
            if (!cov.containsKey(bp)) {
                cov.put(bp, 2 * del.cnt);
            }
            if (cov.containsKey(del.start) && cov.get(bp) < cov.get(del.start)) {
                cov.put(bp, cov.get(del.start));
            }
            del.used = true;
            referenceResource.getREF(Region.newModifiedRegion(region, del.mstart - 100, del.mend + 100), chrs, conf, 200, REF);
            markSV(del.mend, del.start, Arrays.asList(svfdel), rlen, conf);
        }
    }

    /**
     * Find DEL SV
     * @param conf Configuration
     * @param hash map of variants on positions
     * @param iHash map of insertions on positions
     * @param cov map of coverage on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param reference reference in a given region
     * @param region region of interest
     * @param bams list of BAM files from configuration
     * @param svfdel list of DEL SVs in forward direction
     * @param svrdel list of DEL SVs in reverse direction
     * @param rlen max read length
     * @param chrs map of chromosome lengths
     * @param sample sample name
     * @param referenceResource object for access to reference map
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @throws IOException if BAM file can't be read
     */
    static void findDEL(Configuration conf,
                        Map<Integer, VariationMap<String, Variation>> hash,
                        Map<Integer, VariationMap<String, Variation>> iHash,
                        Map<Integer, Integer> cov,
                        Map<Integer, Sclip> sclip5,
                        Map<Integer, Sclip> sclip3,
                        Reference reference, ReferenceResource referenceResource,
                        Region region,
                        String bams,
                        List<Sclip> svfdel,
                        List<Sclip> svrdel,
                        int rlen,
                        Map<String, Integer> chrs,
                        String sample,
                        Set<String> splice,
                        String ampliconBasedCalling) throws IOException {
        for (Sclip del : svfdel) {
            if (del.used) {
                continue;
            }
            if (del.cnt < conf.minr) {
                continue;
            }
            List<Tuple.Tuple2<Integer, Integer>> soft = del.soft.entrySet().stream()
                    .map(entry -> Tuple.tuple(entry.getKey(), entry.getValue()))
                    .sorted(comparing(tuple -> tuple._2, reverseOrder()))
                    .collect(toList());

            int softp = soft.isEmpty() ? 0 : soft.get(0)._1;
            if (conf.y) System.err.printf("%n%nWorking DEL 5' %d mate cluster cnt: %d%n", softp, del.cnt);
            if (softp != 0) {
                if (!sclip3.containsKey(softp)) {
                    continue;
                }
                Sclip scv = sclip3.get(softp);
                if (scv.used) {
                    continue;
                }
                String seq = findconseq(scv, conf, 0);
                if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                    continue; // next unless( $seq && length($seq) >= $SEED2 );
                }
                if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                    referenceResource.getREF(Region.newModifiedRegion(region, del.mstart, del.mend), chrs, conf, 300, reference);
                    parseSAM(Region.newModifiedRegion(region, del.mstart - 200, del.mend + 200),
                            bams, chrs, sample, splice, ampliconBasedCalling, rlen, reference, referenceResource, conf, hash,
                            iHash, cov, sclip3, sclip5, true);
                }
                Tuple.Tuple2<Integer, String> match = findMatch(conf, seq, reference, softp, 1);
                int bp = match._1;
                String extra = match._2;
                if (bp == 0) {
                    continue;
                }
                if (!(bp - softp > 30 && isOverlap(softp, bp, Tuple.tuple(del.end, del.mstart), rlen))) {
                    continue;
                }
                bp--;
                int dellen = bp - softp + 1;
                Map<Integer, Character> ref = reference.referenceSequences;
                while(ref.containsKey(bp) && ref.containsKey(softp - 1) && ref.get(bp).equals(ref.get(softp - 1))) {
                    bp --;
                    if (bp != 0) {
                        softp--;
                    }
                }

                final Variation variation = getVariation(hash, softp, "-" + dellen);
                variation.cnt = 0;

                SV sv = getSV(hash, softp);
                sv.type = "DEL";
                sv.pairs += del.cnt;
                sv.splits += scv.cnt;
                sv.clusters++;

                if (!(cov.containsKey(softp) && cov.get(softp) > del.cnt)) {
                    cov.put(softp, del.cnt);
                }
                if (cov.containsKey(bp) && cov.get(softp) < cov.get(bp)) {
                    cov.put(softp, cov.get(bp));
                }

                adjCnt(variation, scv, hash.containsKey(softp) && ref.containsKey(softp)
                        ? hash.get(softp).get(ref.get(softp).toString()) : null, conf);

                int mcnt = del.cnt;
                Variation tv = new Variation();
                tv.cnt = mcnt;
                tv.hicnt = mcnt;
                tv.dirPlus = mcnt / 2;
                tv.dirMinus = mcnt - mcnt / 2;
                tv.qmean = del.qmean * mcnt / del.cnt;
                tv.pmean = del.pmean * mcnt / del.cnt;
                tv.Qmean = del.Qmean * mcnt / del.cnt;
                tv.nm = del.nm * mcnt / del.cnt;
                adjCnt(variation, tv, conf);

                del.used = true;
                markSV(softp, bp, Arrays.asList(svrdel), rlen, conf);
                if (conf.y) {
                    System.err.printf("    Found DEL SV from 5' softclip unhappy reads: %d -%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.cnt, mcnt);
                }
            } else { // Look within a read length
                if(!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                    referenceResource.getREF(Region.newModifiedRegion(region, del.mstart, del.mend), chrs, conf, 300, reference);
                }
                if (conf.y) {
                    System.err.printf("%n%nWorking DEL 5' no softp mate cluster cnt: %d%n", del.cnt);
                }

                for (Map.Entry<Integer, Sclip> entry : sclip3.entrySet()) {
                    Integer i = entry.getKey();
                    Sclip scv = entry.getValue();

                    if (scv.used) {
                        continue;
                    }
                    if (!(i >= del.end -3 && i - del.end < 3 * rlen)) {
                        continue;
                    }
                    String seq = findconseq(scv, conf, 0);
                    if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                        continue;
                    }
                    softp = i;
                    Tuple.Tuple2<Integer, String> match = findMatch(conf, seq, reference, softp, 1);
                    int bp = match._1;
                    String EXTRA = match._2;
                    if (bp == 0) {
                        match = findMatch(conf, seq, reference, softp, 1, Configuration.SEED_2, 0);
                        bp = match._1;
                        EXTRA = match._2;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    if (!(bp - softp > 30 && isOverlap(softp, bp, Tuple.tuple(del.end, del.mstart), rlen))) continue;
                    bp--;
                    int dellen = bp - softp + 1;

                    final Variation variation = getVariation(hash, softp, "-" + dellen);
                    variation.cnt = 0;

                    SV sv = getSV(hash, softp);
                    sv.type = "DEL";
                    sv.pairs += del.cnt;
                    sv.splits += scv.cnt;
                    sv.clusters++;

                    if (!(cov.containsKey(softp) && cov.get(softp) > del.cnt)) {
                        cov.put(softp, del.cnt);
                    }
                    if (cov.containsKey(bp) && cov.get(softp) < cov.get(bp)) {
                        cov.put(softp, cov.get(bp));
                    }
                    adjCnt(variation, scv, conf);
                    int mcnt = del.cnt;
                    Variation tv = new Variation();
                    tv.cnt = mcnt;
                    tv.hicnt = mcnt;
                    tv.dirPlus = mcnt / 2;
                    tv.dirMinus = mcnt - mcnt / 2;
                    tv.qmean = del.qmean * mcnt / del.cnt;
                    tv.pmean = del.pmean * mcnt / del.cnt;
                    tv.Qmean = del.Qmean * mcnt / del.cnt;
                    tv.nm = del.nm * mcnt / del.cnt;
                    adjCnt(variation, tv, conf);
                    del.used = true;
                    markSV(softp, bp, Arrays.asList(svrdel), rlen, conf);
                    if (conf.y) {
                        System.err.printf("    Found DEL SV from 5' softclip happy reads: %d -%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.cnt, mcnt);
                    }
                    break;
                }
            }
        }
        for (Sclip del : svrdel) {
            if (del.used) {
                continue;
            }
            if (del.cnt < conf.minr) {
                continue;
            }
            List<Tuple.Tuple2<Integer, Integer>> soft = del.soft.entrySet().stream()
                    .map(entry -> Tuple.tuple(entry.getKey(), entry.getValue()))
                    .sorted(comparing(tuple -> tuple._2, reverseOrder()))
                    .collect(toList());

            int softp = soft.isEmpty() ? 0 : soft.get(0)._1;
            if (softp != 0) {
                if (conf.y) {
                    System.err.printf("%n%nWorking DEL 3' %d mate cluster cnt: %s%n", softp, del.cnt);
                }
                if (!sclip5.containsKey(softp)) {
                    continue;
                }
                Sclip scv = sclip5.get(softp);
                if (scv.used) {
                    continue;
                }
                String seq = findconseq(scv, conf, 0);
                if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                    continue;
                }
                if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                    referenceResource.getREF(Region.newModifiedRegion(region, del.mstart, del.mend), chrs, conf, 300, reference);
                    parseSAM(Region.newModifiedRegion(region, del.mstart - 200, del.mend + 200),
                            bams, chrs, sample, splice, ampliconBasedCalling, rlen, reference, referenceResource, conf,
                            hash, iHash, cov, sclip3, sclip5, true);
                }
                Tuple.Tuple2<Integer, String> match = findMatch(conf, seq, reference, softp, -1);// my ($bp, $EXTRA) = findMatch($seq, $reference, $softp, -1);
                int bp = match._1;
                String EXTRA = match._2;
                if (bp == 0) {
                    match = findMatch(conf, seq, reference, softp, -1, Configuration.SEED_2, 0);
                    bp = match._1;
                    EXTRA = match._2;
                }
                if (bp == 0) {
                    continue;
                }
                if (!(softp - bp > 30 && isOverlap(bp, softp, Tuple.tuple(del.mend, del.start), rlen))) continue;
                bp++;
                softp--;
                int dellen = softp - bp + 1;
                final Variation variation = getVariation(hash, bp, "-" + dellen);
                variation.cnt = 0;

                SV sv = getSV(hash, bp);
                sv.type = "DEL";
                sv.pairs += del.cnt;
                sv.splits += scv.cnt;
                sv.clusters++;

                adjCnt(variation, scv, conf);
                if (!(cov.containsKey(bp) && cov.get(bp) > del.cnt)) {
                    cov.put(bp, del.cnt);
                }
                if (cov.containsKey(softp) && cov.get(softp) > cov.get(bp)) {
                    cov.put(bp, cov.get(softp));
                }
                int mcnt = del.cnt;
                Variation tv = new Variation();
                tv.cnt = mcnt;
                tv.hicnt = mcnt;
                tv.dirPlus = mcnt / 2;
                tv.dirMinus = mcnt - mcnt / 2;
                tv.qmean = del.qmean * mcnt / del.cnt;
                tv.pmean = del.pmean * mcnt / del.cnt;
                tv.Qmean = del.Qmean * mcnt / del.cnt;
                tv.nm = del.nm * mcnt / del.cnt;
                adjCnt(variation, tv, conf);
                del.used = true;
                markSV(bp, softp, Arrays.asList(svfdel), rlen, conf);
                if (conf.y) {
                    System.err.printf("    Found DEL SV from 3' softclip unhappy reads: %d -+%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.cnt, mcnt);
                }
            } else {
                if (conf.y) {
                    System.err.printf("%n%nWorking DEL 3' no softp mate cluster %s %d %d cnt: %d%n", region.chr, del.mstart, del.mend, del.cnt);
                }
                if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                    referenceResource.getREF(Region.newModifiedRegion(region, del.mstart, del.mend), chrs, conf, 300, reference);
                }
                for (Map.Entry<Integer, Sclip> entry : sclip5.entrySet()) {
                    int i = entry.getKey();
                    Sclip scv = entry.getValue();
                    if (scv.used) {
                        continue;
                    }
                    if (!(i <= del.start + 3 && del.start - i < 3 * rlen)) {
                        continue;
                    }
                    String seq = findconseq(scv, conf, 0);
                    if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                        continue;
                    }
                    softp = i;
                    Tuple.Tuple2<Integer, String> match = findMatch(conf, seq, reference, softp, -1);
                    int bp = match._1;
                    String EXTRA = match._2;
                    if (bp == 0) {
                        match = findMatch(conf, seq, reference, softp, - 1, Configuration.SEED_2, 0);
                        bp = match._1;
                        EXTRA = match._2;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    if (!(softp - bp > 30 && isOverlap(bp, softp, Tuple.tuple(del.mend, del.start), rlen))) {
                        continue;
                    }
                    bp++;
                    softp--;
                    int dellen = softp - bp + 1;

                    final Variation variation = getVariation(hash, bp, "-" + dellen);
                    variation.cnt = 0;

                    SV sv = getSV(hash, bp);
                    sv.type = "DEL";
                    sv.pairs += del.cnt;
                    sv.splits += scv.cnt;
                    sv.clusters++;

                    adjCnt(variation, scv, conf);
                    if (!cov.containsKey(bp)) {
                        cov.put(bp, del.cnt);
                    }
                    if (cov.containsKey(softp) && cov.get(softp) > cov.get(bp)) {
                        cov.put(bp, cov.get(softp));
                    }
                    incCnt(cov, bp, scv.cnt);

                    int mcnt = del.cnt;
                    Variation tv = new Variation();
                    tv.cnt = mcnt;
                    tv.hicnt = mcnt;
                    tv.dirPlus = mcnt / 2;
                    tv.dirMinus = mcnt - mcnt / 2;
                    tv.qmean = del.qmean * mcnt / del.cnt;
                    tv.pmean = del.pmean * mcnt / del.cnt;
                    tv.Qmean = del.Qmean * mcnt / del.cnt;
                    tv.nm = del.nm * mcnt / del.cnt;
                    adjCnt(variation, tv, conf);
                    del.used = true;
                    markSV(bp, softp, Arrays.asList(svfdel), rlen, conf);
                    if (conf.y) {
                        System.err.printf("    Found DEL SV from 3' softclip happy reads: %d -%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.cnt, mcnt);
                    }
                    break;
                }
            }
        }
    }

    /**
     * Find INV SV for all structural variants structures
     * @param conf Configuration
     * @param hash map of variants on positions
     * @param iHash map of insertions on positions
     * @param cov map of coverage on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param REF reference in a given region
     * @param region region of interest
     * @param bams list of BAM files from configuration
     * @param svfinv3 list of INV SVs in forward direction in 3'
     * @param svrinv3 list of INV SVs in reverse direction in 3'
     * @param svfinv5 list of INV SVs in forward direction in 5'
     * @param svrinv5 list of INV SVs in reverse direction in 5'
     * @param rlen max read length
     * @param chrs map of chromosome lengths
     * @param referenceResource object for access to reference map
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @throws IOException if BAM file can't be read
     */
    static void findINV(Configuration conf,
                        Map<Integer, VariationMap<String, Variation>> hash,
                        Map<Integer, VariationMap<String, Variation>> iHash,
                        Map<Integer, Integer> cov,
                        Map<Integer, Sclip> sclip5,
                        Map<Integer, Sclip> sclip3,
                        Reference REF, ReferenceResource referenceResource,
                        Region region,
                        String bams,
                        Iterable<Sclip> svfinv3,
                        Iterable<Sclip> svrinv3,
                        Iterable<Sclip> svfinv5,
                        Iterable<Sclip> svrinv5,
                        int rlen,
                        Map<String, Integer> chrs,
                        String sample,
                        Set<String> splice,
                        String ampliconBasedCalling) throws IOException {
        findINVsub(conf, hash, iHash, cov, sclip5, sclip3, REF, referenceResource, region, bams, svfinv5, 1 , Side._5, rlen, chrs, sample, splice, ampliconBasedCalling);
        findINVsub(conf, hash, iHash, cov, sclip5, sclip3, REF, referenceResource, region, bams, svrinv5, -1 , Side._5, rlen, chrs, sample, splice, ampliconBasedCalling);
        findINVsub(conf, hash, iHash, cov, sclip5, sclip3, REF, referenceResource, region, bams, svfinv3, 1 , Side._3, rlen, chrs, sample, splice, ampliconBasedCalling);
        findINVsub(conf, hash, iHash, cov, sclip5, sclip3, REF, referenceResource, region, bams, svrinv3, -1 , Side._3, rlen, chrs, sample, splice, ampliconBasedCalling);
    }

    /**
     * Find INV SV
     * Would only consider those with supports from both orientations
     * (svfinv5) --&gt; | &lt;-- (svrinv5) ..... (svfinv3) --&gt; | &lt;-- (svrinv3)
     * @param conf Configuration
     * @param hash map of variants on positions
     * @param iHash map of insertions on positions
     * @param cov map of coverage on positions
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param reference reference in a given region
     * @param region chromosome in region of interest
     * @param bams list of BAM files from configuration
     * @param svref list of INV SVs
     * @param dir direction specified in findINV()
     * @param side 3 or 5 end
     * @param rlen max read length
     * @param chrs map of chromosome lengths
     * @param referenceResource object for access to reference map
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @throws IOException if BAM file can't be read
     * @return created variation
     */
    static Variation findINVsub(Configuration conf,
                                Map<Integer, VariationMap<String, Variation>> hash,
                                Map<Integer, VariationMap<String, Variation>> iHash,
                                Map<Integer, Integer> cov,
                                Map<Integer, Sclip> sclip5,
                                Map<Integer, Sclip> sclip3,
                                Reference reference, ReferenceResource referenceResource,
                                Region region,
                                String bams,
                                Iterable<Sclip> svref,
                                int dir,
                                Side side,
                                int rlen,
                                Map<String, Integer> chrs,
                                String sample,
                                Set<String> splice,
                                String ampliconBasedCalling) throws IOException {
        // dir = 1 means 3' soft clipping
        for (Sclip inv : svref) {
            if (inv.used) {
                continue;
            }
            if (inv.cnt < conf.minr) {
                continue;
            }
            List<Tuple.Tuple2<Integer, Integer>> soft = inv.soft.entrySet().stream()
                    .map(entry -> Tuple.tuple(entry.getKey(), entry.getValue()))
                    .sorted(comparing(tuple -> tuple._2, reverseOrder()))
                    .collect(toList());
            int softp = soft.isEmpty() ? 0 : soft.get(0)._1;
            Map<Integer, Sclip> sclip = dir == 1 ? sclip3 : sclip5;
            if (conf.y) System.err.printf("%n%nWorking INV %d %d %s pair_cnt: %d%n", softp, dir, side, inv.cnt);
            if (!isLoaded(region.chr, inv.mstart, inv.mend, reference)) {
                referenceResource.getREF(Region.newModifiedRegion(region, inv.mstart, inv.mend), chrs, conf, 500, reference);
                parseSAM(Region.newModifiedRegion(region, inv.mstart - 200, inv.mend + 200),
                        bams, chrs, sample, splice, ampliconBasedCalling, rlen, reference, referenceResource, conf, hash,
                        iHash, cov, sclip3, sclip5, true);
            }
            int bp = 0;
            Sclip scv = new Sclip();
            String seq = "";
            String extra = "";
            if (softp != 0) {
                if (!sclip.containsKey(softp)) {
                    continue;
                }
                scv = sclip.get(softp);
                if (scv.used) {
                    continue;
                }
                seq = findconseq(scv, conf, 0);
                if (seq.isEmpty()) {
                    continue;
                }
                Tuple.Tuple2<Integer, String> matchRev = findMatchRev(conf, seq, reference, softp, dir);
                bp = matchRev._1;
                extra = matchRev._2;
                if (bp == 0 ) {
                    matchRev = findMatchRev(conf, seq, reference, softp, dir, Configuration.SEED_2, 0);
                    bp = matchRev._1;
                    extra = matchRev._2;
                }
                if (bp == 0) {
                    continue;
                }
            } else {
                // Look within 100bp to see whether a soft cliping can be found but not associated with discordant pairs
                int sp = dir == 1 ? inv.end : inv.start; // starting position
                for (int i = 1; i <= 2 * rlen; i++) {
                    int cp = sp + i * dir;
                    if (!sclip.containsKey(cp)) {
                        continue;
                    }
                    scv = sclip.get(cp);
                    if (scv.used) {
                        continue;
                    }
                    seq = findconseq(scv, conf, 0);
                    if (seq.isEmpty()) {
                        continue;
                    }
                    Tuple.Tuple2<Integer, String> matchRev = findMatchRev(conf, seq, reference, cp, dir);
                    bp = matchRev._1;
                    extra = matchRev._2;
                    if (bp == 0 ) {
                        matchRev = findMatchRev(conf, seq, reference, cp, dir, Configuration.SEED_2, 0);
                        bp = matchRev._1;
                        extra = matchRev._2;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    softp = cp;
                    if ((dir == 1 && abs(bp - inv.mend) < Configuration.MINSVCDIST * rlen)
                            || (dir == -1 && abs(bp - inv.mstart) < Configuration.MINSVCDIST * rlen)) break;
                }
                if (bp == 0) {
                    continue;
                }
            }
            if (conf.y) {
                System.err.printf("    %d %d %d %s %s pair_cnt: %d soft_cnt: %d%n", softp, bp, dir, side, seq, inv.cnt, scv.cnt);
            }
            if (side == Side._5) {
                if (dir == -1) {
                    bp--;
                }
            } else {
                if (dir == 1) {
                    bp++;
                    if (bp != 0) {
                        softp--;
                    }
                } else {
                    softp--;
                }
            }
            if (side == Side._3) {
                int tmp = bp;
                bp = softp;
                softp = tmp;
            }
            Map<Integer, Character> ref = reference.referenceSequences;
            if ((dir == -1 && side == Side._5) || dir == 1 && side == Side._3) {
                while (ref.containsKey(softp) && ref.containsKey(bp)
                        && ref.get(softp) == complement(ref.get(bp))) {
                    softp++;
                    if (softp != 0) {
                        bp--;
                    }
                }
            }
            while(ref.containsKey(softp - 1) && ref.containsKey(bp + 1)
                    && ref.get(softp - 1) == complement(ref.get(bp + 1))) {
                softp--;
                if (softp != 0) {
                    bp++;
                }
            }
            if (bp > softp && bp - softp > 150 && (bp - softp) / (double) abs(inv.mlen) < 1.5) {
                int len = bp - softp + 1;
                String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp - Configuration.SVFLANK + 1, bp));
                String ins3 = SequenceUtil.reverseComplement(joinRef(ref, softp, softp + Configuration.SVFLANK - 1));
                String ins = ins5 + "<inv" + (len - 2 * Configuration.SVFLANK) + ">" + ins3;
                if (len - 2 * Configuration.SVFLANK <= 0) {
                    ins = SequenceUtil.reverseComplement(joinRef(ref, softp, bp));
                }
                if (dir == 1 && !extra.isEmpty()) {
                    extra = SequenceUtil.reverseComplement(extra);
                    ins = extra + ins;
                } else if (dir == -1 && !extra.isEmpty()) {
                    ins = ins + extra;
                }
                String gt = "-" + len + "^" + ins;
                // TODO: 26/07/18 Seem that next statement is redundant since it assigns cnt to 0 if it is already 0
                // $hash->{ $softp }->{ $gt }->{ cnt } = 0 unless( $hash->{ $softp }->{ $gt }->{ cnt } );

                final Variation vref = getVariation(hash, softp, gt);
                inv.used = true;
                vref.pstd = true;
                vref.qstd = true;

                SV sv = getSV(hash, softp);
                sv.type = "INV";
                sv.splits += scv.cnt;
                sv.pairs += inv.cnt;
                sv.clusters++;

                Variation vrefSoftp = dir == -1
                        ? (hash.containsKey(softp) && ref.containsKey(softp) ? hash.get(softp).get(ref.get(softp).toString()) : null )
                        : null;
                adjCnt(vref, scv, vrefSoftp, conf);
                Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
                Map<String, Integer> map = new HashMap<>();
                map.put(gt, inv.cnt);
                dels5.put(softp, map);
                cov.put(softp, cov.containsKey(softp - 1) ? cov.get(softp - 1) : inv.cnt);
                scv.used = true;
                realigndel(hash, dels5, cov, sclip5, sclip3, reference, region, chrs, rlen, bams, conf);
                if (conf.y) {
                    System.err.printf(
                            "  Found INV SV: %s %d %s BP: %d cov: %d Cnt: %d EXTRA: %s %d %d %d cnt: %d %d\t DIR: %d Side: %s%n",
                            seq, softp, gt, bp, cov.get(softp), inv.cnt, extra, inv.mstart, inv.mend, inv.mlen, scv.cnt, (bp - softp) / abs(inv.mlen), dir, side
                    );
                }
                return vref;
            }
        }
        return null;
    }

    /**
     * Output remaining soft-clipped reads that haven't been used
     * @param sclip5 map of soft clips 5' on positions
     * @param sclip3 map of soft clips 3' on positions
     * @param conf Configuration
     */
    static void outputClipping(Map<Integer, Sclip> sclip5,
                                      Map<Integer, Sclip> sclip3,
                                      Configuration conf) {
        System.err.println("5' Remaining clipping reads");
        for (Map.Entry<Integer, Sclip> entry: sclip5.entrySet()) {
            int position_p = entry.getKey();
            Sclip sclip_sc = entry.getValue();

            if (sclip_sc.used) {
                continue;
            }
            if (sclip_sc.cnt < conf.minr) {
                continue;
            }
            String seq = findconseq(sclip_sc, conf, 0);
            if (!seq.isEmpty() && seq.length() > Configuration.SEED_2) {
                seq = new StringBuilder(seq).reverse().toString();
                System.err.printf("  P: %s Cnt: %s Seq: %s\n", position_p, sclip_sc.cnt, seq);
            }
        }
        System.err.println("3' Remaining clipping reads");
        for (Map.Entry<Integer, Sclip> entry: sclip3.entrySet()) {
            int position_p = entry.getKey();
            Sclip sclip_sc = entry.getValue();
            if (sclip_sc.used) {
                continue;
            }
            if (sclip_sc.cnt < conf.minr) {
                continue;
            }
            String seq = findconseq(sclip_sc, conf, 0);
            if (!seq.isEmpty() && seq.length() > Configuration.SEED_2) {
                System.err.printf("  P: %s Cnt: %s Seq: %s\n", position_p, sclip_sc.cnt, seq);
            }
        }
    }
}
