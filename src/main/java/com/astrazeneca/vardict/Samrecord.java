package com.astrazeneca.vardict;

/**
 * Class for holding aligned segment data from SAM/BAM file
 */
public class Samrecord {

    /**
     * Query template name
      */
    final String qname; //0

    /**
     * Bitwise flag
     */
    final Flag flag;//1

    /**
     * Reference sequence name
     */
    final String rname;//2

    /**
     * 1-based leftmost mapping position of the first matching base. Set to 0 for unmapped read without coordinate
     */
    final int position;//3

    /**
     * Mapping quality. Equal to -10*log_10 Pr(mapping position is wrong), rounded to the nearest integer.
     * Value 255 indicates that the mapping quality is not available.
     */
    final int mapq;//4

    /**
     * CIGAR string (see Samtools documentation)
     */
    final String cigar; //5

    /**
     * Reference sequence name of the primary alignment of the NEXT read in the template.
     * Set to '*' when information is not available, and to '=' if identical to rname.
     */
    final String mrnm; //6

    /**
     * Position of the primary alignment of the NEXT read in the template. Set as 0 when the information is unavailable.
     */
    final int mpos; //7

    /**
     * Signed observed Template LENgth (see Samtools documentation).
     */
    final int isize; //8

    /**
     * Segment sequence
     */
    final String querySeq; //9

    /**
     * ASCII of base QUALity plus 33 (same as the quality string in the Sanger FASTQ format).
     * A base quality is the phred-scaled base error probability which equals -10*log10(Pr{base is wrong})
     * This field can be a '*' when quality is not stored.
     */
    final String queryQual; //10

    int recordSize;

    public Samrecord(String record) {
        String[] a = record.split("\t");
        recordSize = a.length;

        qname = a[0];
        flag = new Flag(Integer.parseInt(a[1]));
        rname = a.length > 2 ? a[2] : null;
        position = a.length > 3 ? Integer.parseInt(a[3]) : 0;
        mapq = a.length > 4 ? Integer.parseInt(a[4]) : 0;
        cigar = a.length > 5 ? a[5] : "";
        mrnm = a.length > 6 ? a[6] : "";
        mpos = a.length > 7 ? Integer.parseInt(a[7]) : 0;
        isize = a.length > 8 ? Integer.parseInt(a[8]) : 0;
        querySeq = a.length > 9 ? a[9] : "";
        queryQual = a.length > 10 ? a[10] : "";
    }

    /**
     * Is record field defined
     * @param num - record index
     * @return Is record defined
     */
    public boolean isDefined(int num) {
        return num < recordSize;
    }


    /**
     * Class for holding bitwise flag in SAM/BAM record
     */
    public static class Flag {
        private int flag;

        public Flag(int flag) {
            this.flag = flag;
        }

        /**
         * @return indicates that the corresponding alignment line is part of a chimeric alignment.
         */
        public boolean isSupplementaryAlignment() {
            return (flag & 0x800) != 0;
        }

        /**
         *
         * @return next segment in the template is unmapped
         */
        public boolean isUnmappedMate() {
            return (flag & 0x8) != 0;
        }

        /**
         *
         * @return sequence is reverse complemented
         */
        public boolean isReverseStrand() {
            return (flag & 0x10) != 0;
        }

        /**
         * Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use are aware of this bit.
         * It is typically used to flag alternative mappings when multiple mappings are presented in a SAM.
         * @return Is secondary alignment
         */
        public boolean isNotPrimaryAlignment() {
            return (flag & 0x100) != 0;
        }

    }


}