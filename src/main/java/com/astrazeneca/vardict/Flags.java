package com.astrazeneca.vardict;

/**
 * Class for holding bitwise flag in SAM/BAM record
 */
public class Flags {
    private int flag;

    public Flags(int flag) {
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