package com.astrazeneca.vardict;

public class Samrecord {

    final String qname; //0
    final Flag flag;//1
    final String rname;//2
    final int position;//3
    final int mapq;//4
    final String cigar; //5
    final String mrnm; //6
    final int mpos; //7
    final int isize; //8
    final String querySeq; //9
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

    public boolean isDefined(int num) {
        return num < recordSize;
    }


    public static class Flag {
        private int flag;

        public Flag(int flag) {
            this.flag = flag;
        }

        public boolean isSupplementaryAlignment() {
            return (flag & 0x800) != 0;
        }

        public boolean isUnmappedMate() {
            return (flag & 0x8) != 0;
        }

        public boolean isReverseStrand() {
            return (flag & 0x10) != 0;
        }

        public boolean isNotPrimaryAlignment() {
            return (flag & 0x100) != 0;
        }

    }


}