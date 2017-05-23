package com.astrazeneca.vardict.integrationtests.utils;

public class VarDictInput {

    private static final String DELIMITER = ";";
    private static final String REGION_DELIMITER = "-";
    private static final String CSV_DELIMITER = ",";

    public final String fasta;
    public final String bam;
    public final String chr;
    public final String start;
    public final String end;
    public final String args;
    public final String startAmp;
    public final String endAmp;
    public final String testName;

    public VarDictInput(String testName, String fasta, String bam, String chr, String start, String end, String startAmp, String endAmp, String args) {
        this.testName = testName;
        this.fasta = fasta;
        this.bam = bam;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.startAmp = startAmp;
        this.endAmp = endAmp;
        this.args = args;
    }

    @Override
    public String toString() {
        String[] bams = bam.split("\\|");
        return new StringBuilder()
                .append(testName).append(DELIMITER)
                .append(fasta).append(DELIMITER)
                .append(bams.length > 1 ? (bams[0] + REGION_DELIMITER + bams[1]) : bam).append(DELIMITER)
                .append(chr).append(DELIMITER)
                .append(start).append(REGION_DELIMITER)
                .append(end).append(DELIMITER)
                .append(startAmp == null ? "" : startAmp + REGION_DELIMITER)
                .append(endAmp == null ? "" : endAmp + DELIMITER)
                .append(args.length() == 0 ? "" : args + DELIMITER)
                .toString()
                .replaceAll(":", "");
    }

    public static VarDictInput fromCSVLine(String testCaseFromCSV) {
        try {
            String[] csvValues = testCaseFromCSV.split(CSV_DELIMITER);
            checkBamFileName(csvValues[2]);
            if (csvValues.length < 8) {
                return new VarDictInput(
                        csvValues[0],
                        csvValues[1],
                        csvValues[2],
                        csvValues[3],
                        csvValues[4],
                        csvValues[5],
                        null,
                        null,
                        csvValues.length == 7 ? csvValues[6] : ""
                );
            } else {
                return new VarDictInput(
                        csvValues[0],
                        csvValues[1],
                        csvValues[2],
                        csvValues[3],
                        csvValues[4],
                        csvValues[5],
                        csvValues[6],
                        csvValues[7],
                        csvValues.length == 9 ? csvValues[8] : ""
                );
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IllegalArgumentException("Wrong test case: " + testCaseFromCSV);
        }
    }

    private static void checkBamFileName(String bamName) {
        if (bamName.contains(":")) {
            throw new RuntimeException("Test cases should contains only one bam file except paired mode (bam file names delimited by \"|\")");
        }
    }
}
