package com.astrazeneca.vardict.integrationtests.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class CSVReferenceManager {

    private static Map<String, CSVReferenceReader> readers = new HashMap<>();

    public static void init() {
        File directory = new File("testdata/fastas");
        for (File f : directory.listFiles()) {
            if (f.getPath().endsWith(".csv")) {
                readers.put(f.getName().replace(".csv", ""), new CSVReferenceReader(f.getAbsolutePath()));
            }
        }
    }

    private CSVReferenceManager() {}

    public static CSVReferenceReader getReader(String fastaFileName) {
        return readers.get(fastaFileName);
    }

    public static class CSVReferenceReader {

        private final String filePath;
        private final static String CSV_DELIMITER = ",";
        private Map<String, NavigableMap<Integer, String>> thinFasta;

        private CSVReferenceReader(String filePath) {
            this.filePath = filePath;
            thinFasta = readThinFastaFromFile();
        }

        private Map<String, NavigableMap<Integer, String>> readThinFastaFromFile() {
            Map<String, NavigableMap<Integer, String>> thinFasta = new HashMap<>();
            try (Scanner in = new Scanner(new File(filePath))) {
                while (in.hasNext()) {
                    String[] line = in.nextLine().split(CSV_DELIMITER);
                    NavigableMap<Integer, String> tupleMap = thinFasta.computeIfAbsent(line[0], k -> new TreeMap<>());
                    tupleMap.put(Integer.parseInt(line[1]), line[3]);
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            return thinFasta;
        }

        public String[] queryFasta(String chr, int start, int end) {
            String region = ">" + chr + ":" + start + "-" + end;
            NavigableMap<Integer, String> refByStart = thinFasta.get(chr);
            if (refByStart != null) {
                Map.Entry<Integer, String> entry = refByStart.floorEntry(start);
                if (entry != null) {
                    int csvStart = entry.getKey();
                    String seq = entry.getValue();
                    int csvEnd = csvStart + seq.length() - 1 ;
                    if (csvEnd >= end) {
                        return new String[]{region, getSubSeq(csvStart, csvEnd, start, end, seq)};
                    }
                }
            }
            throw new RuntimeException("No region: " + region + " in file " + filePath);
        }

        private String getSubSeq(int csvStart, int csvEnd, int start, int end, String seq) {
            int startPosition = csvStart == start ? 0 : start - csvStart;
            int endPosition = csvEnd == end ? seq.length() : seq.length() - csvEnd + end;
            return seq.substring(startPosition, endPosition);
        }
    }
}
