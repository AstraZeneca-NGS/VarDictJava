package com.astrazeneca.vardict;

import java.io.IOException;
import java.util.*;

import static com.astrazeneca.utils.Utils.toInt;
import static java.util.Collections.singletonList;

class RegionBuilder {

    static final BedRowFormat DEFAULT_BED_ROW_FORMAT = new BedRowFormat(2, 6, 7, 9, 10, 12);
    private static final BedRowFormat CUSTOM_BED_ROW_FORMAT = new BedRowFormat(0, 1, 2, 3, 1, 2);
    private static final BedRowFormat AMP_BED_ROW_FORMAT = new BedRowFormat(0, 1, 2, 6, 7, 3);
    private final static Comparator<Region> INSERT_START_COMPARATOR = Comparator.comparingInt(o -> o.variantsStart);

    private static final String CHR_LABEL = "chr";

    private final Configuration config;
    private final Map<String, Integer> chromosomesLengths;

    RegionBuilder(Map<String, Integer> chromosomesLengths, Configuration config) {
        this.chromosomesLengths = chromosomesLengths;
        this.config = config;
    }

    List<List<Region>> buildRegions(List<String> segRaws, Boolean zeroBased) throws IOException {
        boolean isZeroBased = config.isZeroBasedDefined() ? config.zeroBased : false;
        List<List<Region>> segs = new LinkedList<>();
        BedRowFormat format = config.bedRowFormat;
        for (String seg : segRaws) {
            String[] split = seg.split(config.delimiter);
            if (!config.isColumnForChromosomeSet() && split.length == 4) {
                try {
                    int a1 = toInt(split[1]);
                    int a2 = toInt(split[2]);
                    if (a1 <= a2) {
                        format = CUSTOM_BED_ROW_FORMAT;
                        if (zeroBased == null) {
                            isZeroBased = true;
                        }
                    }
                } catch (NumberFormatException ignored) {}
            }
            String chr = split[format.chrColumn];
            chr = correctChromosome(chromosomesLengths, chr);
            int start = toInt(split[format.startColumn]);
            int end = toInt(split[format.endColumn]);
            String gene = format.geneColumn < split.length ? split[format.geneColumn] : chr;

            String[] thickStarts = split[format.thickStartColumn].split(",");
            String[] thickEnds = split[format.thickEndColumn].split(",");
            List<Region> thickRegions = new LinkedList<>();
            for (int i = 0; i < thickStarts.length; i++) {
                int thickStart = toInt(thickStarts[i]);
                int thickEnd = toInt(thickEnds[i]);
                if (start > thickEnd) {
                    continue;
                }
                if (end > thickEnd) {
                    break;
                }
                if (thickStart < start)
                    thickStart = start;
                if (thickEnd > end)
                    thickEnd = end;
                thickStart -= config.numberNucleotideToExtend;
                thickEnd += config.numberNucleotideToExtend;
                if (isZeroBased && thickStart < thickEnd) {
                    thickStart++;
                }
                thickRegions.add(new Region(chr, thickStart, thickEnd, gene));
            }
            segs.add(thickRegions);
        }
        return segs;
    }

    List<List<Region>> buildAmpRegions(List<String> segRaws, boolean zeroBased) throws IOException {
        List<List<Region>> segs = new LinkedList<>();
        Map<String, List<Region>> tsegs = new HashMap<>();
        for (String string : segRaws) {
            String[] split = string.split(config.delimiter);
            String chr = split[AMP_BED_ROW_FORMAT.chrColumn];
            chr = correctChromosome(chromosomesLengths, chr);
            int start = toInt(split[AMP_BED_ROW_FORMAT.startColumn]);
            int end = toInt(split[AMP_BED_ROW_FORMAT.endColumn]);
            String gene = split[AMP_BED_ROW_FORMAT.geneColumn];
            int insertStart = toInt(split[AMP_BED_ROW_FORMAT.thickStartColumn]);
            int insertEnd = toInt(split[AMP_BED_ROW_FORMAT.thickEndColumn]);
            if (zeroBased && start < end) {
                start++;
                insertStart++;
            }
            List<Region> chrRegions = tsegs.get(chr);
            if (chrRegions == null) {
                chrRegions = new ArrayList<>();
                tsegs.put(chr, chrRegions);
            }
            chrRegions.add(new Region(chr, start, end, gene, insertStart, insertEnd));
        }
        List<Region> regions = new LinkedList<>();
        segs.add(regions);
        int previousEnd = -1;
        for (Map.Entry<String, List<Region>> entry : tsegs.entrySet()) {
            List<Region> chrRegions = entry.getValue();
            chrRegions.sort(INSERT_START_COMPARATOR);
            String previousChr = null;
            for (Region region : chrRegions) {
                if (previousEnd != -1 && (!region.chr.equals(previousChr) || region.variantsStart > previousEnd)) {
                    regions = new LinkedList<>();
                    segs.add(regions);
                }
                regions.add(region);
                previousChr = region.chr;
                previousEnd = region.variantsEnd;
            }
        }

        return segs;
    }

    /**
     * Create region from command-line option
     */
    List<List<Region>> buildRegionFromConfiguration() {
        String[] split = config.regionOfInterest.split(":");
        String chr = split[0];
        chr = correctChromosome(chromosomesLengths, chr);
        String gene = split.length < 3 ? chr : split[2];
        String[] range = split[1].split("-");
        int start = toInt(range[0].replaceAll(",", ""));
        int end = range.length < 2 ? start : toInt(range[1].replaceAll(",", ""));
        start -= config.numberNucleotideToExtend;
        end += config.numberNucleotideToExtend;
        boolean zeroBased = config.isZeroBasedDefined() ? config.zeroBased : false;
        if (zeroBased && start < end) {
            start++;
        }
        if (start > end)
            start = end;

        return singletonList(singletonList(new Region(chr, start, end, gene)));
    }

    private static String correctChromosome(Map<String, Integer> chromosomesLengths, String chr) {
        if (!chromosomesLengths.containsKey(chr)) {
            if (chr.startsWith(CHR_LABEL)) {
                chr = chr.substring(CHR_LABEL.length());
            } else {
                chr = CHR_LABEL + chr;
            }
        }
        return chr;
    }

    static class BedRowFormat {
        final int chrColumn;
        final int startColumn;
        final int endColumn;
        final int thickStartColumn;
        final int thickEndColumn;
        final int geneColumn;

        BedRowFormat(int chrColumn, int startColumn, int endColumn, int thickStartColumn, int thickEndColumn, int geneColumn) {
            this.chrColumn = chrColumn;
            this.startColumn = startColumn;
            this.endColumn = endColumn;
            this.thickStartColumn = thickStartColumn;
            this.thickEndColumn = thickEndColumn;
            this.geneColumn = geneColumn;
        }
    }
}
