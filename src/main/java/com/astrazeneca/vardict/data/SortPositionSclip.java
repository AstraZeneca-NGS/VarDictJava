package com.astrazeneca.vardict.data;

import com.astrazeneca.vardict.variations.Sclip;

/**
 * Temp structure from hash tables of insertion or deletions positions to description string
 * and counts of variation.
 */
public class SortPositionSclip {
    public int position;
    /**
     * description string for variant
     */
    public Sclip softClip;
    /**
     * SoftClip count
     */
    public int count;

    public SortPositionSclip(int position, Sclip softClip, int count) {
        this.position = position;
        this.softClip = softClip;
        this.count = count;
    }
}
