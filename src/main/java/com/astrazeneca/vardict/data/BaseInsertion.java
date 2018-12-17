package com.astrazeneca.vardict.data;

/**
 * Contains data about base insertion for realignment process (position of insertion, insertion sequence and position
 * without extra sequence)
 */
public class BaseInsertion {
    /**
     * Starting position of insert
     */
    public Integer baseInsert; //$bi
    /**
     * Insertion sequence
     */
    public String insertionSequence; //$ins
    /**
     * base position without extra sequence
     */
    public Integer baseInsert2; //$bi2

    public BaseInsertion(int baseInsert, String insertionSequence, int baseInsert2) {
        this.baseInsert = baseInsert;
        this.insertionSequence = insertionSequence;
        this.baseInsert2 = baseInsert2;
    }
}