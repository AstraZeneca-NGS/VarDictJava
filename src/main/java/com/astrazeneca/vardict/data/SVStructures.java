package com.astrazeneca.vardict.data;

import com.astrazeneca.vardict.variations.Sclip;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Structure for storage together all the collections needed for structural variants analysis.
 */
public class SVStructures {
    // for structural variant: Deletion
    public int svdelfend;
    public int svdelrend;
    public List<Sclip> svfdel = new ArrayList<>();
    public List<Sclip> svrdel = new ArrayList<>();

    // for structural variant: Duplication
    public int svdupfend;
    public int svduprend;
    public List<Sclip> svfdup = new ArrayList<>();
    public List<Sclip> svrdup = new ArrayList<>();

    // for structural variant: Inversion
    public int svinvfend3;
    public int svinvrend3;
    public int svinvfend5;
    public int svinvrend5;

    public List<Sclip> svfinv3 = new ArrayList<>();
    public List<Sclip> svfinv5 = new ArrayList<>();
    public List<Sclip> svrinv3 = new ArrayList<>();
    public List<Sclip> svrinv5 = new ArrayList<>();

    public Map<String, List<Sclip>> svffus = new HashMap<>();
    public Map<String, List<Sclip>> svrfus = new HashMap<>();
    public Map<String, Integer> svfusfend = new HashMap<>();
    public Map<String, Integer> svfusrend = new HashMap<>();

    // for structural variant: Insertion. Not used.
    // public int svinsfend = 0;
    // public int svinsrend = 0;
    // List<Sclip> svfins = new ArrayList<>();
    // List<Sclip> svrins = new ArrayList<>();
}
