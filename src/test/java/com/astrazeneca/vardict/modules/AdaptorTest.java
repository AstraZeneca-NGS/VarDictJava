package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.VariationUtils;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class AdaptorTest {
    @AfterMethod
    public void cleanUp () {
        GlobalReadOnlyScope.clear();
    }

    @Test
    public void findConseqWithoutAdapter() {
        Configuration config = new Configuration();

        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), new HashMap<>());
        Sclip sclip = initSclip();
        String actualSequence = VariationUtils.findconseq(sclip, 0);
        String expectedSequence = "CTAAATC";

        Assert.assertEquals(actualSequence, expectedSequence);
    }

    @Test
    public void findConseqWithAdapter3Strand() {
        Configuration config = new Configuration();

        Map<String, Integer> adaptorForward = new HashMap<>();
        adaptorForward.put("CTAAAT", 1);

        GlobalReadOnlyScope.init(config, null, null, null, "", adaptorForward, new HashMap<>());
        Sclip sclip = initSclip();
        String actualSequence = VariationUtils.findconseq(sclip, 3);
        String expectedSequence = "";

        Assert.assertEquals(actualSequence, expectedSequence);
    }

    @Test
    public void findConseqWithAdapter5Strand() {
        Configuration config = new Configuration();

        Map<String, Integer> adaptorReverse = new HashMap<>();
        adaptorReverse.put("TAAATC", 1);

        GlobalReadOnlyScope.init(config, null, null, null, "", new HashMap<>(), adaptorReverse);
        Sclip sclip = initSclip();
        String actualSequence = VariationUtils.findconseq(sclip, 5);
        String expectedSequence = "";

        Assert.assertEquals(actualSequence, expectedSequence);
    }

    public Sclip initSclip() {
        Sclip sclip = new Sclip();
        sclip.varsCount = 2;
        sclip.nt.put(0, new TreeMap<Character, Integer>() {{ put('C', 2); }} );
        sclip.nt.put(1, new TreeMap<Character, Integer>() {{ put('T', 2); }} );
        sclip.nt.put(2, new TreeMap<Character, Integer>() {{ put('A', 2); }} );
        sclip.nt.put(3, new TreeMap<Character, Integer>() {{ put('A', 2); }} );
        sclip.nt.put(4, new TreeMap<Character, Integer>() {{ put('A', 2); }} );
        sclip.nt.put(5, new TreeMap<Character, Integer>() {{ put('T', 2); }} );
        sclip.nt.put(6, new TreeMap<Character, Integer>() {{ put('C', 2); }} );
        sclip.seq.put(0, new TreeMap<Character, Variation>() {{put('C', new Variation() {{ meanQuality=82.0; }}); }} );
        sclip.seq.put(1, new TreeMap<Character, Variation>() {{put('T', new Variation() {{ meanQuality=82.0; }}); }} );
        sclip.seq.put(2, new TreeMap<Character, Variation>() {{put('A', new Variation() {{ meanQuality=78.0; }}); }} );
        sclip.seq.put(3, new TreeMap<Character, Variation>() {{put('A', new Variation() {{ meanQuality=53.0; }}); }} );
        sclip.seq.put(4, new TreeMap<Character, Variation>() {{put('A', new Variation() {{ meanQuality=78.0; }}); }} );
        sclip.seq.put(5, new TreeMap<Character, Variation>() {{put('T', new Variation() {{ meanQuality=73.0; }}); }} );
        sclip.seq.put(6, new TreeMap<Character, Variation>() {{put('C', new Variation() {{ meanQuality=78.0; }}); }} );
        return sclip;
    }
}