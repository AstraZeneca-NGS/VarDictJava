package com.astrazeneca.vardict.variations;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Variants for position
 */
public class Vars {
    /**
     * Reference variant
     */
    public Variant referenceVariant;

    /**
     * List of all variants except reference variant
     */
    public List<Variant> variants = new ArrayList<>();

    /**
     * Map of all variants except reference variant.
     * Key - variant description string, value - variant
     */
    public Map<String, Variant> varDescriptionStringToVariants = new HashMap<>();

    public String sv = "";
}
