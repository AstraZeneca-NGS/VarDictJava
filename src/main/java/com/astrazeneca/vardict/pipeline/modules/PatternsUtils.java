package com.astrazeneca.vardict.pipeline.modules;

import java.util.regex.Pattern;

/**
 * Created by Nikolai_Karulin on 1/30/2017.
 */
public class PatternsUtils {
    public static final Pattern BEGIN_DIGITS = Pattern.compile("^(\\d+)");
    public static final Pattern HASH_GROUP_CARET_GROUP = Pattern.compile("#(.+)\\^(.+)");
    public static final Pattern PLUS_NUMBER = Pattern.compile("\\+(\\d+)");
    public static final Pattern AMP_ATGC = Pattern.compile("&([ATGC]+)");
    public static final Pattern CARET_ATGNC = Pattern.compile("\\^([ATGNC]+)");
    public static final Pattern START_DIG = Pattern.compile("^-(\\d+)");
}
