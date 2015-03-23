package com.astrazeneca.vardict;

import java.util.*;

public final class Utils {

    private Utils() {}

    public static <E> String toString(Collection<E> collection) {
        Iterator<E> it = collection.iterator();
        if (! it.hasNext())
            return "";

        StringBuilder sb = new StringBuilder();
        for (;;) {
            E e = it.next();
            sb.append(e == collection ? "(this Collection)" : e);
            if (! it.hasNext())
                return sb.toString();
            sb.append(' ');
        }
    }

    public static String join(String delim, Object... args) {
        if (args.length == 0) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            sb.append(args[i]);
            if (i + 1 != args.length) {
                sb.append(delim);
            }
        }
        return sb.toString();
    }

    public static String joinNotNull(String delim, Object... args) {
        if (args.length == 0) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            if (args[i] == null) {
                if (i + 1 != args.length && args[i + 1] != null) {
                    sb.append(delim);
                }
                continue;
            }
            sb.append(args[i]);
            if (i + 1 != args.length && args[i + 1] != null) {
                sb.append(delim);
            }
        }
        return sb.toString();
    }

    public static <K, V> V getOrElse(Map<K, V> map, K key, V or) {
        V v = map.get(key);
        if (v == null) {
            v = or;
            map.put(key, v);
        }
        return v;
    }

    public static double round(double value, int dp) {
        double mf = Math.pow(10, dp);
        double d = value * mf;
        return Math.round(d) / mf;

    }

    public static int toInt(String intStr) {
        return Integer.parseInt(intStr);
    }

    public static String substr(String string, int idx) {
        if (idx >= 0) {
            return string.substring(idx);
        } else {
            return string.substring(Math.max(0, string.length() + idx));
        }
    }

    public static String substr(String string, int begin, int len) {
        if (begin < 0) {
            begin = string.length() + begin;
        }
        if (len > 0) {
            return string.substring(begin, Math.min(begin + len, string.length()));
        } else if (len == 0) {
            return "";
        } else {
            int end = string.length() + len;
            if (end < begin) {
                return "";
            }
            return string.substring(begin, end);
        }
    }

    public static char charAt(String str, int index) {
        if (index < 0) {
            int i = str.length() + index;
            if (i < 0)
                return (char)-1; //
            return str.charAt(i);
        }
        return str.charAt(index);
    }

    public static char charAt(StringBuilder str, int index) {
        if (index < 0) {
            int i = str.length() + index;
            if (i < 0)
                return (char)-1;//
            return str.charAt(str.length() + index);
        }
        return str.charAt(index);
    }

    public static int sum(Collection<?> list) {
        int result = 0;
        for (Object object : list) {
            result += toInt(String.valueOf(object));
        }
        return result;
    }

    public static List<String> globalFind(jregex.Pattern alignedLength, String string) {
        List<String> result = new LinkedList<>();
        jregex.Matcher matcher = alignedLength.matcher(string);
        while (matcher.find()) {
            result.add(matcher.group(1));
        }
        return result;

    }


}
