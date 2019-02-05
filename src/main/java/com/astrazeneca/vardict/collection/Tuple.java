package com.astrazeneca.vardict.collection;

/**
 * Util class for creating tuple structures of different size
 */
public final class Tuple {

    private Tuple() {
    }

    public static <T1, T2> Tuple2<T1, T2> tuple(T1 f, T2 s) {
        return new Tuple2<T1, T2>(f, s);
    }

    public static <T1, T2, T3> Tuple3<T1, T2, T3> tuple(T1 f, T2 s, T3 t) {
        return new Tuple3<T1, T2, T3>(f, s, t);
    }

    public static class Tuple2<T1, T2> {

        public final T1 _1;
        public final T2 _2;

        public Tuple2(T1 f, T2 s) {
            _1 = f;
            _2 = s;
        }

        public static <T1, T2> Tuple2<T1, T2> newTuple(T1 f, T2 s) {
            return new Tuple2<T1, T2>(f, s);
        }
    }

    public static class Tuple3<T1, T2, T3> {

        public final T1 _1;
        public final T2 _2;
        public final T3 _3;

        public Tuple3(T1 f, T2 s, T3 t) {
            _1 = f;
            _2 = s;
            _3 = t;
        }

        public static <T1, T2, T3> Tuple3<T1, T2, T3> newTuple(T1 f, T2 s, T3 t) {
            return new Tuple3<T1, T2, T3>(f, s, t);
        }
    }
}
