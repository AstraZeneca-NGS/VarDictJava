package com.astrazeneca.vardict;


public final class Tuple {

    private Tuple() { }

    public static <T1, T2> Tuple2<T1, T2> tuple(T1 f, T2 s) {
        return new Tuple2<T1, T2>(f, s);
    }

    public static <T1, T2, T3> Tuple3<T1, T2, T3> tuple(T1 f, T2 s, T3 t) {
        return new Tuple3<T1, T2, T3>(f, s, t);
    }

    public static <T1, T2, T3, T4> Tuple4<T1, T2, T3, T4> tuple(T1 f, T2 s, T3 t, T4 v) {
        return new Tuple4<T1, T2, T3, T4>(f, s, t, v);
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

    public static class Tuple4<T1, T2, T3, T4> {

        public final T1 _1;
        public final T2 _2;
        public final T3 _3;
        public final T4 _4;

        public Tuple4(T1 f, T2 s, T3 t, T4 v) {
            _1 = f;
            _2 = s;
            _3 = t;
            _4 = v;
        }

        public static <T1, T2, T3, T4> Tuple4<T1, T2, T3, T4> newTuple(T1 f, T2 s, T3 t, T4 v) {
            return new Tuple4<T1, T2, T3, T4>(f, s, t, v);
        }
    }
}
