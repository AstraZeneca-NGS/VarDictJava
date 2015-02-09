package com.epam;


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

        private final T1 _f;
        private final T2 _s;

        public Tuple2(T1 f, T2 s) {
            _f = f;
            _s = s;
        }

        public T1 _1() {
            return _f;
        }

        public T2 _2() {
            return _s;
        }

        public static <T1, T2> Tuple2<T1, T2> newTuple(T1 f, T2 s) {
            return new Tuple2<T1, T2>(f, s);
        }
    }

    public static class Tuple3<T1, T2, T3> {

        private final T1 _f;
        private final T2 _s;
        private final T3 _t;

        public Tuple3(T1 f, T2 s, T3 t) {
            _f = f;
            _s = s;
            _t = t;
        }

        public T1 _1() {
            return _f;
        }

        public T2 _2() {
            return _s;
        }

        public T3 _3() {
            return _t;
        }

        public static <T1, T2, T3> Tuple3<T1, T2, T3> newTuple(T1 f, T2 s, T3 t) {
            return new Tuple3<T1, T2, T3>(f, s, t);
        }
    }

    public static class Tuple4<T1, T2, T3, T4> {

        private final T1 _f;
        private final T2 _s;
        private final T3 _t;
        private final T4 _v;

        public Tuple4(T1 f, T2 s, T3 t, T4 v) {
            _f = f;
            _s = s;
            _t = t;
            _v = v;
        }

        public T1 _1() {
            return _f;
        }

        public T2 _2() {
            return _s;
        }

        public T3 _3() {
            return _t;
        }

        public T4 _4() {
            return _v;
        }

        public static <T1, T2, T3, T4> Tuple4<T1, T2, T3, T4> newTuple(T1 f, T2 s, T3 t, T4 v) {
            return new Tuple4<T1, T2, T3, T4>(f, s, t, v);
        }
    }
}
