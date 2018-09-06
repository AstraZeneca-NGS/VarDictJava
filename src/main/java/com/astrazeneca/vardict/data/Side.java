package com.astrazeneca.vardict.data;

public enum Side {
    _3, _5, UNKNOWN;

    public static Side valueOf(int side) {
        return side == 3 ? _3 : side == 5 ? _5 : UNKNOWN;
    }
}

