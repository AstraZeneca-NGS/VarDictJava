package com.astrazeneca.vardict;

import org.apache.commons.cli.ParseException;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws ParseException, IOException {
        Configuration config = CmdParser.parseParams(args);
        new VarDictLauncher(config).start();
    }
}
