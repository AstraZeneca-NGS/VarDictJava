package com.astrazeneca.vardict;

import java.io.IOException;

import com.astrazeneca.vardict.data.ReferenceResource;
import org.apache.commons.cli.*;

public class Main {
    /**
     * Method to build options from command line
     * @param args array of arguments from command line
     * @throws ParseException
     * @throws IOException
     */
    public static void main(String[] args) throws ParseException, IOException {
        Configuration config = CmdParser.parseParams(args);
        ReferenceResource referenceResource = new ReferenceResource();
        new VarDictLauncher(referenceResource).start(config);
    }
}
