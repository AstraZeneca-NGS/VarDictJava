package com.astrazeneca.vardict.integrationtests.utils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class OutputFileCreator {

    public static File fileCreator(String testcaseId, String moduleClass, String json) {
        File jsonFile = null;
        PrintWriter writerActual = null;
        try {
            jsonFile = new File(createDir(Paths.get("build/reports/" + testcaseId)).toString(),
                    moduleClass + ".txt");
            writerActual = new PrintWriter(jsonFile);
            writerActual.println(json);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (writerActual != null) {
                writerActual.close();
            }
        }
        return jsonFile;
    }


    private static Path createDir(Path path) throws IOException {
        if (path.endsWith("build/")) {
            return path;
        }
        createDir(path.getParent());
        return !Files.isDirectory(path) ? Files.createDirectory(path): path;
    }
}
