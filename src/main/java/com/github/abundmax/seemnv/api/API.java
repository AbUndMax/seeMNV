package com.github.abundmax.seemnv.api;

import com.github.abundmax.seemnv.core.Core;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public final class API {

    public static void writeCMNVsToOutFile(Path vcfPath,
                                           Path codingSequencesPath,
                                           Path outputPath
    ) {

        String command = "Runned from seeMNV API within external Tool. Inputs: input-VCF-file=" +
                vcfPath.toString() + " coding-sequence-fasta=" + codingSequencesPath.toString() +
                " output-VCF-file=" + outputPath.toString();

        Core.runAlgorithm(vcfPath, codingSequencesPath, outputPath, command);

    }

    public static List<String> getCMNVEntries(Path vcfPath,
                                              Path codingSequencesPath
    ) throws IOException {

        return Core.getCMNVEntries(vcfPath, codingSequencesPath);

    }

    // Anbindung an HTSJDK um VCF objekte direkt zu übergeben und nicht über files arbeiten zu müssen
}
