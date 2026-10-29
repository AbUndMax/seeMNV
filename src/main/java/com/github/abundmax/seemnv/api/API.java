package com.github.abundmax.seemnv.api;

import com.github.abundmax.seemnv.core.Core;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public final class API {

    /**
     * Processes a VCF file and a coding sequences file to compute CMNV (Complex Multi-Nucleotide Variant) entries
     * and writes the output to a specified file.
     *
     * @param vcfPath The path to the input VCF (Variant Call Format) file containing genetic variants.
     * @param codingSequencesPath The path to the coding sequences file in FASTA format, used for sequence annotation.
     * @param outputPath The path to the output file where the processed CMNV data will be written.
     */
    public static void writeCMNVsToOutFile(Path vcfPath,
                                           Path codingSequencesPath,
                                           Path outputPath
    ) {

        String command = "Runned from seeMNV API within external Tool. Inputs: input-VCF-file=" +
                vcfPath.toString() + " coding-sequence-fasta=" + codingSequencesPath.toString() +
                " output-VCF-file=" + outputPath.toString();

        Core.runAlgorithm(vcfPath, codingSequencesPath, outputPath, command);

    }

    /**
     * Retrieves a list of CMNV (Complex Multi-Nucleotide Variant) entries by processing a VCF file
     * and a coding sequences file.
     *
     * @param vcfPath The path to the input VCF (Variant Call Format) file containing genetic variants.
     * @param codingSequencesPath The path to the coding sequences file in FASTA format, used for sequence annotation.
     * @return A list of strings representing the computed CMNV entries.
     * @throws IOException If an I/O error occurs while reading the input files.
     */
    public static List<String> getCMNVEntries(Path vcfPath,
                                              Path codingSequencesPath
    ) throws IOException {

        return Core.getCMNVEntries(vcfPath, codingSequencesPath);

    }

    // Anbindung an HTSJDK um VCF objekte direkt zu übergeben und nicht über files arbeiten zu müssen
}
