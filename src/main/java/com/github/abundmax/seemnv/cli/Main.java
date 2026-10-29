package com.github.abundmax.seemnv.cli;
import ArgsParser.ArgsParser;
import ArgsParser.ParameterTypes.PthParameter;
import com.github.abundmax.seemnv.core.Core;


public final class Main {

    static void main(String[] args) {

        ArgsParser parser = new ArgsParser("Identify and note in-codon MNVs in VCF files");
        PthParameter vcfPath = parser.addParameter(new PthParameter("vcfPath",
                                             "vp",
                                             "Path to the VCF file.",
                                             true, true));
        PthParameter codingSeqPath = parser.addParameter(new PthParameter("codingSequencePath",
                                                                          "csp",
                                                                          "Path to the coding Sequences FASTA file.",
                                                                          true, true));
        PthParameter outPath = parser.addParameter(new PthParameter("outPath",
                                             "o",
                                             "Path where the output VCF file should be saved.",
                                             true, false));

        parser.parse(args);

        String command = String.join(" ", args);

        System.out.println("Running cMNV detection!");
        Core.runAlgorithm(vcfPath.getArgument(), codingSeqPath.getArgument(), outPath.getArgument(), command);
    }
}