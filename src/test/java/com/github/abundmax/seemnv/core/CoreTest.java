package com.github.abundmax.seemnv.core;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class CoreTest {

    @TempDir
    Path tempDir;

    @BeforeEach
    void clearFieldConfigMap() {
        Core.fieldToConfig.clear();
    }

    @Test
    void testRunAlgorithmWritesExpectedHeaderLines() throws IOException {
        Path inputVcf = tempDir.resolve("input.vcf");
        Path codingSequences = tempDir.resolve("coding_sequences.fasta");
        Path outputVcf = tempDir.resolve("output.vcf");

        Files.writeString(inputVcf, minimalVcfWithSingleEntry(
                "NC_021508",
                "114",
                "A",
                "G"
        ));

        Files.writeString(codingSequences, minimalCodingSequence());

        Core.runAlgorithm(
                inputVcf,
                codingSequences,
                outputVcf,
                "seeMNV test command"
        );

        String output = Files.readString(outputVcf);

        assertTrue(output.contains("##fileformat=VCFv4.2"));
        assertTrue(output.contains("##seeMNVCommand=seeMNV test command"));
        assertTrue(output.contains("##INFO=<ID=cMNV"));
        assertTrue(output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1"));
    }

    @Test
    void testRunAlgorithmParsesFormatHeadersIntoFieldConfigMap() throws IOException {
        Path inputVcf = tempDir.resolve("input.vcf");
        Path codingSequences = tempDir.resolve("coding_sequences.fasta");
        Path outputVcf = tempDir.resolve("output.vcf");

        Files.writeString(inputVcf, minimalVcfWithSingleEntry(
                "NC_021508",
                "114",
                "A",
                "G"
        ));

        Files.writeString(codingSequences, minimalCodingSequence());

        Core.runAlgorithm(
                inputVcf,
                codingSequences,
                outputVcf,
                "seeMNV test command"
        );

        assertTrue(Core.fieldToConfig.containsKey("AD"));
        assertTrue(Core.fieldToConfig.containsKey("AF"));
        assertTrue(Core.fieldToConfig.containsKey("DP"));
        assertTrue(Core.fieldToConfig.containsKey("GQ"));
        assertTrue(Core.fieldToConfig.containsKey("PL"));
    }

    @Test
    void testRunAlgorithmWritesStarAltEntryUnchanged() throws IOException {
        Path inputVcf = tempDir.resolve("input.vcf");
        Path codingSequences = tempDir.resolve("coding_sequences.fasta");
        Path outputVcf = tempDir.resolve("output.vcf");

        Files.writeString(inputVcf, minimalVcfWithSingleEntry(
                "NC_021508",
                "114",
                "A",
                "*"
        ));

        Files.writeString(codingSequences, minimalCodingSequence());

        Core.runAlgorithm(
                inputVcf,
                codingSequences,
                outputVcf,
                "seeMNV test command"
        );

        String output = Files.readString(outputVcf);

        assertTrue(output.contains(
                "NC_021508\t114\t.\tA\t*\t1986.54\tPASS\tAC=1\tGT:AD:AF:DP:GQ:PL\t0:222,0:0:222:99:0,1800"
        ));
    }

    @Test
    void testRunAlgorithmWritesMultiBaseAltEntryUnchanged() throws IOException {
        Path inputVcf = tempDir.resolve("input.vcf");
        Path codingSequences = tempDir.resolve("coding_sequences.fasta");
        Path outputVcf = tempDir.resolve("output.vcf");

        Files.writeString(inputVcf, minimalVcfWithSingleEntry(
                "NC_021508",
                "114",
                "A",
                "GA"
        ));

        Files.writeString(codingSequences, minimalCodingSequence());

        Core.runAlgorithm(
                inputVcf,
                codingSequences,
                outputVcf,
                "seeMNV test command"
        );

        String output = Files.readString(outputVcf);

        assertTrue(output.contains(
                "NC_021508\t114\t.\tA\tGA\t1986.54\tPASS\tAC=1\tGT:AD:AF:DP:GQ:PL\t0:222,0:0:222:99:0,1800"
        ));
    }

    private String minimalVcfWithSingleEntry(String chrom, String pos, String ref, String alt) {
        return String.join("\n",
                           "##fileformat=VCFv4.2",
                           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                           "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
                           "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fraction\">",
                           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
                           "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
                           "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">",
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                           chrom + "\t" + pos + "\t.\t" + ref + "\t" + alt + "\t1986.54\tPASS\tAC=1\tGT:AD:AF:DP:GQ:PL\t0:222,0:0:222:99:0,1800"
        );
    }

    private String minimalCodingSequence() {
        return """
                >lcl|NC_021508.1_cds_WP_010881451.1_1 [gene=dnaA] [locus_tag=TPASS_RS00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=WP_010881451.1] [location=4..1398] [gbkey=CDS]
                ATGGACGCAGTAGGGTATGAAGTATTCTGGAACGAGACACTCAGCCAGATACGGAGTGAATCGACCGAAG
                CAGAATTTAACATGTGGTTTGCTCATTTGTTCTTTATCGCATCTTTTGAAAACGCTATCGAAATAGCAGT
                ACCTTCAGACTTTTTCCGAATACAGTTTAGCCAAAAATATCAAGAAAAGCTTGAGCGCAAGTTCCTCGAA
                CTTTCTGGACACCCCATTAAACTTTTGTTTGCCGTTAAAAAAGGCACCCCTCATGGAAATACTGCTCCCC
                CCAAACACGTGCATACCTACCTGGAGAAAAACTCTCCTGCAGAGGTTCCTTCCAAAAAGAGCTTTCACCC
                CGACCTGAACAGAGACTATACCTTCGAGAACTTTGTATCCGGAGAAGAAACCAAATTCAGCCATAGCGCT
                GCTATCTCCGTATCAAAAAACCCAGGCACTTCCTACAATCCGTTACTTATCTACGGTGGAGTGGGACTAG
                GAAAAACCCACCTTATGCAGGCTATTGGACACGAGATCTACAAGACAACAGACCTGAACGTCATATACGT
                CACTGCGGAGAATTTTGGAAATGAATTCATTTCCACATTACTCAATAAAAAGACCCAGGATTTTAAAAAA
                AAATACCGCTACACCGCGGATGTACTTCTTATAGATGACATTCATTTTTTTGAAAACAAAGACGGATTAC
                AAGAAGAGCTTTTCTATACGTTCAACGAACTTTTCGAGAAAAAAAAACAAATTATCTTTACCTGCGACAG
                GCCTGTACAAGAATTGAAAAATCTCTCTTCTCGCTTACGCTCGAGGTGCTCCCGAGGGCTTAGCACTGAT
                CTGAATATGCCATGTTTTGAAACGCGCTGTGCTATCTTGATTAAAAAAATACAAAACTATAACAGCACCT
                ATCCTCACAAAGCCATCCACATTTCAGACGATGTTGTCCGACTTGTTTCTGAAAACATTTCTTCAAATAT
                CAGGGATCTTGAGGGGGCATTAACAAAAATTATCGCTTTCATTGAAGTGTCGGGATCCATCACGATAGAT
                ATCGTTCCCTCTCTCCTAAAAGAGTTCTTCCTCTCTGCAAGGCCAAAACACATCACAGTAGAAACTATTC
                TTCATGTAGTTGCAGATCACTTTAACATTTCGTATTCAGATCTAAAGGGTAAGAAACGCAATAAAAGCGT
                TGTTTATCCTCGGCAAATCGCTATGTTTCTCTCAAAGGAACTGACAGAGCTCTCCACTACTGAACTTGGT
                ATCGAATTTGGTGGCAGAGATCATTCAACCGTCATTTACGGATGTCAAAAAATAGAAGGAGAAATTCTCA
                CTAATCCTTCGTTACAGGCAAATCTTGATTTGCTGAAAAGTAAAGTTCAAGATTCAATCCGCTAG
                >lcl|NC_021508.1_cds_WP_010881452.1_2 [gene=dnaN] [locus_tag=TPASS_RS00010] [protein=DNA polymerase III subunit beta] [protein_id=WP_010881452.1] [location=1641..2756] [gbkey=CDS]
                GTGAAAATCCTATGCGAGAAAGAAGCCTTTCTGAAGGAAATAAGCACAGCACAAGAGGTTATTTCAAATA
                AAAAAAACACGTCTATTTTTTCGAACGTCCTATTAGCTGCTCAAGGAGCCCTGCTTACCATCAGAGCAAC
                CGACACAAAAGTTACCTTTGAAACTAGCATTCCCGTCAATGTTCTCGCCGAAGGCACAACGACAGTTTTT
                TGCGACAAACTTGTGAATGTTGTTTCTGCCCTTCCAACAAAAGAAATCGAATTAACGTTATGTGAAGAAC
                AACTTGTCATTACCCCTCCAAACAAAAAGATAAGCTTTCAGCTCAGAACCCTCTCGCATGAGAGTTTTCC
                ATGTTTCCCTCAAAATGAAGGAGGCGTCTCTCTTGCTGTGCCTACCTCCGATCTTAGAAACATGATTAAC
                CATACCGTTTTTGCAGTTTCAGAAGACAGTACGCGCCATTTTATCAATGGCGTACACGTTGATTTTCAGT
                ATGGAAATATTATTTGTGTTTCAACAGATGGAAAGCGGCTCGCCTATATAGAAAAAAAGGGAGAATCCTC
                TCCCCAATCCTTTTCGGGTGTTATTGTGCCAACTAAGATCTTAGGCATAGTAAACCGTAAGCTTACCCCT
                GAAGGATCAGTGACGCTATGCATTACGTCGCAGCACGTTTACTTTTTTTTCGGTGGATATAAGTTTTCTT
                CTGTGCTTATTGAGGGGCAATTTCCTAATTACAAAAGAGTAATCCCTGATCATCAGGAGCGTTCTTTTTG
                TGTTGGACGTGTGGAGCTAATGGAGGCACTTAAACGAGTCTCGTTGTTGGTAGAACAAAAATCTCACAGG
                ATATTTATTACCATACAGCAGGGTTTGTTGACTTTAAGCTCAAAAGCTCACACTCAAGAAAATGAAATAG
                GTGATGCTCAGGAAGAAATAGCCTGTGCTTATACAGGAGAAAGTGAGGTCATAGCTCTTAACTATCTATA
                CCTTGAAGAACCGCTTAAGGTTTTTACTTCGAAGGAGGTTCAAGTGGAATTTACCGATCCTGCAAAAGCA
                CTCACGCTTCGTGCTGTACCAAACACGGACTGCTTTCACATCATTATGCCTATGCAAACGGAGTGA""";
    }
}