package com.github.abundmax.seemnv.core;

import com.github.abundmax.seemnv.core.fieldArray.FieldConfig;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Core {

    // reading header fills fieldConfig Map -> FieldID to FieldConfig object
    public static final HashMap<String, FieldConfig> fieldToConfig = new HashMap<>();

    /**
     * Main loop that reads the input files, finds codonMNVs and writes them to the out file.
     *
     * @param vcfPath path to the vcfFile to be read from
     * @param codingSequencesPath path to the coding Sequence fasta file
     * @param outputPath
     */
    public static void runAlgorithm(Path vcfPath, Path codingSequencesPath, Path outputPath, String command) {

        String CALL_STRING = "##seeMNVCommand=" + command + "\n";
        String INFO_STRING = "##INFO=<ID=cMNV,Number=1,Type=String,Description=\"" +
                "Computed Multiple Nucleotide Variant within one Codon." +
                "Format: Num_alt_in_codon | Codon_idx_in_CDS | Original_codon" +
                "| Variant_codon | affected_positions | samples_with_alt_per_position" +
                "| intersection_of_samples_with_alt\">\n";

        try (
                BufferedReader vcfReader = Files.newBufferedReader(vcfPath);
                BufferedReader codingSeqReader = Files.newBufferedReader(codingSequencesPath);
                BufferedWriter vcfWriter = Files.newBufferedWriter(outputPath);
                ) {

            String vcfReaderLine;
            String codingSeqReaderLine;

            // parse Header lines and built fieldToConfig Map
            while ((vcfReaderLine = vcfReader.readLine()).startsWith("##")) {
                vcfWriter.write(vcfReaderLine);
                vcfWriter.newLine();

                if (vcfReaderLine.startsWith("##FORMAT")) {
                    FieldConfig fieldConfig = FieldConfig.buildFromHeaderString(vcfReaderLine);
                    if (fieldConfig != null) fieldToConfig.put(fieldConfig.id(), fieldConfig);
                }
            }

            vcfWriter.write(CALL_STRING);
            vcfWriter.write(INFO_STRING);
            // write column header #
            vcfWriter.write(vcfReaderLine);
            vcfWriter.newLine();





            // begin codon scanning and MNV detection
            long currentCDSStart = 0L;
            long currentCDSEnd = 0L;
            String currentSequenceHeader = codingSeqReader.readLine();
            String currentCDSSequence = "";

            Codon currentCodon = new Codon(0L, "init");

            while((vcfReaderLine = vcfReader.readLine()) != null) {

                // write original line into out VCF
                vcfWriter.write(vcfReaderLine);
                vcfWriter.newLine();

                String[] entry = vcfReaderLine.split("\t");
                String[] sampleColumns = Arrays.copyOfRange(entry, 9, entry.length);
                VCFRecord vcfRecord = new VCFRecord(
                        entry[0],
                        Long.parseLong(entry[1]),
                        entry[2],
                        entry[3],
                        entry[4],
                        entry[5],
                        entry[6],
                        entry[7],
                        entry[8].split(":"),
                        sampleColumns);

                // ignore frame shifts or continued deletions (* values in ALT column)
                if (vcfRecord.alt().equals("*") || vcfRecord.alt().length() > 1 || vcfRecord.pos() < currentCDSStart) {
                    continue;
                }




                // if the current VCFentry position exceeds the current CDS, read the next CDS
                if (currentCDSEnd < vcfRecord.pos()) {
                    Pattern pattern = Pattern.compile("\\[location=(\\d+)\\.\\.(\\d+)]");
                    Matcher matcher = pattern.matcher(currentSequenceHeader);

                    if (matcher.find()) {
                        currentCDSStart = Long.parseLong(matcher.group(1));
                        currentCDSEnd = Long.parseLong(matcher.group(2));
                    }

                    // read current CDS sequence
                    StringBuilder sequenceBuilder = new StringBuilder();
                    while((codingSeqReaderLine = codingSeqReader.readLine()) != null
                            && !codingSeqReaderLine.startsWith(">")) {
                        sequenceBuilder.append(codingSeqReaderLine);
                    }
                    currentCDSSequence = sequenceBuilder.toString();
                    currentSequenceHeader = codingSeqReaderLine;
                }




                long codonIdxOfRecord = (vcfRecord.pos() - currentCDSStart) / 3;

                // if a new Codon is reached, evaluate the current Codon in regard of MNVs, write them
                // then initialize a new Codon
                if (codonIdxOfRecord > currentCodon.getCodonIdx()) {

                    List<ParsedVCFEntry> mnvEntries = currentCodon.calculateMNVEntries();
                    for (ParsedVCFEntry mnvEntry : mnvEntries) {
                        vcfWriter.write(mnvEntry.toString());
                    }

                    // generate the codon sequence from the CDS sequence
                    String codonSeq = "";
                    int codonStart = Math.toIntExact(vcfRecord.pos() - currentCDSStart);
                    for (int i = codonStart; i < codonStart + 3; i++) {
                        codonSeq = codonSeq + currentCDSSequence.charAt(i);
                    }

                    currentCodon = new Codon(codonIdxOfRecord, codonSeq);
                }

                int inCodonPosition = Math.toIntExact((vcfRecord.pos() - currentCDSStart) % 3);
                currentCodon.addRecordAtCodonPos(inCodonPosition, vcfRecord);
            }

        } catch (IOException e) {
            // TODO
        }
    }

}
