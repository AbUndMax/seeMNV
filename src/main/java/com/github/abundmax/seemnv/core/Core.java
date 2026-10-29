package com.github.abundmax.seemnv.core;

import com.github.abundmax.seemnv.core.fieldArray.FieldConfig;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Core {

    // reading header fills fieldConfig Map -> FieldID to FieldConfig object
    public static final HashMap<String, FieldConfig> fieldToConfig = new HashMap<>();

    /**
     * Functional interface representing a consumer for writing lines of text.
     * Designed to abstract the writing mechanism, allowing for more flexible
     * handling of line-based outputs.
     * <p>
     * The method signature includes an {@code IOException} to account for
     * the possibility of IO-related issues during the writing process.
     */
    @FunctionalInterface
    private interface LineWriter {
        void write(String line) throws IOException;
    }

    /**
     * Reads, processes, and analyzes a VCF (Variant Call Format) file along with coding sequence
     * data to detect and compute Multiple Nucleotide Variants (MNVs) that occur within a codon.
     * Outputs the processed information to a specified output file.
     * <p>
     *
     * <b>THE INPUT FILE REMAINS UNTOUCHED!</b>
     *
     * @param vcfPath              Path to the input VCF file containing variant information.
     * @param codingSequencesPath  Path to the file containing coding sequence data used for analysis.
     * @param outputPath           Path to the output file where processed data is written.
     * @param command              Command string to describe the execution context stored in the output file metadata.
     */
    public static void runAlgorithm(
            Path vcfPath,
            Path codingSequencesPath,
            Path outputPath,
            String command
    ) {

        String CALL_STRING = "##seeMNVCommand=" + command + "\n";
        String INFO_STRING = "##INFO=<ID=cMNV,Number=1,Type=String,Description=\"" +
                "Computed Multiple Nucleotide Variant within one Codon." +
                "Format: Num_alt_in_codon | Codon_idx_in_CDS | Original_codon | Variant_codon " +
                "| QUAL_per_pos | FILTER_per_pos | affected_positions_in_genome | samples_with_alt_per_position" +
                "| intersection_of_samples_with_alt\">\n";

        try (
                BufferedReader vcfReader = Files.newBufferedReader(vcfPath);
                BufferedReader codingSeqReader = Files.newBufferedReader(codingSequencesPath);
                BufferedWriter vcfWriter = Files.newBufferedWriter(outputPath)
        ) {

            String vcfReaderLine;

            // parse Header lines and built fieldToConfig Map
            while ((vcfReaderLine = vcfReader.readLine()) != null && vcfReaderLine.startsWith("##")) {
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



            processVariants(
                    vcfReader,
                    codingSeqReader,
                    line -> {
                        vcfWriter.write(line);
                        vcfWriter.newLine();
                    },
                    line -> {
                        vcfWriter.write(line);
                        vcfWriter.newLine();
                    }
            );

        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    /**
     * Reads and processes a VCF (Variant Call Format) file and coding sequence data to identify
     * and extract codon multi-nucleotide variation (cMNV) entries. This method analyzes the
     * provided data and returns a list of entries as strings.
     *
     * @param vcfPath            Path to the VCF file that contains variant data.
     * @param codingSequencesPath Path to the coding sequences file used for analysis.
     * @return A list of strings representing the identified cMNV entries.
     * @throws IOException If an I/O error occurs while reading the input files.
     */
    public static List<String> getCMNVEntries(
            Path vcfPath,
            Path codingSequencesPath
    ) throws IOException {
        List<String> entries = new ArrayList<>();

        try (
                BufferedReader vcfReader = Files.newBufferedReader(vcfPath);
                BufferedReader codingSeqReader = Files.newBufferedReader(codingSequencesPath)
        ) {
            String vcfReaderLine;

            // Skip VCF metadata header lines.
            while ((vcfReaderLine = vcfReader.readLine()) != null && vcfReaderLine.startsWith("##")) {
                // Header line skipped.
            }

            processVariants(vcfReader,
                            codingSeqReader,
                            _ -> {},
                            entries::add);

        }

        return entries;
    }

    /**
     * Processes VCF (Variant Call Format) and coding sequence data to detect and handle
     * codon variations and multi-nucleotide variations (MNVs). This method reads, analyzes,
     * and writes processed data to the provided consumers for further use.
     *
     * @param vcfReader             BufferedReader for reading VCF file entries line by line.
     * @param codingSeqReader       BufferedReader for reading coding sequences line by line.
     * @param originalLineConsumer  LineWriter functional interface for consuming original VCF lines.
     * @param cMNVConsumer          LineWriter functional interface for consuming computed MNV entries.
     * @throws IOException          If an I/O error occurs during reading or writing operations.
     */
    private static void processVariants(
            BufferedReader vcfReader,
            BufferedReader codingSeqReader,
            LineWriter originalLineConsumer,
            LineWriter cMNVConsumer
    ) throws IOException {

        String vcfReaderLine;
        String codingSeqReaderLine;

        // begin codon scanning and MNV detection
        long currentCDSStart = 0L;
        long currentCDSEnd = 0L;
        String currentSequenceHeader = codingSeqReader.readLine();
        String currentCDSSequence = "";

        Codon currentCodon = new Codon(0L, "init");

        while((vcfReaderLine = vcfReader.readLine()) != null) {

            // write original line into out VCF
            originalLineConsumer.write(vcfReaderLine);

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
            // TODO: check for same chromosome!!!!
            if (codonIdxOfRecord > currentCodon.getCodonIdx()) {

                List<VCFRecord> mnvEntries = calculateMNVEntriesInCodon(currentCodon);
                for (VCFRecord mnvEntry : mnvEntries) {
                    cMNVConsumer.write(mnvEntry.toString());
                }

                // generate the codon sequence from the CDS sequence
                StringBuilder codonSeq = new StringBuilder();
                int codonStart = Math.toIntExact(vcfRecord.pos() - currentCDSStart);
                for (int i = codonStart; i < codonStart + 3; i++) {
                    codonSeq.append(currentCDSSequence.charAt(i));
                }

                currentCodon = new Codon(codonIdxOfRecord, codonSeq.toString());
            }

            int inCodonPosition = Math.toIntExact((vcfRecord.pos() - currentCDSStart) % 3);
            currentCodon.addRecordAtCodonPos(inCodonPosition, vcfRecord);
        }

        // eval last codon
        List<VCFRecord> mnvEntries = calculateMNVEntriesInCodon(currentCodon);

        for (VCFRecord mnvEntry : mnvEntries) {
            cMNVConsumer.write(mnvEntry.toString());
        }
    }

    /**
     * Represents the participation of haplotypes in a given multiple nucleotide variation (MNV) candidate.
     * An MNV refers to a group of nucleotide changes that occur together within a codon.
     * <p>
     * This record encapsulates two components:
     * 1. The {@link MNVCandidate} object, representing the MNV candidate that includes the relevant
     *    VCF (Variant Call Format) entries and codon positions.
     * 2. An array of {@link BitSet} objects, referred to as participation masks. These masks identify
     *    the contribution of individual haplotypes (samples) to the MNV candidate by indicating
     *    which samples exhibit the alteration (ALT) across all relevant entries.
     * <p>
     * Participation is primarily determined through the use of the `calculateParticipationMasks` method
     * from the containing class. Each {@link BitSet} in the array represents one haplotype's participation.
     */
    private record MNVParticipation(MNVCandidate candidate, BitSet[] participationMasks) {}

    /**
     * generates all valid MNVentries from a codon.
     * @param codon
     * @return
     */
    public static List<VCFRecord> calculateMNVEntriesInCodon(Codon codon) {
        ArrayList<VCFRecord> mnvEntries = new ArrayList<>();

        List<MNVCandidate> candidates = codon.generateMNVCandidates();

        if (candidates.isEmpty()) {
            return mnvEntries;
        }

        ArrayList<MNVCandidate> tripleCandidates = new ArrayList<>();
        ArrayList<MNVCandidate> pairCandidates = new ArrayList<>();

        for (MNVCandidate candidate : candidates) {
            if (candidate.size() == 3) {
                tripleCandidates.add(candidate);
            } else if (candidate.size() == 2) {
                pairCandidates.add(candidate);
            }
        }

        ArrayList<MNVParticipation> validTripleParticipations = new ArrayList<>();

        for (MNVCandidate tripleCandidate : tripleCandidates) {
            int ploidy = tripleCandidate.entries()[0].getPloidy();

            BitSet[] participationMasks = calculateParticipationMasks(tripleCandidate.entries(), ploidy);

            if (!hasAnyParticipation(participationMasks)) {
                continue;
            }

            MNVParticipation participation = new MNVParticipation(tripleCandidate, participationMasks);

            validTripleParticipations.add(participation);

            mnvEntries.add(MNVRecordBuilder.generate(tripleCandidate.entries(),
                                                     tripleCandidate.codonPositions(),
                                                     participationMasks,
                                                     codon.getCodonSeq(),
                                                     codon.getCodonIdx())
            );
        }

        for (MNVCandidate pairCandidate : pairCandidates) {
            int ploidy = pairCandidate.entries()[0].getPloidy();

            BitSet[] participationMasks = calculateParticipationMasks(pairCandidate.entries(), ploidy);

            for (MNVParticipation tripleParticipation : validTripleParticipations) {
                if (tripleParticipation.candidate().containsAll(pairCandidate)) {
                    subtractParticipationMasksInPlace(participationMasks, tripleParticipation.participationMasks());
                }
            }

            if (!hasAnyParticipation(participationMasks)) {
                continue;
            }

            mnvEntries.add(MNVRecordBuilder.generate(pairCandidate.entries(),
                                                     pairCandidate.codonPositions(),
                                                     participationMasks,
                                                     codon.getCodonSeq(),
                                                     codon.getCodonIdx())
            );
        }

        return mnvEntries;
    }

    /**
     * Calculates per-haplotype participation masks for a given combination of VCF entries.
     * Each returned BitSet contains the samples that carry ALT on all given entries
     * for the corresponding haplotype index.
     *
     * @param entries the entries to test
     * @param ploidy the expected ploidy
     * @return one BitSet per haplotype
     */
    private static BitSet[] calculateParticipationMasks(ParsedVCFEntry[] entries, int ploidy) {
        // TODO: include phasing logic!
        BitSet[] result = new BitSet[ploidy];

        for (int haplotypeIdx = 0; haplotypeIdx < ploidy; haplotypeIdx++) {
            BitSet intersection = (BitSet) entries[0].getHaplotypeMasks()[haplotypeIdx].clone();

            for (int entryIdx = 1; entryIdx < entries.length; entryIdx++) {
                intersection.and(entries[entryIdx].getHaplotypeMasks()[haplotypeIdx]);

                if (intersection.isEmpty()) {
                    break;
                }
            }

            result[haplotypeIdx] = intersection;
        }

        return result;
    }

    private static void subtractParticipationMasksInPlace(BitSet[] masks, BitSet[] masksToSubtract) {
        for (int haplotypeIdx = 0; haplotypeIdx < masks.length; haplotypeIdx++) {
            masks[haplotypeIdx].andNot(masksToSubtract[haplotypeIdx]);
        }
    }

    private static boolean hasAnyParticipation(BitSet[] masks) {
        for (BitSet mask : masks) {
            if (!mask.isEmpty()) {
                return true;
            }
        }

        return false;
    }
}
