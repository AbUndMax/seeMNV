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
     * Executes an algorithm to detect and process multiple in-codon nucleotide variants (cMNVs)
     * based on the provided VCF file, coding sequences, and command string. The output is written
     * to the specified output file.
     *
     * @param vcfPath Path to the input VCF file containing variant data.
     * @param codingSequencesPath Path to the file containing coding sequences.
     * @param outputPath Path to the output file where processed data will be written.
     * @param command The command string to annotate in the output file, specifying the method of analysis.
     */
    public static void runAlgorithm(Path vcfPath, Path codingSequencesPath, Path outputPath, String command) {

        String CALL_STRING = "##seeMNVCommand=" + command + "\n";
        String INFO_STRING = "##INFO=<ID=cMNV,Number=1,Type=String,Description=\"" +
                "Computed Multiple Nucleotide Variant within one Codon." +
                "Format: QUAL_per_pos | FILTER_per_pos | Num_alt_in_codon | Codon_idx_in_CDS | Original_codon" +
                "| Variant_codon | affected_positions_in_genome | samples_with_alt_per_position" +
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
                // TODO: check for same chromosome!!!!
                if (codonIdxOfRecord > currentCodon.getCodonIdx()) {

                    List<VCFRecord> mnvEntries = calculateMNVEntries(currentCodon);
                    for (VCFRecord mnvEntry : mnvEntries) {
                        vcfWriter.write(mnvEntry.toString());
                        vcfWriter.newLine();
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

    /**
     * generates all valid MNVentries from a codon.
     * @param codon
     * @return
     */
    public static List<VCFRecord> calculateMNVEntries(Codon codon) {
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

            mnvEntries.add(calculateSingleMNV(tripleCandidate, participationMasks, codon.getCodonSeq(), codon.getCodonIdx()));
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

            mnvEntries.add(calculateSingleMNV(pairCandidate, participationMasks, codon.getCodonSeq(), codon.getCodonIdx()));
        }

        return mnvEntries;
    }

    private record MNVParticipation(MNVCandidate candidate, BitSet[] participationMasks) {}

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

    /**
     * Calculates a MNV from a given candidate. Does no logical checking, only combines
     * it into a print-ready MNV.
     *
     * @param candidate the candidate to combine
     * @param participationMasks masks describing participating samples per haplotype
     * @return the aggregated MNV record
     */
    private static VCFRecord calculateSingleMNV(MNVCandidate candidate,
                                                BitSet[] participationMasks,
                                                String codonSequence,
                                                long codon_idx) {
        return ParsedVCFEntry.aggregateToMNV(candidate.entries(), candidate.codonPositions(), participationMasks, codonSequence, codon_idx);
    }

}
