package com.github.abundmax.seemnv.core;


import com.github.abundmax.seemnv.core.fieldArray.*;

import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.StringJoiner;

public class ParsedVCFEntry {

    // saves the FORMAT specification from the VCF-Entry
    private final VCFRecord vcfRecord;

    private final Genotype[] gtArray;
    private int ploidy;
    // Bitsets serving as Masks to calculate possible combinations later on
    private BitSet[] haplotypeTracks;

    private final HashMap<String, FieldArray<?>> fieldIDToArray = new HashMap<>();


    /**
     * Takes the a VCFrecord and creates a full functional Entry that parses the Sample Field values.
     * @param vcfRecord the record that should be processed to a functional Entry
     */
    public ParsedVCFEntry(VCFRecord vcfRecord) {
        this.vcfRecord =vcfRecord;

        boolean gtContained = false;
        for (String fieldID : vcfRecord.format()) {
            if (fieldID.equalsIgnoreCase("GT")) {
                gtContained = true;
                break;
            }
        }

        if (!gtContained) throw new IllegalArgumentException("VCFentry: " + vcfRecord.chrom() + vcfRecord.pos() + "does not contain GT values!");

        gtArray = new Genotype[vcfRecord.numSamples()];
        buildFieldArrays(vcfRecord.samples());
    }


    private void buildFieldArrays(String[] sampleColumns) {

        String[][] formattedSampleColumns = reformatSampleColumns(sampleColumns);

        for (int fieldIdx = 0; fieldIdx < numFields(); fieldIdx++) {
            String fieldID = vcfRecord.format()[fieldIdx];

            if (fieldID.equalsIgnoreCase("GT")) {
                fillGtArray(formattedSampleColumns[fieldIdx]);
                continue;
            }

            FieldConfig fieldConfig = Core.fieldToConfig.get(fieldID);
            FieldArray<?> fieldArray = FieldArrayFabricator.get(fieldConfig).apply(formattedSampleColumns[fieldIdx]);
            fieldIDToArray.put(fieldID, fieldArray);
        }
    }

    /**
     * Transfers the per sample information into a per field matrix:
     *
     * input:  [GT(sample1):AF(sample1):DP(sample1), GT(sample2):Af(sample2):DP(sample2), ...]
     * Output: [[GT(sample1), GT(sample2), ...], [AF(sample1), AF(sample2), ...], [DP(sample1), DP(sample2), ...]]
     *
     * Since th input follows the specified FORMAT order in the vcf entry, the order of Fields in the resulting
     * Matrix is the same.
     * @param sampleColumns
     */
    private String[][] reformatSampleColumns(String[] sampleColumns) {

        String[][] result = new String[numFields()][sampleColumns.length];

        for (int sampleIdx = 0; sampleIdx < sampleColumns.length; sampleIdx++) {
            String[] split = sampleColumns[sampleIdx].split(":");

            for (int fieldIdx = 0; fieldIdx < split.length; fieldIdx++) {
                result[fieldIdx][sampleIdx] = split[fieldIdx];
            }
        }

        return result;
    }

    /**
     * Populates the internal genotype array and initializes haplotype tracking for the samples.
     *
     * @param stringGtArray An array of genotype strings, where each element corresponds to a sample's genotype.
     */
    private void fillGtArray(String[] stringGtArray) {
        // set first Genotype and lookup the ploidy
        int firstIdx = 0;
        Genotype firstGT = new Genotype(stringGtArray[firstIdx]);
        gtArray[firstIdx] = firstGT;
        this.ploidy = firstGT.ploidy();

        // initialize BitSets
        haplotypeTracks = new BitSet[ploidy];
        for (int i = 0; i < ploidy; i++) {
            haplotypeTracks[i] = new BitSet(numSamples());
        }

        setHaploTracks(firstGT, firstIdx);


        for (int i = 1; i < numSamples(); i++) {
            Genotype genotype = new Genotype(stringGtArray[i]);
            gtArray[i] = genotype;
            setHaploTracks(genotype, i);
        }
    }

    /**
     * Sets the boolean value in the BitSets for the given genotype at the specified sample Index
     * @param genotype
     * @param sampleIdx
     */
    private void setHaploTracks(Genotype genotype, int sampleIdx) {
        Integer[] alleles = genotype.getAlleles();

        for (int alleleIdx = 0; alleleIdx < haplotypeTracks.length; alleleIdx++) {
            if (alleles[alleleIdx] != null && alleles[alleleIdx] > 0) {
                haplotypeTracks[alleleIdx].set(sampleIdx);
            }
        }
    }

    public FieldArray<?> getFieldArrayAt(String fieldID) {
        return fieldIDToArray.get(fieldID);
    }

    /**
     * Uses the formatID that is located at the index of formatIdx inside the FORMAT array to look up
     * the respective fieldArray associated with it.
     *
     * @param formatIdx The index inside FORMAT of the wanted formatID.
     * @return the {@code FieldArray<?>} of field FORMAT[formatIdx].
     */
    public FieldArray<?> getFieldArrayAt(int formatIdx) {
        String formatID = this.vcfRecord.format()[formatIdx];
        return fieldIDToArray.get(formatID);
    }

    public String[] getFORMAT() {
        return vcfRecord.format();
    }

    public String getCHROM() {
        return vcfRecord.chrom();
    }

    public long getPOS() {
        return vcfRecord.pos();
    }

    public String getID() {
        return vcfRecord.id();
    }

    public String getREF() {
        return vcfRecord.ref();
    }

    public String getALT() {
        return vcfRecord.alt();
    }

    public String getQUAL() {
        return vcfRecord.qual();
    }

    public String getFILTER() {
        return vcfRecord.filter();
    }

    public String getINFO() {
        return vcfRecord.info();
    }

    public Genotype[] getGtArray() {
        return gtArray;
    }

    // ideally copies should be returned, but to save on comp-power, teh direct value is returned
    // and not a deep copy.
    public BitSet[] getHaplotypeMasks() {
        return haplotypeTracks;
    }

    public int getPloidy() {
        return this.ploidy;
    }

    public int numSamples() {
        return vcfRecord.numSamples();
    }

    public int numFields() {
        return vcfRecord.numFields();
    }


    /**
     * Takes several {@code ParsedVCFEntries} and aggregates them to a MNVEntry.
     *
     * @param entries the {@code ParsedVCFEntries} to aggregate
     * @param codonPositions the affected codon positions
     * @param participationMasks sample participation masks per haplotype
     * @return the fully setup MNV entry
     */
    public static VCFRecord aggregateToMNV(
            ParsedVCFEntry[] entries,
            int[] codonPositions,
            BitSet[] participationMasks,
            String codonSequence,
            long codonIdx
    ) {
        if (entries == null || entries.length == 0) {
            throw new IllegalArgumentException("At least one VCFentry is required");
        }

        ParsedVCFEntry reference = entries[0];

        for (int i = 1; i < entries.length; i++) {
            // TODO: implement dynamic fallback for unordered FORMAT -> but all fields have to contain the same Fieldids
            if (!Arrays.equals(reference.getFORMAT(), entries[i].getFORMAT())) {
                throw new IllegalArgumentException("Cannot combine VCF entries with different FORMAT fields");
            }
        }

        String[] samples = combineEntriesIntoSampleStrings(entries, participationMasks);

        // Format: Num_alt_in_codon | Codon_idx_in_CDS | Original_codon | Variant_codon |QUAL_per_pos | FILTER_per_pos | affected_positions_genome | samples_with_alt_per_position| intersection_of_samples_with_alt
        // samples_with_alt_per_position -> actual number of samples showing the Alt at that position
        // intersection_of_samples_with_alt -> that number of samples which had ALL ALTs over the affected positions for this MNV
        int ploidy = participationMasks.length;


        // ------------------------------------ num_alt_in_codon
        int num_alt_in_codon = codonPositions.length;


        // ------------------------------------ variant_codon
        char[] variant_codon_prep = codonSequence.toCharArray();
        for (int i = 0; i < codonPositions.length; i++) {
            int codonPosition = codonPositions[i];
            variant_codon_prep[codonPosition] = entries[i].getALT().charAt(0);
        }
        String variant_codon = String.valueOf(variant_codon_prep);


        // ------------------------------------ qual_per_pos
        String[] qual_per_pos_prep = {"/", "/", "/"};
        for (int i = 0; i < codonPositions.length; i++) {
            int codonPosition = codonPositions[i];
            qual_per_pos_prep[codonPosition] = entries[i].getQUAL();
        }
        String qual_per_pos = Arrays.toString(qual_per_pos_prep);


        // ------------------------------------ filter_per_pos
        String[] filter_per_pos_prep = {"/", "/", "/"};
        for (int i = 0; i < codonPositions.length; i++) {
            int codonPosition = codonPositions[i];
            filter_per_pos_prep[codonPosition] = entries[i].getFILTER();
        }
        String filter_per_pos = Arrays.toString(filter_per_pos_prep);


        // ------------------------------------ affected_positions_in_genome
        String[] affected_positions_prep = {"/", "/", "/"};
        for (int i = 0; i < codonPositions.length; i++) {
            int codonPosition = codonPositions[i];
            affected_positions_prep[codonPosition] = String.valueOf(entries[i].getPOS());
        }
        String affected_positions_in_genome = Arrays.toString(affected_positions_prep);


        // ------------------------------------ samples_with_alt_per_position
        String[] samples_with_alt_per_position_prep = {"/", "/", "/"};
        for (int i = 0; i < codonPositions.length; i++) {
            int codonPosition = codonPositions[i];
            BitSet[] masks = entries[i].getHaplotypeMasks();

            StringJoiner joiner = new StringJoiner(", ", "[", "]");

            for (int j = 0; j < ploidy; j++) {
                joiner.add(String.valueOf(masks[j].cardinality()));
            }

            samples_with_alt_per_position_prep[codonPosition] = joiner.toString();
        }
        String samples_with_alt_per_position = Arrays.toString(samples_with_alt_per_position_prep);


        // ------------------------------------ intersection_of_samples_with_alt
        int[] intersection_of_samples_with_alt_prep = new int[ploidy];
        for (int i = 0; i < ploidy; i++) {
            intersection_of_samples_with_alt_prep[i] = participationMasks[i].cardinality();
        }
        String intersection_of_samples_with_alt = Arrays.toString(intersection_of_samples_with_alt_prep);


        // ------------------------------------ final INFO string
        String INFO = "cMNV=" + num_alt_in_codon + "|" + codonIdx + "|" + codonSequence + "|" + variant_codon + "|" +
                qual_per_pos + "|" + filter_per_pos + "|" + affected_positions_in_genome + "|" +
                samples_with_alt_per_position + "|" + intersection_of_samples_with_alt;


        return new VCFRecord(
                reference.getCHROM(),
                entries[0].getPOS(), // MNV entry shows first POS
                reference.getID(),
                codonSequence.substring(codonPositions[0], codonPositions[codonPositions.length - 1] + 1),
                buildALTFromEntries(entries),
                ".",
                ".",
                INFO,
                reference.getFORMAT(),
                samples
        );
    }

    /**
     * Combines multiple {@code ParsedVCFEntry} objects into a single string representation
     * that aggregates sample field values across the provided entries. Each field is processed
     * per sample, and relevant fields such as "GT" are handled with specific logic.
     *
     * @param entries an array of {@code ParsedVCFEntry} objects to combine. The first entry in
     *                the array is used as the reference for determining field structures.
     * @return a array of string where each combined value for a single sample is separated
     *         by a colon.
     */
    private static String[] combineEntriesIntoSampleStrings(ParsedVCFEntry[] entries, BitSet[] participationMasks) {
        ParsedVCFEntry reference = entries[0];

        // Holds arrays of Fields per sample:
        // [[AD(sample1), GT(sample1), DP(sample1)...], [AD(sample2), GT(sample2), DP(sample2)...], ...]
        String[][] perSampleFields = new String[reference.numSamples()][reference.numFields()];

        // for each field, collect the arrays of that field over all entries and combine them into a single one
        // save the combined values to its correct position in the perSampleFields as String
        for (int fieldIdx = 0; fieldIdx < reference.numFields(); fieldIdx++) {
            String fieldID = reference.getFORMAT()[fieldIdx];

            // extra logic for GT fields
            if (fieldID.equalsIgnoreCase("GT")) {

                // for each Sample, iterate over the VCF entries and combine them into a single
                // Genotype value for that sample, write it to perSampleFields at corresponding sampleIdx & fieldIdx
                for (int sampleIdx = 0; sampleIdx < reference.numSamples(); sampleIdx++) {
                    perSampleFields[sampleIdx][fieldIdx] = Genotype.toStringFromParticipationMasks(
                            participationMasks,
                            sampleIdx
                    );
                }

                continue;
            }

            FieldArray<?>[] toCombine = new FieldArray[entries.length];

            for (int entryIdx = 0; entryIdx < entries.length; entryIdx++) {
                toCombine[entryIdx] = entries[entryIdx].getFieldArrayAt(fieldID);
            }

            FieldArray<?> combined = FieldArray.combineWildcard(toCombine);

            for (int sampleIdx = 0; sampleIdx < reference.numSamples(); sampleIdx++) {
                perSampleFields[sampleIdx][fieldIdx] = combined.valueAtAsString(sampleIdx);
            }
        }

        String[] samples = new String[reference.numSamples()];

        for (int sampleIdx = 0; sampleIdx < reference.numSamples(); sampleIdx++) {
            samples[sampleIdx] = String.join(":", perSampleFields[sampleIdx]);
        }

        return samples;
    }

    private static String buildALTFromEntries(ParsedVCFEntry... entries) {
        StringBuilder ALT = new StringBuilder();

        for (ParsedVCFEntry entry : entries) {
            ALT.append(entry.getALT());
        }

        return ALT.toString();
    }
}
