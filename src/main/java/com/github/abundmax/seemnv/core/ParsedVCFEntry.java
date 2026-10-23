package com.github.abundmax.seemnv.core;


import com.github.abundmax.seemnv.core.fieldArray.*;

import java.util.Arrays;
import java.util.HashMap;

public class ParsedVCFEntry {

    // saves the FORMAT specification from the VCF-Entry
    private final VCFRecord vcfRecord;

    private final Genotype[] gtArray;

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

        for (int fieldIdx = 0; fieldIdx < vcfRecord.format().length; fieldIdx++) {
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

        String[][] result = new String[vcfRecord.format().length][sampleColumns.length];

        for (int sampleIdx = 0; sampleIdx < sampleColumns.length; sampleIdx++) {
            String[] split = sampleColumns[sampleIdx].split(":");

            for (int fieldIdx = 0; fieldIdx < split.length; fieldIdx++) {
                result[fieldIdx][sampleIdx] = split[fieldIdx];
            }
        }

        return result;
    }

    private void fillGtArray(String[] stringGtArray) {
        for (int i = 0; i < gtArray.length; i++) {
            gtArray[i] = new Genotype(stringGtArray[i]);
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

    public int numSamples() {
        return vcfRecord.numSamples();
    }

    public int numFields() {
        return vcfRecord.numFields();
    }


    /**
     * Takes several {@code ParsedVCFEntries} and aggregate them to a MNVEntry
     * @param entries the {@code ParsedVCFEntries} to aggregate
     * @return the fully setup MNV entry
     */
    public static VCFRecord aggregateToMNV(ParsedVCFEntry... entries) {
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

        String sampleString = combineEntriesIntoSampleString(entries);





        //TODO create new vcfRecord

        return null;
    }


    /**
     * Combines multiple {@code ParsedVCFEntry} objects into a single string representation
     * that aggregates sample field values across the provided entries. Each field is processed
     * per sample, and relevant fields such as "GT" are handled with specific logic.
     *
     * @param entries an array of {@code ParsedVCFEntry} objects to combine. The first entry in
     *                the array is used as the reference for determining field structures.
     * @return a tab-separated string where each combined value for a single sample is separated
     *         by a colon, and the aggregated values for all samples are concatenated.
     */
    private static String combineEntriesIntoSampleString(ParsedVCFEntry... entries) {
        // Combine FieldArrays

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
                    Genotype[] genotypeOverEntriesAtSample = new Genotype[entries.length];
                    for (int vcfEntry = 0; vcfEntry < genotypeOverEntriesAtSample.length; vcfEntry++) {
                        genotypeOverEntriesAtSample[vcfEntry] = entries[vcfEntry].getGtArray()[sampleIdx];
                    }
                    perSampleFields[sampleIdx][fieldIdx] = Genotype.combine(genotypeOverEntriesAtSample).toString();
                }
                continue;
            }

            // loop over VCF entries
            FieldArray<?>[] toCombine = new FieldArray[entries.length];
            for (int vcfEntry = 0; vcfEntry < toCombine.length; vcfEntry++) {
                toCombine[vcfEntry] = entries[vcfEntry].getFieldArrayAt(fieldID);
            }

            FieldArray<?> combined = FieldArray.combineWildcard(toCombine);
            for (int sampleIdx = 0; sampleIdx < reference.numSamples(); sampleIdx++) {
                perSampleFields[sampleIdx][fieldIdx] = combined.valueAtAsString(sampleIdx);
            }
        }

        // combine the separat values into a single String per Sample
        String[] singleStringPerSample = new String[reference.numSamples()];
        for (int sampleIdx = 0; sampleIdx < reference.numSamples(); sampleIdx++) {
            singleStringPerSample[sampleIdx] = String.join(":", perSampleFields[sampleIdx]);
        }

        return String.join("\t", singleStringPerSample);
    }

}
