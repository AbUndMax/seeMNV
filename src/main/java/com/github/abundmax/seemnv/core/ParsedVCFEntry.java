package com.github.abundmax.seemnv.core;


import com.github.abundmax.seemnv.core.fieldArray.*;

import java.util.Arrays;

public class ParsedVCFEntry {

    // saves the FORMAT specification from the VCF-Entry
    private final VCFRecord vcfRecord;
    private final String[] noGtFormat;

    private final Genotype[] gtArray;

    // they match exactkly the order of FORMAT
    // i.e. FORMAT= [GT, AD, FP] -> fieldArrays= [GTarray, ADarray, FParray]
    private final FieldArray<?>[] fieldArrays;


    /**
     * Takes the a VCFrecord and creates a full functional Entry that parses the Sample Field values.
     * @param vcfRecord the record that should be processed to a functional Entry
     */
    public ParsedVCFEntry(VCFRecord vcfRecord) {
        this.vcfRecord =vcfRecord;

        // separat FORMAT contained that resembles the order of Fields in the FieldArray[]
        noGtFormat = new String[vcfRecord.format().length - 1];

        boolean gtContained = false;
        int j = 0;
        for (int i = 0; i < vcfRecord.format().length; i++) {
            if (vcfRecord.format()[i].equalsIgnoreCase("GT")) {
                gtContained = true;
                i++;
            }
            noGtFormat[j++] = vcfRecord.format()[i];
        }

        if (!gtContained) throw new IllegalArgumentException("VCFentry: " + vcfRecord.chrom() + vcfRecord.pos() + "does not contain GT values!");

        fieldArrays = new FieldArray[noGtFormat.length];
        gtArray = new Genotype[vcfRecord.samples().length];
        buildFieldArrays(vcfRecord.samples());
    }


    private void buildFieldArrays(String[] sampleColumns) {

        String[][] formattedSampleColumns = reformatSampleColumns(sampleColumns);

        int fieldArrayIdx = 0;
        for (int fieldIdx = 0; fieldIdx < vcfRecord.format().length; fieldIdx++) {
            String fieldID = vcfRecord.format()[fieldIdx];

            if (fieldID.equalsIgnoreCase("GT")) {
                fillGtArray(formattedSampleColumns[fieldIdx]);
                continue;
            }

            FieldConfig fieldConfig = Core.fieldToConfig.get(fieldID);
            fieldArrays[fieldArrayIdx++] = FieldArrayFabricator.get(fieldConfig).apply(formattedSampleColumns[fieldIdx]);
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

    public FieldArray<?>[] getFieldArrays() {
        return fieldArrays;
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

    public String[] getNoGtFormat() {
        return noGtFormat;
    }

    public Genotype[] getGtArray() {
        return gtArray;
    }

    public static ParsedVCFEntry combine(ParsedVCFEntry... entries) {
        if (entries == null || entries.length == 0) {
            throw new IllegalArgumentException("At least one VCFentry is required");
        }

        String[] referenceFormat = entries[0].getFORMAT();

        for (int i = 1; i < entries.length; i++) {
            // TODO: implement dynamic fallback for unordered FORMAT -> but all fields have to be the same
            if (!Arrays.equals(referenceFormat, entries[i].getFORMAT())) {
                throw new IllegalArgumentException("Cannot combine VCF entries with different FORMAT fields");
            }
        }

        FieldArray<?>[] combinedFieldArrays = new FieldArray[referenceFormat.length];

        // loop over FORMAT fields / FieldArrays
        for (int i = 0; i < referenceFormat.length; i++) {

            FieldArray<?> base = entries[0].getFieldArrays()[i];
            FieldArray<?>[] others = new FieldArray[entries.length - 1];

            // loop over remaining VCF entries
            for (int j = 1; j < entries.length; j++) {
                others[j - 1] = entries[j].getFieldArrays()[i];
            }

            // combinedFieldArrays[i] = base.combineWithWildcard(others);
        }



        // VCFentry combinedEntry = new VCFentry();

        // TODO: set combinedEntry fields properly
        // combinedEntry.FORMAT = referenceFormat;
        // combinedEntry.fieldArrays = combinedFieldArrays;

        return null; // combinedEntry;
    }

}
