package com.github.abundmax.seemnv.core;


import com.github.abundmax.seemnv.core.fieldArray.*;

import java.util.BitSet;
import java.util.HashMap;

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
}
