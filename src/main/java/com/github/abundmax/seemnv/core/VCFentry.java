package com.github.abundmax.seemnv.core;


import com.github.abundmax.seemnv.core.fieldArray.FieldArray;
import com.github.abundmax.seemnv.core.fieldArray.FieldArrayFabricator;
import com.github.abundmax.seemnv.core.fieldArray.FieldConfig;
import com.github.abundmax.seemnv.core.fieldArray.Genotype;

import java.util.Arrays;

public class VCFentry {

    // saves the FORMAT specification from the VCF-Entry
    private final String CHROM;
    private final Integer POS;
    private final String ID;
    private final String REF;
    private final String ALT;
    private final String QUAL;
    private final String FILTER;
    private final String INFO;
    private final String[] FORMAT;
    private final String[] noGtFormat;

    private final Genotype[] gtArray;

    // they match exactkly the order of FORMAT
    // i.e. FORMAT= [GT, AD, FP] -> fieldArrays= [GTarray, ADarray, FParray]
    private final FieldArray<?>[] fieldArrays;


    /**
     * Takes the full entry string and splits it in its respective columns by setting the fields of this class
     * @param entryString
     */
    public VCFentry(String entryString) {

        // NC_021508	114	.	A	G	1986.54	PASS	AC=1;AF=0.0003937;AN=2539;DP=274532;FS=0;MQ=41.86;QD=25.47;SOR=0.855	GT:AD:AF:DP:GQ:PL	0:222,0:0:222:99:0,1800
        String[] entry = entryString.split("\t");
        CHROM = entry[0];
        POS = Integer.valueOf(entry[1]);
        ID = entry[2];
        REF = entry[3];
        ALT = entry[4];
        QUAL = entry[5];
        FILTER = entry[6];
        INFO = entry[7];
        FORMAT = entry[8].split(":");


        // GT arrays are held seperatly since most of the cMNV combination logic is based on them!
        String[] formatSplit = entry[8].split(":");

        // separat FORMAT contained that resembles the order of Fields in the FieldArray[]
        noGtFormat = new String[formatSplit.length - 1];

        boolean gtContained = false;
        int j = 0;
        for (int i = 0; i < formatSplit.length; i++) {
            if (formatSplit[i].equalsIgnoreCase("GT")) {
                gtContained = true;
                i++;
            }
            noGtFormat[j++] = formatSplit[i];
        }

        if (!gtContained) throw new IllegalArgumentException("VCFentry: " + CHROM + POS + "does not contain GT values!");

        String[] sampleColumns = Arrays.copyOfRange(entry, 9, entry.length);
        fieldArrays = new FieldArray[noGtFormat.length];
        gtArray = new Genotype[sampleColumns.length];
        buildFieldArrays(sampleColumns);
    }


    private void buildFieldArrays(String[] sampleColumns) {

        String[][] formattedSampleColumns = reformatSampleColumns(sampleColumns);

        int fieldArrayIdx = 0;
        for (int fieldIdx = 0; fieldIdx < this.FORMAT.length; fieldIdx++) {
            String fieldID = FORMAT[fieldIdx];

            if (fieldID.equalsIgnoreCase("GT")) {
                fillGtArray(formattedSampleColumns[fieldIdx]);
                continue;
            }

            FieldConfig fieldConfig = Core.fieldToConfig.get(fieldID);
            fieldArrays[fieldArrayIdx++] = FieldArrayFabricator.get(fieldConfig).apply(formattedSampleColumns[fieldIdx]);
        }

        //TODO...
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

        String[][] result = new String[FORMAT.length][sampleColumns.length];

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
        return FORMAT;
    }

    public String getCHROM() {
        return CHROM;
    }

    public Integer getPOS() {
        return POS;
    }

    public String getID() {
        return ID;
    }

    public String getREF() {
        return REF;
    }

    public String getALT() {
        return ALT;
    }

    public String getQUAL() {
        return QUAL;
    }

    public String getFILTER() {
        return FILTER;
    }

    public String getINFO() {
        return INFO;
    }

    public String[] getNoGtFormat() {
        return noGtFormat;
    }

    public Genotype[] getGtArray() {
        return gtArray;
    }

    public static VCFentry combine(VCFentry... entries) {
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



        VCFentry combinedEntry = new VCFentry("dummy");

        // TODO: set combinedEntry fields properly
        // combinedEntry.FORMAT = referenceFormat;
        // combinedEntry.fieldArrays = combinedFieldArrays;

        return combinedEntry;
    }

}
