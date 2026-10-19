package com.github.abundmax.seemnv.core;


import com.github.abundmax.seemnv.core.fieldArray.FieldArray;
import com.github.abundmax.seemnv.core.fieldArray.FieldArrayFabricator;

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

        fieldArrays = new FieldArray[FORMAT.length];
        buildFieldArrays(Arrays.copyOfRange(entry, 9, entry.length));
    }

    private void buildFieldArrays(String[] sampleColumns) {

        for (int i = 0; i < this.FORMAT.length; i++) {
            String fieldID = FORMAT[i];
            // fieldArrays[i] = FieldArrayFabricator.get(Core.fieldToConfig.get(fieldID));
        }

        //TODO...
    }

    public FieldArray<?>[] getFieldArrays() {
        return fieldArrays;
    }

    public String[] getFORMAT() {
        return FORMAT;
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
