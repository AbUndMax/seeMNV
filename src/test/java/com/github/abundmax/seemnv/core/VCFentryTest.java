package com.github.abundmax.seemnv.core;

import com.github.abundmax.seemnv.core.fieldArray.FieldArray;
import com.github.abundmax.seemnv.core.fieldArray.FieldConfig;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class VCFentryTest {

    @BeforeAll
    static void setup() {
        Core.fieldToConfig.put("AF", new FieldConfig("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">"));
        Core.fieldToConfig.put("DP", new FieldConfig("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">"));
        Core.fieldToConfig.put("GQ", new FieldConfig("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
        // Core.fieldToConfig.put("GT", new FieldConfig("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
        Core.fieldToConfig.put("PL", new FieldConfig("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"));
        Core.fieldToConfig.put("AD", new FieldConfig("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"));
    }

    @Test
    void testVCFentryCreation() {
        String chrom = "NC_021508";
        String pos = "114";
        String id = ".";
        String ref = "A";
        String alt = "G";
        String qual = "1986.54";
        String filter = "PASS";
        String info = "AC=1;AF=0.0003937;AN=2539;DP=274532;FS=0;MQ=41.86;QD=25.47;SOR=0.855";
        String format = "GT:AD:AF:DP:GQ:PL";

        String sample1 = "0:222,0:0:222:99:0,1800";
        String sample2 = "0:193,0:0:193:99:0,1350";
        String sample3 = "0:203,0:0:203:99:0,1350";
        String sample4 = "0:206,1:0.004831:207:99:0,1575";

        String vcfEntryString = String.join("\t", chrom, pos, id, ref, alt, qual, filter, info, format,
                                            sample1, sample2, sample3, sample4
        );

        VCFentry testEntry = new VCFentry(vcfEntryString);

        assertEquals(chrom, testEntry.getCHROM());
        assertEquals(Integer.valueOf(pos), testEntry.getPOS());
        assertEquals(id, testEntry.getID());
        assertEquals(ref, testEntry.getREF());
        assertEquals(alt, testEntry.getALT());
        assertEquals(qual, testEntry.getQUAL());
        assertEquals(filter, testEntry.getFILTER());
        assertEquals(info, testEntry.getINFO());

        assertArrayEquals(format.split(":"), testEntry.getFORMAT());
        assertArrayEquals(new String[]{"AD", "AF", "DP", "GQ", "PL"}, testEntry.getNoGtFormat());

        assertEquals(4, testEntry.getGtArray().length);
        assertEquals(sample1.split(":")[0], testEntry.getGtArray()[0].toString());
        assertEquals(sample2.split(":")[0], testEntry.getGtArray()[1].toString());
        assertEquals(sample3.split(":")[0], testEntry.getGtArray()[2].toString());
        assertEquals(sample4.split(":")[0], testEntry.getGtArray()[3].toString());

        assertEquals(5, testEntry.getFieldArrays().length);

        FieldArray<?> adArray = testEntry.getFieldArrays()[0];
        FieldArray<?> afArray = testEntry.getFieldArrays()[1];
        FieldArray<?> dpArray = testEntry.getFieldArrays()[2];
        FieldArray<?> gqArray = testEntry.getFieldArrays()[3];
        FieldArray<?> plArray = testEntry.getFieldArrays()[4];

        assertArrayEquals(new Integer[]{222, 0}, adArray.values[0]);
        assertArrayEquals(new Integer[]{193, 0}, adArray.values[1]);
        assertArrayEquals(new Integer[]{203, 0}, adArray.values[2]);
        assertArrayEquals(new Integer[]{206, 1}, adArray.values[3]);

        assertArrayEquals(new Float[]{0.0f}, afArray.values[0]);
        assertArrayEquals(new Float[]{0.0f}, afArray.values[1]);
        assertArrayEquals(new Float[]{0.0f}, afArray.values[2]);
        assertArrayEquals(new Float[]{0.004831f}, afArray.values[3]);

        assertArrayEquals(new Integer[]{222}, dpArray.values[0]);
        assertArrayEquals(new Integer[]{193}, dpArray.values[1]);
        assertArrayEquals(new Integer[]{203}, dpArray.values[2]);
        assertArrayEquals(new Integer[]{207}, dpArray.values[3]);

        assertArrayEquals(new Integer[]{99}, gqArray.values[0]);
        assertArrayEquals(new Integer[]{99}, gqArray.values[1]);
        assertArrayEquals(new Integer[]{99}, gqArray.values[2]);
        assertArrayEquals(new Integer[]{99}, gqArray.values[3]);

        assertArrayEquals(new Integer[]{0, 1800}, plArray.values[0]);
        assertArrayEquals(new Integer[]{0, 1350}, plArray.values[1]);
        assertArrayEquals(new Integer[]{0, 1350}, plArray.values[2]);
        assertArrayEquals(new Integer[]{0, 1575}, plArray.values[3]);
    }
}
