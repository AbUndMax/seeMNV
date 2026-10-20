package com.github.abundmax.seemnv.core.fieldArray;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class GenotypeTests {

    @Test
    public void GenotypeCreationTest() {
        String gtField = "0/12/.";

        Genotype gt = new Genotype(gtField);

        assertEquals(gtField, gt.toString());
    }

    @Test
    public void GenotypePhasingAlterationTest() {
        String gtField = "0/.|13";
        Genotype gt = new Genotype(gtField);
        assertEquals(gtField, gt.toString());

        String gtField2 = "|2/1|0";
        Genotype gt2 = new Genotype(gtField2);
        assertEquals(gtField2, gt2.toString());

        String gtField3 = "2|1|0";
        Genotype gt3 = new Genotype(gtField3);
        assertEquals(gtField3, gt3.toString());
    }
}
