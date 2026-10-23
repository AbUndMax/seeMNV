package com.github.abundmax.seemnv.core;

import com.github.abundmax.seemnv.core.fieldArray.FieldConfig;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.time.Duration;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTimeoutPreemptively;

public class CodonTest {

    @BeforeEach
    void setup() {
        Core.fieldToConfig.clear();
        Core.fieldToConfig.put(
                "AD",
                FieldConfig.buildFromHeaderString(
                        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
                )
        );
    }

    @Test
    public void testGenerateAllALTCombinationsWithOneEntryAtEachCodonPosition() {
        assertTimeoutPreemptively(Duration.ofMillis(200), () -> {
            Codon codon = new Codon(0L, "ATG");
            VCFRecord dummyRecord = dummyRecord();

            codon.addRecordAtCodonPos(0, dummyRecord);
            codon.addRecordAtCodonPos(1, dummyRecord);
            codon.addRecordAtCodonPos(2, dummyRecord);

            List<ParsedVCFEntry[]> combinations = codon.generateAllALTCombinations();

            assertEquals(1, combinations.size());

            ParsedVCFEntry[] combination = combinations.getFirst();

            assertEquals(3, combination.length);
            assertNotNull(combination[0]);
            assertNotNull(combination[1]);
            assertNotNull(combination[2]);
        });
    }

    @Test
    public void testGenerateAllALTCombinationsWithMissingMiddleCodonPosition() {
        assertTimeoutPreemptively(Duration.ofMillis(200), () -> {
            Codon codon = new Codon(0L, "ATG");
            VCFRecord dummyRecord = dummyRecord();

            codon.addRecordAtCodonPos(0, dummyRecord);
            codon.addRecordAtCodonPos(2, dummyRecord);

            List<ParsedVCFEntry[]> combinations = codon.generateAllALTCombinations();

            assertEquals(1, combinations.size());

            ParsedVCFEntry[] combination = combinations.getFirst();

            assertEquals(3, combination.length);
            assertNotNull(combination[0]);
            assertNull(combination[1]);
            assertNotNull(combination[2]);
        });
    }

    @Test
    public void testGenerateAllALTCombinationsWithMultipleEntriesAtOneCodonPosition() {
        assertTimeoutPreemptively(Duration.ofMillis(200), () -> {
            Codon codon = new Codon(0L, "ATG");
            VCFRecord dummyRecord = dummyRecord();

            codon.addRecordAtCodonPos(0, dummyRecord);
            codon.addRecordAtCodonPos(0, dummyRecord);
            codon.addRecordAtCodonPos(1, dummyRecord);
            codon.addRecordAtCodonPos(2, dummyRecord);

            List<ParsedVCFEntry[]> combinations = codon.generateAllALTCombinations();

            assertEquals(2, combinations.size());

            ParsedVCFEntry[] firstCombination = combinations.get(0);
            ParsedVCFEntry[] secondCombination = combinations.get(1);

            assertEquals(3, firstCombination.length);
            assertNotNull(firstCombination[0]);
            assertNotNull(firstCombination[1]);
            assertNotNull(firstCombination[2]);

            assertEquals(3, secondCombination.length);
            assertNotNull(secondCombination[0]);
            assertNotNull(secondCombination[1]);
            assertNotNull(secondCombination[2]);
        });
    }

    private VCFRecord dummyRecord() {
        return new VCFRecord(
                "NC_021508",
                114L,
                ".",
                "A",
                "G",
                "1986.54",
                "PASS",
                "AC=1",
                new String[]{"GT", "AD"},
                new String[]{"0:222,0"}
        );
    }
}