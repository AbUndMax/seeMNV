package com.github.abundmax.seemnv.core;

import com.github.abundmax.seemnv.core.fieldArray.FieldArray;
import com.github.abundmax.seemnv.core.fieldArray.Genotype;

import java.util.Arrays;
import java.util.BitSet;
import java.util.StringJoiner;

public final class MNVRecordBuilder {

    /**
     * Takes several {@code ParsedVCFEntries} and generates a MNVEntry from them.
     *
     * @param entries the {@code ParsedVCFEntries} to aggregate
     * @param codonPositions the affected codon positions
     * @param participationMasks sample participation masks per haplotype
     * @param codonSequence original codon sequence
     * @param codonIdx index of the codon in the CDS
     * @return the fully setup MNV entry
     */
    public static VCFRecord generate(
            ParsedVCFEntry[] entries,
            int[] codonPositions,
            BitSet[] participationMasks,
            String codonSequence,
            long codonIdx
    ) {
        validate(entries, codonPositions, participationMasks, codonSequence);

        ParsedVCFEntry reference = entries[0];
        String[] samples = combineEntriesIntoSampleStrings(entries, participationMasks);

        String ref = buildRef(codonSequence, codonPositions);
        String alt = buildAlt(codonSequence, entries, codonPositions);
        String info = buildInfo(entries, codonPositions, participationMasks, codonSequence, codonIdx);

        return new VCFRecord(
                reference.getCHROM(),
                entries[0].getPOS(),
                reference.getID(),
                ref,
                alt,
                ".",
                ".",
                info,
                reference.getFORMAT(),
                samples
        );
    }

    private static void validate(
            ParsedVCFEntry[] entries,
            int[] codonPositions,
            BitSet[] participationMasks,
            String codonSequence
    ) {
        if (entries == null || entries.length == 0) {
            throw new IllegalArgumentException("At least one ParsedVCFEntry is required");
        }

        if (codonPositions == null || codonPositions.length != entries.length) {
            throw new IllegalArgumentException("codonPositions must have the same length as entries");
        }

        if (participationMasks == null || participationMasks.length == 0) {
            throw new IllegalArgumentException("participationMasks must not be empty");
        }

        if (codonSequence == null || codonSequence.length() != 3) {
            throw new IllegalArgumentException("codonSequence must have length 3");
        }

        ParsedVCFEntry reference = entries[0];

        for (int i = 0; i < entries.length; i++) {
            if (entries[i] == null) {
                throw new IllegalArgumentException("entries must not contain null values");
            }

            if (codonPositions[i] < 0 || codonPositions[i] > 2) {
                throw new IllegalArgumentException("Codon positions must be between 0 and 2");
            }

            if (!Arrays.equals(reference.getFORMAT(), entries[i].getFORMAT())) {
                throw new IllegalArgumentException("Cannot combine VCF entries with different FORMAT fields");
            }
        }
    }

    private static String buildInfo(
            ParsedVCFEntry[] entries,
            int[] codonPositions,
            BitSet[] participationMasks,
            String codonSequence,
            long codonIdx
    ) {
        int numAltInCodon = codonPositions.length;
        String variantCodon = buildVariantCodon(codonSequence, entries, codonPositions);
        String qualPerPos = buildStringPerCodonPosition(entries, codonPositions, FieldValue.QUAL);
        String filterPerPos = buildStringPerCodonPosition(entries, codonPositions, FieldValue.FILTER);
        String affectedPositionsInGenome = buildStringPerCodonPosition(entries, codonPositions, FieldValue.POS);
        String samplesWithAltPerPosition = buildSamplesWithAltPerPosition(entries, codonPositions, participationMasks.length);
        String intersectionOfSamplesWithAlt = buildIntersectionOfSamplesWithAlt(participationMasks);

        return "cMNV=" + numAltInCodon + "|" +
                codonIdx + "|" +
                codonSequence + "|" +
                variantCodon + "|" +
                qualPerPos + "|" +
                filterPerPos + "|" +
                affectedPositionsInGenome + "|" +
                samplesWithAltPerPosition + "|" +
                intersectionOfSamplesWithAlt;
    }

    private static String buildVariantCodon(
            String codonSequence,
            ParsedVCFEntry[] entries,
            int[] codonPositions
    ) {
        char[] variantCodon = codonSequence.toCharArray();

        for (int i = 0; i < codonPositions.length; i++) {
            variantCodon[codonPositions[i]] = entries[i].getALT().charAt(0);
        }

        return String.valueOf(variantCodon);
    }

    private static String buildRef(String codonSequence, int[] codonPositions) {
        int start = arrayMin(codonPositions);
        int end = arrayMax(codonPositions);

        return codonSequence.substring(start, end + 1);
    }

    private static String buildAlt(
            String codonSequence,
            ParsedVCFEntry[] entries,
            int[] codonPositions
    ) {
        int start = arrayMin(codonPositions);
        int end = arrayMax(codonPositions);

        char[] alt = codonSequence.substring(start, end + 1).toCharArray();

        for (int i = 0; i < codonPositions.length; i++) {
            int relativePosition = codonPositions[i] - start;
            alt[relativePosition] = entries[i].getALT().charAt(0);
        }

        return String.valueOf(alt);
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

    private static String buildStringPerCodonPosition(
            ParsedVCFEntry[] entries,
            int[] codonPositions,
            FieldValue fieldValue
    ) {
        String[] result = {"/", "/", "/"};

        for (int i = 0; i < codonPositions.length; i++) {
            result[codonPositions[i]] = fieldValue.get(entries[i]);
        }

        return Arrays.toString(result);
    }

    private static String buildSamplesWithAltPerPosition(
            ParsedVCFEntry[] entries,
            int[] codonPositions,
            int ploidy
    ) {
        String[] result = {"/", "/", "/"};

        for (int i = 0; i < codonPositions.length; i++) {
            BitSet[] masks = entries[i].getHaplotypeMasks();
            StringJoiner joiner = new StringJoiner(", ", "[", "]");

            for (int haplotypeIdx = 0; haplotypeIdx < ploidy; haplotypeIdx++) {
                joiner.add(String.valueOf(masks[haplotypeIdx].cardinality()));
            }

            result[codonPositions[i]] = joiner.toString();
        }

        return Arrays.toString(result);
    }

    private static String buildIntersectionOfSamplesWithAlt(BitSet[] participationMasks) {
        int[] result = new int[participationMasks.length];

        for (int haplotypeIdx = 0; haplotypeIdx < participationMasks.length; haplotypeIdx++) {
            result[haplotypeIdx] = participationMasks[haplotypeIdx].cardinality();
        }

        return Arrays.toString(result);
    }

    private static int arrayMin(int[] values) {
        int min = values[0];

        for (int value : values) {
            if (value < min) {
                min = value;
            }
        }

        return min;
    }

    private static int arrayMax(int[] values) {
        int max = values[0];

        for (int value : values) {
            if (value > max) {
                max = value;
            }
        }

        return max;
    }

    private enum FieldValue {
        QUAL {
            @Override
            String get(ParsedVCFEntry entry) {
                return entry.getQUAL();
            }
        },
        FILTER {
            @Override
            String get(ParsedVCFEntry entry) {
                return entry.getFILTER();
            }
        },
        POS {
            @Override
            String get(ParsedVCFEntry entry) {
                return String.valueOf(entry.getPOS());
            }
        };

        abstract String get(ParsedVCFEntry entry);
    }
}