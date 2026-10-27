package com.github.abundmax.seemnv.core.fieldArray;

import java.util.BitSet;

public class Genotype implements Comparable<Genotype> {

    private final boolean[] phasing;
    private final Integer[] alleleNum;

    // implicit Phasing is allowed at the start and assumes:
    //  '/' if at least one '/' is contained, else '|'
    private char implicitPhaseIndicator = '0'; // '0' value means explicit phasing

    public Genotype(String gtFieldString) {

        boolean implicitPhasing = gtFieldString.charAt(0) != '|' && gtFieldString.charAt(0) != '/';
        int size = countNumbers(gtFieldString);

        if (implicitPhasing) {
            size++;
            implicitPhaseIndicator = '|';
        }

        phasing = new boolean[size];
        alleleNum = new Integer[size];

        int alleleIndex = 0;
        Integer number = 0;
        boolean hasDigit = false;

        for (int i = 0; i < gtFieldString.length(); i++) {
            char c = gtFieldString.charAt(i);

            if (c >= '0' && c <= '9') {
                number = number * 10 + (c - '0');
                hasDigit = true;

            } else if (c == '/' || c == '|') {
                if (hasDigit) {
                    alleleNum[alleleIndex] = number;
                    alleleIndex++;
                    number = 0;
                    hasDigit = false;
                }

                // suffice the implicaiton rule: if at least one / is contained, an unphased start is implied.
                if (implicitPhasing && c == '/') implicitPhaseIndicator = '/';
                phasing[alleleIndex] = c == '|';

            } else if (c == '.') {
                number = null;
                hasDigit = true;
            }
        }

        alleleNum[alleleIndex] = number;
        if (implicitPhasing) phasing[0] = implicitPhaseIndicator == '|';
    }

    private Genotype(boolean[] phasing, Integer[] alleleNum, char implicitPhaseIndicator) {
        this.phasing = phasing;
        this.alleleNum = alleleNum;
        this.implicitPhaseIndicator = implicitPhaseIndicator;
    }

    /**
     * Counts the numbers contained
     * @param input
     * @return
     */
    private static int countNumbers(String input) {
        int count = 0;

        for (int i = 0; i < input.length(); i++) {
            char c = input.charAt(i);

            if (c == '|' || c == '/') {
                count++;
            }
        }

        return count;
    }

    public int ploidy() {
        return this.alleleNum.length;
    }

    public Integer[] getAlleles() {
        return alleleNum;
    }

    @Override
    public int compareTo(Genotype o) {
        return 0;
    }

    @Override
    public String toString() {
       StringBuilder builder = new StringBuilder();
       for (int i = 0; i < phasing.length; i++) {
           char phase = '/';
           if (phasing[i]) phase = '|';
           builder.append(phase);
           String num = alleleNum[i] == null ? "." : alleleNum[i].toString();
           builder.append(num);
       }

       // if an implicit Phaseindicator at the first position is used, we delete the first.
       if (implicitPhaseIndicator != '0') builder.deleteCharAt(0);

       return builder.toString();
    }

    public static Genotype combine(Genotype... genotypes) {
        if (genotypes == null || genotypes.length == 0) {
            throw new IllegalArgumentException("At least one genotype is required");
        }

        Genotype reference = genotypes[0];
        int ploidy = reference.ploidy();

        for (Genotype genotype : genotypes) {
            if (genotype == null || genotype.ploidy() != ploidy) {
                throw new IllegalArgumentException("All genotypes must have the same ploidy");
            }
        }

        Integer[] combinedAlleles = new Integer[ploidy];

        for (int i = 0; i < ploidy; i++) {
            boolean allAlt = true;
            boolean hasMissing = false;

            for (Genotype genotype : genotypes) {
                Integer allele = genotype.alleleNum[i];

                if (allele == null) {
                    hasMissing = true;
                    break;
                }

                if (allele == 0) {
                    allAlt = false;
                }
            }

            if (hasMissing) {
                combinedAlleles[i] = null;
            } else {
                combinedAlleles[i] = allAlt ? 1 : 0;
            }
        }

        return new Genotype(reference.phasing.clone(), combinedAlleles, reference.implicitPhaseIndicator);
    }

    //TODO implement better phasing mechanisms
    public static String toStringFromParticipationMasks(
            BitSet[] participationMasks,
            int sampleIdx
    ) {
        StringBuilder builder = new StringBuilder(participationMasks.length * 2 - 1);

        for (int haplotypeIdx = 0; haplotypeIdx < participationMasks.length; haplotypeIdx++) {
            if (haplotypeIdx > 0) {
                builder.append('|');
            }

            builder.append(participationMasks[haplotypeIdx].get(sampleIdx) ? '1' : '0');
        }

        return builder.toString();
    }
}
