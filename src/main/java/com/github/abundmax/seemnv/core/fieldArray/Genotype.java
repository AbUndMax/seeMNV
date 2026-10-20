package com.github.abundmax.seemnv.core.fieldArray;

public class Genotype implements Comparable<Genotype> {

    private final boolean[] phasing;
    private final Integer[] alleleNum;

    // implicit Phasing is allowed at the start and assumes:
    //  '/' if at least one '/' is contained, else '|'
    private char implicitPhaseIndicator = '0';

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

    /**
     * Returns an {@code Genotype} object holding the value of the specified {@code String}.
     * @param string the string to be parsed.
     * @return an {@code Genotype} object
     */
   public Genotype valueOf(String string) {
       return new Genotype(string);
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
}
