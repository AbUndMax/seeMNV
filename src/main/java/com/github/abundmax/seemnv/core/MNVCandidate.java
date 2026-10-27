package com.github.abundmax.seemnv.core;

import java.util.Arrays;

public record MNVCandidate(ParsedVCFEntry[] entries, int[] codonPositions) {

    public int size() {
        return entries.length;
    }

    public boolean containsAll(MNVCandidate other) {
        for (ParsedVCFEntry otherEntry : other.entries) {
            boolean found = false;

            for (ParsedVCFEntry entry : entries) {
                if (entry == otherEntry) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                return false;
            }
        }

        return true;
    }

    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        }

        if (!(object instanceof MNVCandidate other)) {
            return false;
        }

        return sameEntriesByIdentity(entries, other.entries);
    }

    @Override
    public int hashCode() {
        int result = 1;

        for (ParsedVCFEntry entry : entries) {
            result = 31 * result + System.identityHashCode(entry);
        }

        return result;
    }

    private static boolean sameEntriesByIdentity(
            ParsedVCFEntry[] first,
            ParsedVCFEntry[] second
    ) {
        if (first.length != second.length) {
            return false;
        }

        for (int i = 0; i < first.length; i++) {
            if (first[i] != second[i]) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String toString() {
        return "MNVCandidate{" +
                "entries=" + Arrays.toString(entries) +
                ", codonPositions=" + Arrays.toString(codonPositions) +
                '}';
    }
}