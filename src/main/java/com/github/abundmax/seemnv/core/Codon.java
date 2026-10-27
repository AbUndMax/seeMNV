package com.github.abundmax.seemnv.core;

import java.util.ArrayList;
import java.util.List;

public class Codon {

    private final long codonIdx;
    private final String codonSeq;

    // Per position multiple Entries can be
    private final List<List<ParsedVCFEntry>> records = List.of(
            new ArrayList<>(), // VCFEnries in first position
            new ArrayList<>(), // VCFEnries in second position
            new ArrayList<>()  // VCFEnries in third position
    );

    public Codon(long codonIdx, String codonSeq) {
        this.codonIdx = codonIdx;
        this.codonSeq = codonSeq;
    }

    public long getCodonIdx() {
        return codonIdx;
    }

    public String getCodonSeq() {
        return codonSeq;
    }

    public void addRecordAtCodonPos(int inCodonPosition, VCFRecord vcfRecord) {
        records.get(inCodonPosition).add(new ParsedVCFEntry(vcfRecord));
    }

    public List<MNVCandidate> generateMNVCandidates() {
        List<MNVCandidate> candidates = new ArrayList<>();

        addTripleCandidates(candidates);

        addPairCandidates(candidates, 0, 1);
        addPairCandidates(candidates, 0, 2);
        addPairCandidates(candidates, 1, 2);

        return candidates;
    }

    /**
     * Generates MNVs of length 3 if all positions in the Codon has an VCFentry!
     * Elsewise the loop won't run and thus return nothing.
     * @param candidates
     */
    private void addTripleCandidates(List<MNVCandidate> candidates) {
        for (ParsedVCFEntry first : records.get(0)) {
            for (ParsedVCFEntry second : records.get(1)) {
                for (ParsedVCFEntry third : records.get(2)) {
                    candidates.add(new MNVCandidate(
                            new ParsedVCFEntry[]{first, second, third},
                            new int[]{0, 1, 2}
                    ));
                }
            }
        }
    }

    private void addPairCandidates(
            List<MNVCandidate> candidates,
            int firstPosition,
            int secondPosition
    ) {
        for (ParsedVCFEntry first : records.get(firstPosition)) {
            for (ParsedVCFEntry second : records.get(secondPosition)) {
                candidates.add(new MNVCandidate(
                        new ParsedVCFEntry[]{first, second},
                        new int[]{firstPosition, secondPosition}
                ));
            }
        }
    }

    /**
     * Returns a list of length 3 arrays. Each entry in the list corresponds to one possible
     * combination of {@code ParsedVCFEntry} over the three positions in this codon.
     *
     * <p>Example:</p>
     * <pre>
     * A codon with:
     * pos 1: ParsedVCFEntry#1, ParsedVCFEntry#2
     * pos 2: /
     * pos 3: ParsedVCFEntry#3
     *
     * will return the combinations:
     * [[ParsedVCFEntry#1, null, ParsedVCFEntry#3],
     * [ParsedVCFEntry#2, null, ParsedVCFEntry#3]]
     * </pre>
     *
     * @return a list of possible {@code ParsedVCFEntry} combinations
     */
    public List<ParsedVCFEntry[]> generateAllALTCombinations() {

        List<ParsedVCFEntry[]> combinations = new ArrayList<>();

        if (records.get(0).isEmpty()
                && records.get(1).isEmpty()
                && records.get(2).isEmpty()) {
            return combinations;
        }

        int firstPositionOptions = Math.max(1, records.get(0).size());
        int secondPositionOptions = Math.max(1, records.get(1).size());
        int thirdPositionOptions = Math.max(1, records.get(2).size());

        for (int i = 0; i < firstPositionOptions; i++) {
            for (int j = 0; j < secondPositionOptions; j++) {
                for (int k = 0; k < thirdPositionOptions; k++) {
                    ParsedVCFEntry[] combination = new ParsedVCFEntry[]{
                            this.getOrNull(0, i),
                            this.getOrNull(1, j),
                            this.getOrNull(2, k)
                    };

                    combinations.add(combination);
                }
            }
        }

        return combinations;
    }

    /**
     * Returns the {@code ParsedVCFEntry} at position posIdx with ALT variant at ALTIdx or null
     * if no {@code ParsedVCFEntry} is found at that position.
     * @param posIdx the position a {@code ParsedVCFEntry} should be looked up
     * @param ALTIdx the ALT variant index of that {@code ParsedVCFEntry}.
     * @return the {@code ParsedVCFEntry} or null.
     */
    private ParsedVCFEntry getOrNull(int posIdx, int ALTIdx) {
        if (records.get(posIdx).isEmpty()) {
            return null;
        }

        return records.get(posIdx).get(ALTIdx);
    }
}
