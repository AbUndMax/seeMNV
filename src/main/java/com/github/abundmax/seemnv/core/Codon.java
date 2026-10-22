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

    public void addRecordAtCodonPos(int inCodonPosition, VCFRecord vcfRecord) {
        records.get(inCodonPosition).add(new ParsedVCFEntry(vcfRecord));
    }

    public List<ParsedVCFEntry> calculateMNVEntries() {
        ArrayList<ParsedVCFEntry> mnvEntries = new ArrayList<>();
        // TODO
        return mnvEntries;
    }
}
