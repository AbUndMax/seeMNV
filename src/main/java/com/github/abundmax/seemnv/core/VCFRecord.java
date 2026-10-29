package com.github.abundmax.seemnv.core;

import java.util.Arrays;

// eg. NC_021508	114	.	A	G	1986.54	PASS	AC=1;AF=0.0003937;AN=2539;DP=274532;FS=0;MQ=41.86;QD=25.47;SOR=0.855	GT:AD:AF:DP:GQ:PL	0:222,0:0:222:99:0,1800
public record VCFRecord(
        String chrom,
        long pos,
        String id,
        String ref,
        String alt,
        String qual,
        String filter,
        String info,
        String[] format,
        String[] samples
) {
    public int numSamples() {
        return samples.length;
    }

    public int numFields() {
        return format().length;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();

        builder.append(chrom).append('\t')
                .append(pos).append('\t')
                .append(id).append('\t')
                .append(ref).append('\t')
                .append(alt).append('\t')
                .append(qual).append('\t')
                .append(filter).append('\t')
                .append(info).append('\t')
                .append(String.join(":", format));

        for (String sample : samples) {
            builder.append('\t').append(sample);
        }

        return builder.toString();
    }

    public static VCFRecord fromString(String vcfEntryLine) {
        String[] entry = vcfEntryLine.split("\t");
        String[] sampleColumns = Arrays.copyOfRange(entry, 9, entry.length);
        return new VCFRecord(
                entry[0],
                Long.parseLong(entry[1]),
                entry[2],
                entry[3],
                entry[4],
                entry[5],
                entry[6],
                entry[7],
                entry[8].split(":"),
                sampleColumns
        );
    }
}
