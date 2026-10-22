package com.github.abundmax.seemnv.core.fieldArray;

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
}
