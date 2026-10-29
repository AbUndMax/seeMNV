package com.github.abundmax.seemnv.core;

import com.github.abundmax.seemnv.core.fieldArray.FieldConfig;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;

public class testCodonToMNV {

   @BeforeAll
   static void setup() {
      Core.fieldToConfig.put("AF", FieldConfig.buildFromHeaderString("##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">"));
      Core.fieldToConfig.put("DP", FieldConfig.buildFromHeaderString("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">"));
      Core.fieldToConfig.put("GQ", FieldConfig.buildFromHeaderString("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
      // Core.fieldToConfig.put("GT", FieldConfig.buildFromHeaderString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
      Core.fieldToConfig.put("PL", FieldConfig.buildFromHeaderString("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"));
      Core.fieldToConfig.put("AD", FieldConfig.buildFromHeaderString("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"));
   }

   @Test
   void testCalculateMNVEntriesInCodon() {

      VCFRecord entry31a = VCFRecord.fromString("NC_021508.1\t31\tdnaA_c10_b1_alt1\tT\tA\t.\tPASS\tGENE=dnaA;CODON=10;CODON_BASE=1\tGT:AD:AF:DP:GQ:PL\t0/1:30,20:0.400:50:99:80,0,100\t0/0:50,0:0.000:50:99:0,60,600\t0/0:50,0:0.000:50:99:0,60,600");
      VCFRecord entry31b = VCFRecord.fromString("NC_021508.1\t31\tdnaA_c10_b1_alt2\tT\tC\t.\tPASS\tGENE=dnaA;CODON=10;CODON_BASE=1\tGT:AD:AF:DP:GQ:PL\t0/0:50,0:0.000:50:99:0,60,600\t0/1:30,20:0.400:50:99:80,0,100\t0/0:50,0:0.000:50:99:0,60,600");
      VCFRecord entry32 =  VCFRecord.fromString("NC_021508.1\t32\tdnaA_c10_b2_alt1\tG\tA\t.\tPASS\tGENE=dnaA;CODON=10;CODON_BASE=2\tGT:AD:AF:DP:GQ:PL\t1/1:30,20:0.400:50:99:80,0,100\t0/0:50,0:0.000:50:99:0,60,600\t0/1:30,20:0.400:50:99:80,0,100");
      VCFRecord entry33 =  VCFRecord.fromString("NC_021508.1\t33\tdnaA_c10_b3_alt1\tG\tA\t.\tPASS\tGENE=dnaA;CODON=10;CODON_BASE=3\tGT:AD:AF:DP:GQ:PL\t1/1:30,20:0.400:50:99:80,0,100\t0/0:50,0:0.000:50:99:0,60,600\t0/0:50,0:0.000:50:99:0,60,600");

      Codon cdn = new Codon(9, "TGG");

      cdn.addRecordAtCodonPos(0, entry31a);
      cdn.addRecordAtCodonPos(0, entry31b);
      cdn.addRecordAtCodonPos(1, entry32);
      cdn.addRecordAtCodonPos(2, entry33);

      List<VCFRecord> mnvEntries = Core.calculateMNVEntriesInCodon(cdn);

      VCFRecord cmnv1 = VCFRecord.fromString("NC_021508.1\t31\tdnaA_c10_b1_alt1\tTGG\tAAA\t.\t.\tcMNV=3|9|TGG|AAA|[., ., .]|[PASS, PASS, PASS]|[31, 32, 33]|[[0, 1], [1, 2], [1, 1]]|[0, 1]\tGT:AD:AF:DP:GQ:PL\t0|1:30,20:0.4:50:99:80,0,100\t0|0:50,0:0.0:50:99:0,60,600\t0|0:30,0:0.0:50:99:80,60,600");
      VCFRecord cmnv2 = VCFRecord.fromString("NC_021508.1\t32\tdnaA_c10_b2_alt1\tGG\tAA\t.\t.\tcMNV=2|9|TGG|TAA|[/, ., .]|[/, PASS, PASS]|[/, 32, 33]|[/, [1, 2], [1, 1]]|[1, 0]\tGT:AD:AF:DP:GQ:PL\t1|0:30,20:0.4:50:99:80,0,100\t0|0:50,0:0.0:50:99:0,60,600\t0|0:30,0:0.0:50:99:80,60,600");

      assertEquals(cmnv1.toString(), mnvEntries.get(0).toString());
      assertEquals(cmnv2.toString(), mnvEntries.get(1).toString());
   }

   @Test
   void createsTripleThenResidualPairForMixedFourSampleCalls() {
      Codon codon = new Codon(42, "TGG");

      codon.addRecordAtCodonPos(0, record(
              "p1", 100, "T", "A", "11", "q10",
              "0/1", "0/1", "0/0", "1/1"
      ));
      codon.addRecordAtCodonPos(1, record(
              "p2", 101, "G", "C", "12", "PASS",
              "0/1", "0/0", "0/1", "1/1"
      ));
      codon.addRecordAtCodonPos(2, record(
              "p3", 102, "G", "T", "13", "LowQual",
              "0/1", "0/0", "0/1", "1/1"
      ));

      List<VCFRecord> mnvs = Core.calculateMNVEntriesInCodon(codon);

      assertEquals(2, mnvs.size());

      assertMnv(
              mnvs.get(0),
              100,
              "p1",
              "TGG",
              "ACT",
              "cMNV=3|42|TGG|ACT|[11, 12, 13]|[q10, PASS, LowQual]|" +
                      "[100, 101, 102]|[[1, 3], [1, 3], [1, 3]]|[1, 2]",
              "0|1", "0|0", "0|0", "1|1"
      );

      // Sample 3 participates in positions 2 and 3, but not position 1.
      // The triple's samples are removed from this pair's participation mask.
      assertMnv(
              mnvs.get(1),
              101,
              "p2",
              "GG",
              "CT",
              "cMNV=2|42|TGG|TCT|[/, 12, 13]|[/, PASS, LowQual]|" +
                      "[/, 101, 102]|[/, [1, 3], [1, 3]]|[0, 1]",
              "0|0", "0|0", "0|1", "0|0"
      );
   }

   @Test
   void createsSeparatePairsForDifferentAlternateRecordsAtOneCodonPosition() {
      Codon codon = new Codon(7, "AAC");

      codon.addRecordAtCodonPos(0, record(
              "base1-G", 200, "A", "G", "17", "PASS",
              "0/1", "0/0", "0/1"
      ));
      codon.addRecordAtCodonPos(0, record(
              "base1-T", 200, "A", "T", "18", "PASS",
              "0/0", "0/1", "0/1"
      ));
      codon.addRecordAtCodonPos(2, record(
              "base3-A", 202, "C", "A", "19", "PASS",
              "0/1", "0/1", "0/1"
      ));

      List<VCFRecord> mnvs = Core.calculateMNVEntriesInCodon(codon);

      assertEquals(2, mnvs.size());

      assertMnv(
              mnvs.get(0),
              200,
              "base1-G",
              "AAC",
              "GAA",
              "cMNV=2|7|AAC|GAA|[17, /, 19]|[PASS, /, PASS]|" +
                      "[200, /, 202]|[[0, 2], /, [0, 3]]|[0, 2]",
              "0|1", "0|0", "0|1"
      );

      assertMnv(
              mnvs.get(1),
              200,
              "base1-T",
              "AAC",
              "TAA",
              "cMNV=2|7|AAC|TAA|[18, /, 19]|[PASS, /, PASS]|" +
                      "[200, /, 202]|[[0, 2], /, [0, 3]]|[0, 2]",
              "0|0", "0|1", "0|1"
      );
   }

   @Test
   void createsTwoBaseMnvWhenOnlyTheFirstTwoCodonPositionsHaveRecords() {
      Codon codon = new Codon(0, "TGA");

      codon.addRecordAtCodonPos(0, record(
              "base1-C", 300, "T", "C", ".", "PASS",
              "1/1"
      ));
      codon.addRecordAtCodonPos(1, record(
              "base2-T", 301, "G", "T", "99", "q20",
              "1/1"
      ));

      List<VCFRecord> mnvs = Core.calculateMNVEntriesInCodon(codon);

      assertEquals(1, mnvs.size());

      assertMnv(
              mnvs.getFirst(),
              300,
              "base1-C",
              "TG",
              "CT",
              "cMNV=2|0|TGA|CTA|[., 99, /]|[PASS, q20, /]|" +
                      "[300, 301, /]|[[1, 1], [1, 1], /]|[1, 1]",
              "1|1"
      );
   }

   private VCFRecord record(
           String id,
           long pos,
           String ref,
           String alt,
           String qual,
           String filter,
           String... genotypes
   ) {
      String[] samples = new String[genotypes.length];

      for (int i = 0; i < genotypes.length; i++) {
         samples[i] = genotypes[i] + ":30,20";
      }

      return new VCFRecord(
              "chr1",
              pos,
              id,
              ref,
              alt,
              qual,
              filter,
              ".",
              new String[]{"GT", "AD"},
              samples
      );
   }

   private void assertMnv(
           VCFRecord actual,
           long expectedPos,
           String expectedId,
           String expectedRef,
           String expectedAlt,
           String expectedInfo,
           String... expectedGenotypes
   ) {
      assertEquals("chr1", actual.chrom());
      assertEquals(expectedPos, actual.pos());
      assertEquals(expectedId, actual.id());
      assertEquals(expectedRef, actual.ref());
      assertEquals(expectedAlt, actual.alt());
      assertEquals(".", actual.qual());
      assertEquals(".", actual.filter());
      assertEquals(expectedInfo, actual.info());

      String[] actualGenotypes = new String[actual.numSamples()];

      for (int i = 0; i < actual.numSamples(); i++) {
         actualGenotypes[i] = actual.samples()[i].split(":")[0];
      }

      assertArrayEquals(expectedGenotypes, actualGenotypes);
   }
}
