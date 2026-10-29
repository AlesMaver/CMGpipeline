version 1.0

# Cohort version of ViralScreen_CRAM: scatters the per-sample viral/pathogen screen over many CRAMs.

import "ViralScreen_CRAM.wdl" as vs

workflow ViralScreen_cohort {
  input {
    Array[String] sample_basenames
    Array[File] input_crams
    Array[File] input_cram_indexes
    File reference_fa
    File reference_fai
    File viral_panel_fa
    File viral_panel_fai
    File viral_panel_amb
    File viral_panel_ann
    File viral_panel_bwt
    File viral_panel_pac
    File viral_panel_sa
    Int threads = 8
  }

  scatter (i in range(length(sample_basenames))) {
    call vs.ViralScreen {
      input:
        sample_basename = sample_basenames[i],
        input_cram = input_crams[i],
        input_cram_index = input_cram_indexes[i],
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        viral_panel_fa = viral_panel_fa,
        viral_panel_fai = viral_panel_fai,
        viral_panel_amb = viral_panel_amb,
        viral_panel_ann = viral_panel_ann,
        viral_panel_bwt = viral_panel_bwt,
        viral_panel_pac = viral_panel_pac,
        viral_panel_sa = viral_panel_sa,
        threads = threads,
        min_aligned_length = 50,
        max_mismatch_rate = 0.05
    }
  }

  output {
    Array[File] viral_screen_tsv = ViralScreen.viral_screen_tsv
    Array[File] viral_summary_tsv = ViralScreen.viral_summary_tsv
    Array[File] viral_bam = ViralScreen.viral_bam
    Array[File] viral_bai = ViralScreen.viral_bai
  }
}
