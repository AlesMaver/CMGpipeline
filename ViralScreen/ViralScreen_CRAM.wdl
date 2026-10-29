version 1.0

# Screen an hg19 exome CRAM for reads from the capture's viral/pathogen spike-in panel
# (HPV types, EBV, CMV, MCV, HSV1/2, VZV, Parvovirus B19, T. pallidum, T. gondii loci).
# Reads unmapped to hg19 (or with an unmapped mate) are realigned to the panel with bwa mem,
# deduplicated and filtered; per-target read counts, RPM, breadth and a call are reported.

workflow ViralScreen_CRAM {
  input {
    String sample_basename
    File input_cram
    File input_cram_index
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
    Int min_aligned_length = 50
    Float max_mismatch_rate = 0.05
  }

  call ViralScreen {
    input:
      sample_basename = sample_basename,
      input_cram = input_cram,
      input_cram_index = input_cram_index,
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
      min_aligned_length = min_aligned_length,
      max_mismatch_rate = max_mismatch_rate
  }

  output {
    File viral_screen_tsv = ViralScreen.viral_screen_tsv
    File viral_summary_tsv = ViralScreen.viral_summary_tsv
    File viral_bam = ViralScreen.viral_bam
    File viral_bai = ViralScreen.viral_bai
  }
}

task ViralScreen {
  input {
    String sample_basename
    File input_cram
    File input_cram_index
    File reference_fa
    File reference_fai
    File viral_panel_fa
    File viral_panel_fai
    File viral_panel_amb
    File viral_panel_ann
    File viral_panel_bwt
    File viral_panel_pac
    File viral_panel_sa
    Int threads
    Int min_aligned_length
    Float max_mismatch_rate
  }

  command <<<
    set -euo pipefail
    T=~{threads}
    S=~{sample_basename}

    # bwa and samtools expect index files next to their fasta
    mkdir -p panel ref
    for f in ~{viral_panel_fa} ~{viral_panel_fai} ~{viral_panel_amb} ~{viral_panel_ann} ~{viral_panel_bwt} ~{viral_panel_pac} ~{viral_panel_sa}; do ln -s $f panel/; done
    ln -s ~{reference_fa} ref/; ln -s ~{reference_fai} ref/
    PANEL=panel/$(basename ~{viral_panel_fa})
    REF=ref/$(basename ~{reference_fa})

    # 1. total primary reads (denominator) + reads with read or mate unmapped
    samtools view -@ $T -c -F 0x900 -T $REF ~{input_cram} > total_reads.txt &
    samtools view -@ $T -T $REF -u -F 0x900 -e 'flag.unmap || flag.munmap' ~{input_cram} \
      | samtools collate -@ $T -O -u - \
      | samtools fastq -@ $T -n -1 r1.fq.gz -2 r2.fq.gz -s single.fq.gz -0 /dev/null - 2> /dev/null
    wait

    # 2. align to the panel, keep primary mapped, mark PCR duplicates
    { bwa mem -t $T -Y -R "@RG\tID:$S\tSM:$S" $PANEL r1.fq.gz r2.fq.gz
      bwa mem -t $T -Y -R "@RG\tID:$S\tSM:$S" $PANEL single.fq.gz | { grep -v '^@' || true; }; } 2> bwa.log \
      | samtools sort -n -@ $T -u - | samtools fixmate -m -u - - \
      | samtools sort -@ $T -u - | samtools markdup -@ $T -u - - \
      | samtools view -@ $T -b -F 0x904 -o all.bam -
    samtools index all.bam

    # 3. filter: aligned length, mismatch rate, not duplicate; drop low-complexity reads
    samtools view -h -F 0x400 all.bam | awk -v minlen=~{min_aligned_length} -v maxnm=~{max_mismatch_rate} 'BEGIN{OFS="\t"}
      /^@/ {print; next}
      { alen=0; c=$6; while (match(c,/[0-9]+[M=X]/)) { alen+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH) }
        nm=0; for(i=12;i<=NF;i++) if ($i ~ /^NM:i:/) nm=substr($i,6)
        if (alen < minlen || nm/alen > maxnm) next
        s=toupper($10); n=length(s); a=gsub(/A/,"A",s); cc=gsub(/C/,"C",s); g=gsub(/G/,"G",s); t=gsub(/T/,"T",s)
        if ((a+t)/n > 0.9 || (g+cc)/n > 0.9 || (a+g)/n > 0.95 || (cc+t)/n > 0.95) next
        print }' | samtools view -b -o $S.viral.bam -
    samtools index $S.viral.bam

    # 4. per-target metrics; DETECTED needs >=10 MAPQ20 reads spread over >=1000 covered bases (guards against single-locus pileups)
    TOTAL=$(cat total_reads.txt)
    samtools coverage $S.viral.bam > cov.tsv
    for c in $(cut -f1 $PANEL.fai); do echo -e "$c\t$(samtools view -c -q 20 $S.viral.bam "$c")"; done > mapq20.tsv
    awk -v tot=$TOTAL -v s=$S 'BEGIN{FS=OFS="\t"; print "sample","target","accession","length","reads","reads_mapq20","RPM","covered_bases","breadth_pct","mean_depth","call"}
      NR==FNR {q[$1]=$2; next}
      FNR==1 {next}
      { split($1,a,"|"); r=$4; br=$6
        call = q[$1]>=3 ? "LOW" : "-"
        if (q[$1]>=10 && $5>=1000) call="DETECTED"
        print s,a[1],a[2],$3,r,q[$1],sprintf("%.3f",r/tot*1e6),$5,sprintf("%.2f",br),sprintf("%.3f",$7),call }' mapq20.tsv cov.tsv > $S.viral_screen.tsv
    printf "sample\ttotal_reads\tpanel_reads\tdetected\n%s\t%s\t%s\t%s\n" $S $TOTAL $(samtools view -c $S.viral.bam) \
      "$(awk -F'\t' '$11=="DETECTED"{printf "%s%s",sep,$2; sep=","}' $S.viral_screen.tsv)" > $S.viral_summary.tsv
    rm -f r1.fq.gz r2.fq.gz single.fq.gz all.bam all.bam.bai
  >>>

  runtime {
    docker: "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"
    cpu: threads
    requested_memory_mb_per_core: 1000
  }

  output {
    File viral_screen_tsv = "~{sample_basename}.viral_screen.tsv"
    File viral_summary_tsv = "~{sample_basename}.viral_summary.tsv"
    File viral_bam = "~{sample_basename}.viral.bam"
    File viral_bai = "~{sample_basename}.viral.bam.bai"
  }
}
