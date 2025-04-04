#!/bin/bash

# Navigate to the root of the three runs
mkdir Res_Decona_all_runs/
for run in * ; do 
  if [ -e "$run/Res_Decona" ] ; then
    (
    cd $run/Res_Decona/
      # Concatenate the clusters from Medaka into an all_medaka_fasta file for each barcode, and add the attribute barcodeXX to each header
      (
      cd data/
        mkdir -p ../result/Medaka
        for folder in * ; do
          if [ -e "$folder/multi-seq" ] ; then
            (
            cd "$folder/multi-seq" || exit
              for con in consensus_medaka_*.fasta ; do
                if [ -e "$con" ] ; then
                  #echo "$con"
                  cat "$con" | awk '{print $1}' > n-"$con"
                  sed -i '1 s/.*/&_'"$con"'/' n-"$con"
                  sed -i 's/>.*consensus_medaka_/>consensus_medaka-/' n-"$con"
                else
                  continue
                fi
              done
              if [ "$(ls -1 n-con*.fasta 2>/dev/null | wc -l)" -gt 0 ] ; then
                rm all_medaka_fastas.fasta
                # concatenating clusters
                cat n-con*.fasta >> all_medaka_fastas.fasta
              fi
              # Add the attribute
              awk -v barcode="$folder" '{ if ($0 ~ /^>/) print ">" barcode "-" substr($0, 2); else print $0 }' all_medaka_fastas.fasta > "all_medaka_fastas_${folder}.fasta"
              # Move the fasta files into a result folder
              cp "all_medaka_fastas_${folder}.fasta" ../../../result/Medaka/"all_medaka_fastas_${folder}.fasta" 
            ) # Exit multi-seq, go to data
          fi
        done
      ) # Exit data, go to $run/Res_Decona/

      # Concatenate the barcodes and add the run attribute
      (
      cd result/Medaka/
        rm All_medaka_fastas*.fasta
        cat all_medaka_fastas_*.fasta >> All_medaka_fastas_all_barcodes.fasta
        awk -v run="$run" '{ if ($0 ~ /^>/) print ">" run "_" substr($0, 2) ; else print $0 }' All_medaka_fastas_all_barcodes.fasta > "All_medaka_fastas_all_barcodes_$run.fasta"
      ) # Exit result/Medaka/, go to $run/Res_Decona
    )
  fi
  cp $run/Res_Decona/result/Medaka/All_medaka_fastas_all_barcodes_$run.fasta Res_Decona_all_runs/ 
done 

cd Res_Decona_all_runs/
rm All_medaka_fastas_all_barcodes_all_runs.fasta
cat All_medaka_fastas_all_barcodes_*.fasta >> All_medaka_fastas_all_barcodes_all_runs.fasta

# Reclustering all barcodes together
cd-hit-est -i All_medaka_fastas_all_barcodes_all_runs.fasta -o 2nd_clust.fasta -c 0.95 -n 5 -d 0 -aS 0.9 -G 0 -M 0 -T 32 -g 1 > output_reclustering 2>&1 

# Taxonomic assignation 
makeblastdb -in $1 -dbtype nucl -parse_seqids
blastn -query 2nd_clust.fasta -db $1 -perc_identity 80 -outfmt "6 qseqid pident length mismatch gapopen evalue bitscore salltitles sallseqid" -max_target_seqs 20 -max_hsps 500 -num_threads 32 > BLAST_out_reclustered.txt ;
# blastn options #

# query: Name of the file containing the query sequence(s), or ‘-‘ if these are provided on standard input.

# db: File name of BLAST database to search the query against.
  # Unless an absolute path is used, the database will be searched relative to the current working directory first,
  # then relative to the value specified by the BLASTDB environment variable,
  # then relative to the BLASTDB configuration value specified in the configuration file.
#
# perc_identity: Minimum percent identity of matches to report

# outfmt: Allows for the specification of the search application’s output format.
  # A listing of the possible format types is available via the search application’s -help option.
  # If a custom output format is desired, this can be specified by providing a quoted string composed of the desired output format (tabular, tabular with comments, or comma-separated value), a space, and a space delimited list of output specifiers.
  # The list of supported output specifiers is available via the -help command line option. Unsupported output specifiers will be ignored.
  # This should be specified using double quotes if there are spaces in the output format specification (e.g.: outfmt "7 sseqid ssac qstart qend sstart send qseq evalue bitscore").
#
# max_target_seqs: Maximum number of aligned sequences to keep from the BLAST database.

# max_hsps: Maximum number of HSPs (alignments) to keep for any single query-subject pair. The HSPs shown will be the best as judged by expect value.
  # This number should be an integer that is one or greater.
  # If this option is not set, BLAST shows all HSPs meeting the expect value criteria. Setting it to one will show only the best HSP for every query-subject pair
#



# Concatenation match with identical bitscore
python script_concatenating_double_OTUs.py BLAST_out_reclustered.txt

# Adding taxonomy level to the blast result
python script_res_blast_summary_to_tax.py
python script_res_blast_summary_to_tax_without_filter.py


# Identifying and adding clusters which do not has a taxonomic assignation at all
grep ">" 2nd_clust.fasta > elements.txt
grep ">" 2nd_clust.fasta | sed 's/>//' > elements_clean.txt
rm not_found.txt
while read -r element; do
  if ! grep -q "$element" BLAST_out_reclustered_summary_tax.txt; then
    echo "$element" >> not_found.txt
  fi
done < elements_clean.txt
# re-BLASTing them (optionnal)
#rm queries.fasta
#while read -r element; do
#  grep -A 1 "$element" 2nd_clust.fasta >> queries.fasta
#done < not_found.txt
#
#blastn -query queries.fasta -db "$BLASTdir" -perc_identity 80 -outfmt "6 qseqid pident length mismatch gapopen evalue bitscore salltitles sallseqid" -max_target_seqs 20 -max_hsps 500 -num_threads 32 > result.blast ;
wc -l not_found.txt
cat not_found.txt >> BLAST_out_reclustered_summary_tax.txt

# Adding consensus sequence to the summary fie. 
python script_adding_seq_to_res_sum_tax.py

# Couting reads in each cluster to make one OTU table pour all samples 
python script_counting_reads.py

