(Ctrl+K V)
# Bioinformatics pipeline
The following command lines allows to generate a unique **OTU table** from `.fast5` files of eDNA metabarcoding data from **Nanopore sequencing technology**.

This repository contains : 
- `ngs_rg_nanopore_run_quality_control-master.tar.gz` : script to check the basecalling and run qualities **voir si je peux le publier ou pas**
- `decona_lucie` : je ne pense pas avoir fait de modification de ce script initial de decona donc il faut peut être mieux que je dise qu'il suffit de cloner le repository git de decona 
- `DB240625_MiFish_Actino_v2_modified.fasta` : **à publier sur Zenodo??**
- `barcode_list.txt` : donner un exemple de format de fichier ? 
- `script_reclustering_241003.sh` : est ce que je le laisse en script comme ça, ou je le met à la suite dans ce readme ? 
- les 5 scripts pythons qui sont appelés par `script_reclustering_241003.sh`.

## Basecalling with guppy 
https://github.com/asadprodhan/GPU-accelerated-guppy-basecalling
``` bash 
# To find guppy config files : 
guppy_basecaller --print_workflows

# Example for a R9.4 Flowcell :
guppy_basecaller --input_path /path/to/$experiment_name/$run_name/$run_id/fast5_pass/ --save_path $xperiement_name/$run_name/basecalled/ -c /opt/ont/guppy/data/dna_r9.4.1_450bps_sup.cfg -x "cuda:0" --compress_fastq

# Example for a R10.4 Flowcell :
guppy_basecaller --input_path  /path/to/$experiment_name/$run_name/$run_id/fast5/ --save_path $xperiement_name/$run_name/basecalled/ -c /opt/ont/guppy/data/dna_r10.3_450bps_sup.cfg -x "cuda:0" --compress_fastq
```

To check the basecalling and run qualities : 
```bash
python -m venv venv #initialize a virtual environment
source venv/bin/activate #create the virtual environment
pip install /path/to/ngs_rg_nanopore_run_quality_control-master.tar.gz #package installation
ngs_rg_nanopore_run_quality_control --help

ngs_rg_nanopore_run_quality_control --sequencing_summary /path/to/sequencing_summary_*.txt --output_directory /path/to/output/directory/

deactivate
```

## Demultiplexing with Porechop 

```bash 
porechop -i //fastq_source --require_two_barcodes -b //fastq_results --threads 16 --verbosity 2 > output 2>&1
# -b : demultiplex the reads into bins based on which barcode was found
# --require_two_barcodes : By default, Porechop only requires a single barcode match to bin a read. If you use the --require_two_barcodes option, then it will be much more stringent and assess the start and end of the read independently. I.e. to be binned, the start of a read must have a good match for a barcode and the end of the read must also have a good match for the same barcode.
# The --verbosity option will change the amount of progress info:
#   --verbosity 0 gives no progress output.
#   --verbosity 1 (default) gives summary info about end adapter trimming and middle adapter splitting.
#   --verbosity 2 shows the actual trimmed/split sequences for each read (described more below).
#   --verbosity 3 shows tons of data (mainly for debugging).
```
###### Reads Length distribution
```bash 

cd ~/fasTmp/$Experiment_name/$Run_name/Res_Porechop/
wc -l *.fastq | head -n -1 | awk '{printf "%.0f\n", $1/4}' | Rscript -e 'lines <- (readLines ("stdin"));data <- data.frame(NbReads = as.numeric(lines));dotchart(data$NbReads,xlab="Number of reads", labels = c("BC01","BC02","BC03","BC04","BC05","BC06","BC07","BC08","BC09", "BC10", "BC11","BC12", "BC13","BC15","BC23","BC24","none"), main="Number of reads in each barcode file")' ; mv Rplots.pdf DotChart_NbReadsPerBC.pdf

for BC in *.fastq ; do
  sed -n '1~4s/^@/>/p;2~4p' $BC > ${BC:0:-6}.fasta
  samtools faidx ${BC:0:-6}.fasta
  cut -f2 "${BC:0:-6}".fasta.fai | Rscript -e 'data <- as.numeric (readLines ("stdin")); summary(data); data <- data[data<1000]; hist(data, breaks = 500, xlab = "Reads length", main =  "Distribution of reads length")'
  mv Rplots.pdf ReadsLength_Distribution_"${BC:0:-6}".pdf
done
```
**examples**

``` bash
for BC in *.fasta ; do
  cut -f2 $BC.fai | awk -v OFS='\t' -v file=${BC:0:-6} '{print $0, file}' >> Size_Distrib.txt
done
```
R script 
```r
library(ggplot2)
library(ggridges)

data <- read.table("Size_Distrib.txt")
names(data) <- c("Size","BC")

ggplot(data, aes(x = Size, y = BC)) + 
  geom_density_ridges(scale = 1.5) + 
  labs(title = "Size Distribution of each sample", subtitle = "<Run_name>") + 
  xlab("Size (pb)") + ylab("Sample")
ggsave(filename = paste0(path,"SizeDistrib.svg"), width = 10, height = 10)
```
**examples**

## Decona 
https://github.com/Saskia-Oosterbroek/decona
```bash 
conda activate decona1.4
export PATH=$PATH:/path/to/minimap2

mkdir -p /path/to/$Experiment_name/$Run_name/Res_Decona/

cd /path/to/$Experiment_name/$Run_name/Res_Porechop/
for BC in BC*.fastq ; 
do 
  mkdir /path/to/$Experiment_name/$Run_name/Res_Decona/barcode${BC:2:-6}/ 
  ln -s /path/to/$Experiment_name/$Run_name/Res_Porechop/$BC /path/to/$Experiment_name/$Run_name/Res_Decona/barcode${BC:2:-6}/
done 
cd /home/eDNA/fasTmp/$Experiment_name/$Run_name/Res_Decona/

decona -l 180 -m 250 -q 10 -c 0.97 -n 10 -k 5 -i /path/to/$Experiment_name/$Run_name/Res_Decona/ -T 32 -fNCAM > /path/to/$Experiment_name/$Run_name/Res_Decona/output_97_n10 2>&1 

# Creating the barcode_list.txt file

script_reclustering_250403.sh /path/to/DB250403_MiFish_Actino_v2_modified.fasta

conda deactivate
```