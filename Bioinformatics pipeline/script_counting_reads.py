"""
Script to count the number of reads in each cluster and create a matrix otu table. 
Need : 
    - first argument : barcode_list.txt : Is a list of barcode names (defining column names of the result file) 
    - second argument : the header of the result file (except barcodes list)
    - "2nd_clust.fasta.clstr" : Is an output of cd-hit-est with the details of each clusters. 
    - "BLAST_out_reclustered_summary_tax_seq.txt" : Is an output format 6 from blast with header and sequence's additionnal column
Resulting in "BLAST_out_reclustered_summary_tax_seq_counts.txt" with header and reads counts matrix
"""

import sys

# calculate the total number of read in Dictionary
def sum_recluster_Dict(recluster_Dict):
    recluster_list = list(recluster_Dict.keys())
    for v in recluster_Dict.values(): 
        recluster_list += v
    sum=0
    for r in recluster_list: 
        #print(r)
        sum += int((r.split("-")[2]))
    print(sum)

# Creates a format string with the correct number of {} separated by tabs, based on the length of the headers list
def write_headers(filout, headers):
    filout.write('\t'.join(['{}'] * len(headers)).format(*headers) + '\n')
    #print('\t'.join(['{}'] * len(headers)).format(*headers) + '\n')

# Write a line to the output file
def write_line(filout, line, nb_reads_Dict, barcode_list):
    # Prepare the read counts for writing, defaulting to 0 if the key is not present
    read_counts = [nb_reads_Dict.get(i, 0) for i in barcode_list]
    filout.write('\t'.join(['{}'] * (1+len(read_counts))).format(line, *read_counts) + '\n')

def display_file_content(file_path):
    with open(file_path, 'r') as filin:
        content = filin.read()
        print(content)



def main():
    barcode_list = []
    with open(sys.argv[1] ,'r') as barcode_list_file:
        barcode_list = [line.strip() for line in barcode_list_file]
    #print(barcode_list)

    header = []
    #headers = ['clusters id', '%ID', 'alignment length', 'mismatches', 'gap opens', 'evalue', 'bit score', 'qcovs', 'DB id','tax id', 'Family', 'Genus', 'Species', 'sequence']
    with open(sys.argv[2] ,'r') as header_file:
        header = [line.strip() for line in header_file]
    header=header+barcode_list
    #print(headers)

    # Dictionary to hold reclustered data
    recluster_Dict = {}
    # List to hold clusters temporarily
    clusters = []

    # Open the cluster file and read its content
    # It willl create a dictionnary with the cluster representative as key and all clusters in the "bigcluster" as value
    with open("2nd_clust.fasta.clstr", "r") as file:
        for line in file:
            line = line.strip()
            #print(line)
            # If the line starts with ">", it's a new cluster
            if line.startswith(">"):
                if clusters:
                    # write the previous cluster in a dictionnary with the clusters list
                    recluster_Dict[cluster_id] = clusters
                    clusters = []
            else:
                # If the line contains "*", it indicates the line is the cluster ID
                if "*" in line:
                    cluster_id = line.split()[2][1:-3]
                    #print(cluster_id)
                    recluster_Dict[cluster_id] = []  # Initialize dictionary key with an empty list
                else:
                    # Add element to the clusters list
                    clusters.append(line.split()[2][1:-3])
        # Add the last clust to the dict
        recluster_Dict[cluster_id] = clusters

    # Display the final dictionary
    #print(recluster_Dict)

    # calculate the total number of read in Dictionary
    sum_recluster_Dict(recluster_Dict)

    # for each line (cluster) in the result file, add how many reads it get in each barcode
    with open("BLAST_out_reclustered_summary_tax_seq.txt", 'r') as filin:
        with open("BLAST_out_reclustered_summary_tax_seq_counts.txt", 'w') as filout:
            write_headers(filout, header)
            compteur = 0
            for line in filin:
                line = line.strip()
                #print(line)
                # Initialize the read counters
                nb_reads_Dict = {}
                for barcode in barcode_list :
                    nb_reads_Dict[barcode] = 0
                # Check the cluster ID and update read counters accordingly
                nb_reads_Dict[line.split("-")[0]] = int(line.split("-")[2])
                compteur += int(line.split("-")[2])
                #print("key : ", line.split("-")[0], "value: ", int(line.split("-")[2]))
                # If the the dictionary is not empty for this key, update reads counters for each cluster and each barcode
                if recluster_Dict[line.split()[0]] != []:
                    for c in recluster_Dict[line.split()[0]]:
                        #print(c)
                        nb_reads_Dict[c.split("-")[0]] += int(c.split("-")[2])
                        compteur += int(c.split("-")[2])
                #print(nb_reads_Dict)
                # Write the processed line to the output file
                write_line(filout, line, nb_reads_Dict, barcode_list)
            print(compteur)

# Execute the script
if __name__ == "__main__":
    main()
    #display_file_content("BLAST_out_reclustered_summary_tax_seq_counts.txt")
