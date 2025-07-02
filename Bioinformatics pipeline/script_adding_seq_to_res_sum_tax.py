"""
Script to add corresponding consensus sequence to each OTU.
Need : 
    - 2nd_clust.fasta : output from cd-hit-est, all cluster ID wit their consensus sequence 
    - BLAST_out_reclustered_summary_tax.txt : output format 6 from blast.
Resulting in BLAST_out_reclustered_summary_tax_seq.txt
"""

# Create a table with OTUs are row lines.
def main() :
    # Read the fasta file to get all sequence associated with their id
    with open("2nd_clust.fasta",'r') as file:
        seq_dict = {}
        for line in file:
            if line[0] == ">":
                seq_id = line[0:-1]
                seq_dict[line[0:-1]] = ""
            else:
                seq_dict[seq_id] = line[0:-2]

    #print(seq_dict)

    with open("BLAST_out_reclustered_summary_tax.txt", 'r') as file_in:
        with open("BLAST_out_reclustered_summary_tax_seq.txt", "w+") as file_out:
            for line in file_in:
                line = line.strip()
                line = line.split()
                new_line = ["","","","","","","","","","","","","",""]
                new_line[0:len(line)] = line
                # Search the sequence by their id in the dictionary. 
                seq_id = '>' + line[0]
                if seq_id in seq_dict :
                    new_line[-1] = seq_dict[seq_id]
                else : 
                    print("cluster_id not found in the fasta file")
                file_out.write('{}\n'.format('\t'.join(new_line)))
                #print('{}\n'.format('\t'.join(new_line)))

if __name__ == "__main__":
    #filin = sys.argv[1]
    main()



