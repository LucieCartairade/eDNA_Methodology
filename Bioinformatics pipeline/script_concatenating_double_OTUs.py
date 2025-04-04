"""
Script to concatenate the reference database id between each OTUs that have the same biscore.
Need a file which is the result format 6 of blast. Input this file name in argument
The first column is an uniq name of a cluster, and the 7th is the bitscore, 8th a tax ID.
Result in a new file : 'BLAST_out_reclustered_best_res.txt' with uniq cluster_id and possible muliple tax_id separated by '_'
"""

import sys

def main(filin) : 
    with open(filin, 'r') as filin:
        with open('BLAST_out_reclustered_best_res.txt', 'w') as filout:
            previous_line = "0\t0\t0\t0\t0\t0\t0\t1" # define a random line for the first iteration 
            tax_id = "tax_id" # same
            first = True
            for line in filin:
                line = line.strip()
                # if the OTUs id is the same than the previous line
                if (line.split()[0] == previous_line.split()[0]): # check if the cluster_id is different to the revious one
                    if line.split()[6] == previous_line.split()[6] : # if bitscore are the same
                        tax_id = previous_line.split()[-1] +  '_' + line.split()[-1] # adding new tax id to the previous tax_id
                else : # if the OTUs is a new one :
                    if first != True : # We have not to consider the fisrt linebecause we print the previous line.
                        filout.write('{}\t{}\n'.format('\t'.join(previous_line.split()[:-1]),tax_id))
                        #print('{}\t{}\t{}\n'.format('\t'.join(previous_line.split()[:-2]),tax_id, previous_line.split()[8]))
                    tax_id = line.split()[-1] # define a new tax_id
                    previous_line = line #set the line as the previous one for the next loop
                first = False
            # write the last cluster
            filout.write('{}\t{}\n'.format('\t'.join(previous_line.split()[:-1]),tax_id))


if __name__ == "__main__":
    filin = sys.argv[1]
    main(filin)














