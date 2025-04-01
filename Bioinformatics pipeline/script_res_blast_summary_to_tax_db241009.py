""" 
Script to retrieve taxonomic lineage (Familly Genus, species) from one (or multiple) tax ID.
Need : 
	- /home/edna/src/rankedlineage_tabRemoved.dmp : see how o dowload it at the begining of this script 
	- BLAST_out_reclustered_best_res.txt : output format 6 from blast 
Resulting in BLAST_out_reclustered_summary_tax.txt
"""

# je supprime cette premiere partie car deja telecharge rankedlineage_tabRemoved.dmp

# define commands for downloading taxonomy database from NCBI and saving it in a folder with the current date
# if folder exists, it overwrites it
# etLineage = '''
# ="$(date +"%d-%m-%Y")"
# mkdir ${d}_taxdump
# cd ${d}_taxdump
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
# unzip new_taxdump.zip
# mv rankedlineage.dmp ..
# cd ..
# tr -d "\t" < rankedlineage.dmp > rankedlineage_tabRemoved.dmp
# rm -f ${d}_taxdump/*
# rm -f rankedlineage.dmp
# '''
# create a child bash process
# process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE,
#                          stdout=subprocess.PIPE, universal_newlines=True)
# execute getLineage command and print output on the screen
# out, err = process.communicate(getLineage)
# out

# def taxonomy_dictionary(ncbi_taxonomy):

taxonomy_Dict = {}

# This will parse the ncbi taxonomy file to get the column of interest
# with open("/env/cns/proj/projet_CZK/Database_MiFish/rankedlineage_tabRemoved.dmp", "r") as file:
with open("/home/edna/src/rankedlineage_tabRemoved.dmp", "r") as file:
	for line in file:
		l3 = line.strip().split('|')
		tax_id = l3[0]
		domain = l3[9]
		phylum = l3[7]
		classif = l3[6]
		order = l3[5]
		family = l3[4]
		genus = l3[3]
		species = l3[1]
		if "tax_id" not in line:
			taxonomy_Dict[tax_id] = [domain, phylum,classif, order, family, genus, species]

with open("BLAST_out_reclustered_best_res.txt",'r') as filin:
	with open("BLAST_out_reclustered_summary_tax.txt",'w+') as filout:
		for line in filin:
			line = line.rstrip('\n')
			if line.strip() != '':
				taxon = line.split()[-1].split("_")[0]
				family = taxonomy_Dict[taxon][-3]
				genus = taxonomy_Dict[taxon][-2]
				species = taxonomy_Dict[taxon][-1].split()[1]
				#print(species)
				for taxon in line.split()[-1].split("_") :
					#print(taxon)
					# if the %identity between the sequence and the reference is less thant 88.2, we can't know. 
					if float(line.split()[1])<88.2 :
						genus = ""
						species = ""
						family = ""
					else : 
						# if the %identity between the sequence and the reference is between 88.2 and 96.4, we assign at Family level only. 
						if float(line.split()[1])<96.4 :
							species = ""
							genus = ""
							if taxonomy_Dict[taxon][-3] != family and taxonomy_Dict[taxon][-3] not in family.split("-") : 
								family = family + "-" + taxonomy_Dict[taxon][-3]
						# if the % of identity between the sequence and the reference is between 96.4 and 97, we print Family, Genus name. 
						else : 
							if float(line.split()[1])<97 :
								species = ""
								if taxonomy_Dict[taxon][-2] != genus and taxonomy_Dict[taxon][-2] not in genus.split("-") : 
									genus = genus + "-" + taxonomy_Dict[taxon][-2]
								if taxonomy_Dict[taxon][-3] != family and taxonomy_Dict[taxon][-3] not in family.split("-") : 
									family = family + "-" + taxonomy_Dict[taxon][-3]
							else : # if %ID > 97, we keep print Family, Genus and species names
								if taxonomy_Dict[taxon][-3] != family and taxonomy_Dict[taxon][-3] not in family.split("-") : 
									family = family + "-" + taxonomy_Dict[taxon][-3]
								if taxonomy_Dict[taxon][-2] != genus and taxonomy_Dict[taxon][-2] not in genus.split("-") : 
									genus = genus + "-" + taxonomy_Dict[taxon][-2]
								if taxonomy_Dict[taxon][-1].split()[1] != species and taxonomy_Dict[taxon][-1].split()[1] not in species.split("-") : 
									#print(species)
									species = species + "-" + taxonomy_Dict[taxon][-1].split()[1]
				filout.write('{}\t{}\t{}\t{}\n'.format(line,family,genus,species))
