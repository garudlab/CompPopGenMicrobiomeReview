from __future__ import division
import os
import config
import parse_midas_data
import core_gene_utils

#intermediate_filename_template = '%s.txt.gz'
patric_directory = '%spatric_db/' % (parse_midas_data.data_directory)
tree_directory = '%stree/' % (parse_midas_data.data_directory)

frn_path = "%sSSU_rRNA_all_species.frn" % (tree_directory)
frn_aligned_path = "%sSSU_rRNA_all_species_muscle.frn" % (tree_directory)
frn_aligned_clean_path = "%sSSU_rRNA_all_species_muscle_clean.frn" % (tree_directory)

outgroup_path = "%soutgroup.frn" % (tree_directory)

def get_patric_annotation_files():
    #pathway_directory = '%spatric_pathway/' % (parse_midas_data.data_directory)

    os.system("mkdir -p %s" % patric_directory)
    #os.system("mkdir -p %s" % pathway_directory)

    #intermediate_filename = intermediate_filename_template % (pairwise_directory, species_name)

    # get a list of specis to run this script on.
    good_species_list = parse_midas_data.parse_good_species_list()

    for species_name in good_species_list:

        core_genes = core_gene_utils.parse_core_genes(species_name)

        genome_ids = set([".".join(core_gene.split(".", 2)[:2]) for core_gene in core_genes])

        for genome_id in genome_ids:

            print(species_name, genome_id)


            cmnd_subsystem = "curl ftp://ftp.patricbrc.org/genomes/%s/%s.PATRIC.subsystem.tab -o %s/%s.PATRIC.subsystem.tab" % (genome_id, genome_id, patric_directory, genome_id)
            cmnd_pathway= "curl ftp://ftp.patricbrc.org/genomes/%s/%s.PATRIC.pathway.tab -o %s/%s.PATRIC.pathway.tab" % (genome_id, genome_id, patric_directory, genome_id)
            cmnd_features= "curl ftp://ftp.patricbrc.org/genomes/%s/%s.PATRIC.features.tab -o %s/%s.PATRIC.features.tab" % (genome_id, genome_id, patric_directory, genome_id)

            cmnd_cat = "cat %s/%s.PATRIC.features.tab | cut -f6,21 > %s/%s.kegg.txt" % (patric_directory, genome_id, patric_directory, genome_id)
            cmnd_bzip2 = "bzip2 -k %s/%s.kegg.txt" % (patric_directory, genome_id)

            cmnd_rRNAs = "curl ftp://ftp.patricbrc.org/genomes/%s/%s.PATRIC.frn -o %s/%s.PATRIC.frn" % (genome_id, genome_id, patric_directory, genome_id)
            #ftp://ftp.patricbrc.org/genomes/742726.3/

            #os.system(cmnd_subsystem)
            #os.system(cmnd_pathway)
            #os.system(cmnd_features)
            #os.system(cmnd_cat)
            #os.system(cmnd_bzip2)
            #os.system(cmnd_rRNAs)





class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.frn') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list



def get_16S_fasta():

    good_species_list = parse_midas_data.parse_good_species_list()


    frn = open(frn_path, 'w')

    for species_name in good_species_list:

        core_genes = core_gene_utils.parse_core_genes(species_name)

        genome_ids = list(set([".".join(core_gene.split(".", 2)[:2]) for core_gene in core_genes]))
        genome_id = genome_ids[0]

        species_name_frn_path = "%s/%s.PATRIC.frn" %  (patric_directory, genome_id)

        species_name_frn = classFASTA(species_name_frn_path).readFASTA()

        counted_rRNA = False

        for species_name_frn_name, species_name_frn_seq in species_name_frn:

            if 'ssuRNA' not in species_name_frn_name:
                continue

            if len(species_name_frn_seq) < 1200:
                continue

            if counted_rRNA == False:
                #print(species_name, species_name_frn_name, len(species_name_frn_seq))

                species_name_frn_name_split = species_name_frn_name.split()
                frn_header = species_name+'|'+species_name_frn_name_split[0].split('|')[1]

                frn.write('>%s\n' % frn_header)

                seq_split = [species_name_frn_seq[i:i+80] for i in range(0, len(species_name_frn_seq), 80)]

                for seq in seq_split:

                    frn.write('%s\n' % seq)

                frn.write('\n')

                counted_rRNA = True

            else:
                continue



    # same for outgroup

    outgroup = classFASTA(outgroup_path).readFASTA()
    print(outgroup[0][0])
    frn.write('>%s\n' % outgroup[0][0])
    seq_outgroup_split = [outgroup[0][1][i:i+80] for i in range(0, len(outgroup[0][1]), 80)]

    for seq in seq_outgroup_split:

        frn.write('%s\n' % seq)

    frn.write('\n')



    frn.close()

    os.system('muscle -in %s -out %s' % (frn_path, frn_aligned_path))



def clean_16S_fasta(max_fraction_empty=0.9):

    frn_aligned = classFASTA(frn_aligned_path).readFASTA()

    n = len(frn_aligned)

    frn_aligned_seqs = [x[1] for x in frn_aligned]
    frn_aligned_seqs_names = [x[0] for x in frn_aligned]

    frns = []

    frn_aligned_clean = open(frn_aligned_clean_path, 'w')

    for site in zip(*frn_aligned_seqs):

        fraction_empty = site.count('-')/n

        if fraction_empty > max_fraction_empty:
            continue

        frns.append(site)


    clean_sites_list = zip(*frns)

    for clean_sites_idx, clean_sites in enumerate(clean_sites_list):
        clean_sites_species = frn_aligned_seqs_names[clean_sites_idx]
        clean_sites_seq = "".join(clean_sites)

        frn_aligned_clean.write('>%s\n' % clean_sites_species)

        clean_sites_seq_split = [clean_sites_seq[i:i+80] for i in range(0, len(clean_sites_seq), 80)]

        for seq in clean_sites_seq_split:

            frn_aligned_clean.write('%s\n' % seq)

        frn_aligned_clean.write('\n')


    frn_aligned_clean.close()






#get_patric_annotation_files()
# raxml/8.2.4
#get_16S_fasta()
#clean_16S_fasta()
