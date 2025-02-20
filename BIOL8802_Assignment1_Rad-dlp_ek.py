import argparse # Take in arguments
import pyfaidx # Read fasta files
import pdb # Debug
#import vcf # Read in vcf files ##see below for in-terminal call
from Bio.Restriction.Restriction_Dictionary import rest_dict # Library of restriction enzymes so users can enter enzyme names instead of motifs

#in terminal
#conda install -c bioconda pyvcf3

def my_parse_args():
    parser = argparse.ArgumentParser(description='this script performs either SingleRad or ddRad sequencing on DNA from an input FASTA file and compared the results to variants from an input VCF file')
    #TODO: Complete this function as described in the docstring for GenomeFile, VCFFile, Mode, RE1, RE2
    parser.add_argument('GenomeFile', type=str, help='the path to the .fasta file containing the genome')
    parser.add_argument('VCFFile', type=str, help='the path to the .vcf file with the DNA polymorphisms')
    parser.add_argument('Mode', type=str, choices=['SingleRad', 'ddRad'], help='whether to run SingleRad or ddRad')
    parser.add_argument('RE1', type=str, help='the name (not the motif) of the restriction enzyme')
    parser.add_argument('-re2', '--RE2', type=str, help='the name (not the motif) of the second restriction enzyme (only used in ddRad mode)')
    parsed_args = parser.parse_args()
    return parsed_args

def read_fasta(path_to_fasta: str, chromosome:str) -> str:
    #TODO: complete this function as described in the docstring. You can remove 'dna=None" statement it's just a placeholder
    genome = pyfaidx.Fasta(path_to_fasta)
    dna = genome.records[chromosome][:].seq
    return dna

def find_motifs(dna:str, motif:str) -> list[int]:
    positions = []
    #TODO: complete this function as described in the docstring.
    start = 0
    while True:
        start = dna.find(motif, start)
        if start ==1:
            break
        positions.append(start)
        start += len(motif)
    return positions

def run_single_rad(dna:str, re1:str) -> list[tuple[int, int]]:
    seq_length = 100
    #TODO: use an assert statement to verify that the RE1 site is in teh rest_dict provided by biopython
    assert re1 in rest_dict
    #TODO: use the rest_dict to look up the motif associated with teh enzyme
    motif1 = rest_dict[re1]['site']
    #TODO: use the find_motifs function to find the cut sites
    re_sites = find_motifs(dna,motif1)
    #TODO: use list comprehension (preferably) or a loop (acceptable) to create the desired list of sequence sites
    sequenced_sites = [(x-seq_length, x+seq_length) for x in re_sites]
    return sequenced_sites

def run_ddrad(dna:str, re1:str, re2:str):
    min_size = 300
    max_size = 700
    seq_length = 100
    #TODO: Use an assert statement to verify that RE1 and RE2 are both in the rest_dict provided by biopython
    assert re1 in rest_dict
    assert re2 in rest_dict
    #TODO: use the rest_dict to look up the motifs associated with the two enzyme
    motif1 = rest_dict[re1]['site']
    motif2 = rest_dict[re2]['site']
    #TODO: use your find_motifs function to find the sequencing sites for RE1 and RE2
    re1_sites = find_motifs(dna, motif1)
    re2_sites = find_motifs(dna, motif2)

    #TODO: Complete the below loop, which is currently incomplete but has hints as comments
    sequenced_sites = []
    for r1 in re1_sites: #Loop through each RE1 site
        #TODO: Look left of the current site
        left_hits_r1 = [x for x in re1_sites if x < r1] #Replace None with a list comprehension to find all re1 sites left of the current site
        left_hits_r2 = [x for x in re2_sites if x < r1] #Replace None with a list comprehension to find all re2 sites left of the current site

        #TODO: create a series of nested if statements that check:
        #1) that left_hits_r2 is not empty
        if left_hits_r2 != []:
        #2) that re2 site is between min and max size
            if r1 - left_hits_r2[-1] > min_size and r1-left_hits_r2[-1] < max_size:
        #3) that re2 site is closer re1 site
                if left_hits_r2[0] < left_hits_r1[0]:
        #TODO: if all statements pass, add the two resulting sequence sites to the sequenced_sites list
                    sequenced_sites.append((r1, r1+seq_length))
                    sequenced_sites.append((left_hits_r1[0]-seq_length, left_hits_r1[0]))
                    return sequenced_sites
        #TODO: Now do the same thing but looking right
        right_hits_r1 = [x for x in re1_sites if x > r1] #Replace None with a list of comprehension to find all re1 sites right of the current site
        right_hits_r2 = [x for x in re2_sites if x > r1] #Replace None with a list comprehension to find all re2 sites right of the current site

        if right_hits_r2 != []:
            if right_hits_r2[0] -r1 > min_size and right_hits_r2[0] -r1 < max_size:
                if right_hits_r2[0] < right_hits_r1[0]:
                    sequenced_sites.append((r1, r1+seq_length))
                    sequenced_sites.append((right_hits_r1[0]-seq_length, right_hits_r1[0]))
                    return sequenced_sites

def find_variable_sites(vcf_file_path:str, sequenced_sites:list[tuple[int,int]]):
    #TODO: read in the VCF file using the VCFReader class
    vcf_obj = vcf.VCFReader(filename=vcf_file_path)
    #TODO: complete the below loop
    sequenced_sites_variable = [] #This list will keep track of all the sequenced DNA pieces that also contain variation
    for record in vcf_obj: #Create loop to read through every record
        if record.CHROM != 'NC_036780.1': #These two lines will make sure that you are only looking at record from Chromosome 1
            break
        #TODO:add an if statement here to ensure that both samples have called genotypes (num_called attribute is your friend)
        if record.num_called != 2:
            continue
        #TODO:add another if statement to ensure that one genotype is 0/0 and the other is 1/1 (use gt_nums)
        if (record.samples[0].gt_nums == '0/0' and record.samples[1].gt_nums == '1/1') or (record.samples[1].gt_nums == '0/0' and record.samples[0].gt_nums == '1/1'):
        #TODO:if both conditions pass, Use list comprehension to determine if this variants fall within sequenced DNA
            matches = [x for x in sequenced_sites if x[0] < record.POS and x[1] > record.POS]
        #TODO: if the variant falls with sequenced DNA, add it to the sequenced_sites_variables list
            if len(matches) > 0:
                sequenced_sites_variable.extend(matches)
        return sequenced_sites_variable

# ###TESTING CODE###
# # This section of code is provided to help you test your functions as you develop them. Feel free to add additional
# # debugging code to this section, but ensure that all code required to run the script remains in the function
# # definitions above. As you progress, leave the completed sections uncommented, as some variables are reused. Note
# # also that, in much of this testing code, your command-line input (other than file locations) is ignored, but you
# # will likely still need to provide a value for each argument to get your code to run.
#
# # Task 1: Parsing Command-line Arguments. Uncomment this section when you have completed the my_parse_args function
# # and ensure the args that print match your expectations. This section also checks that the RE1 argument you provide
# # is a valid restriction enzyme name from the rest_dict
print('testing my_parse_args function')
args = my_parse_args()
print(f'the following arguments were received: {args}')
if args.RE1 not in rest_dict:
    print(f'no restriction enzyme named {args.RE1} found in rest_dict. Run rest_dict.keys() to see valid names')
print('\n')
#
#
# # Task2: Reading the Fasta File. Uncomment this section when you have completed the read_fasta function.
print('testing read_fasta function')
target_chromosome = 'NC_036780.1'
dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)
print(f'fasta file read. Length of {target_chromosome} is {len(dna_string)}')
print(f'length matches expected length of 45605991: {len(dna_string)==45605991}')
print('\n')
#
# # Task 3: Finding Motifs. Uncomment this section when you have completed the find_motifs function.
print(f'testing find_motifs function')
motif_locations = find_motifs('GGTTAAAGATCGGCGAGCCAATGGATCGACGATCA','GATC')
print(f'find_motifs returned {motif_locations}')
print(f'motif locations match expected locations of [7,23,30]: {motif_locations==[7,23,30]}')
print('\n')
#
# # Task 4: Implementing SingleRad. Uncomment this section when you have completed the run_single_rad function. Note that
# # this section hard-codes the RE1 argument for consistency, rather than using args.RE1 to get a value from the command
# # line
print('testing run_single_rad function')
srad_seq_sites = run_single_rad(dna=dna_string, re1='AanI')
print(f'located {len(srad_seq_sites)} singlerad sites for re1=AanI')
print(f'number of sites matches expected number (19280 sites): {len(srad_seq_sites) == 19280}')
print(f'first sequence site determined to be {srad_seq_sites[0]}')
print(f'first sequence site matches expected site (44450, 44650): {srad_seq_sites[0] == (44450, 44650)}')
print('\n')
#
#
# # Task 5: Implementing ddRad. Uncomment this section when you have completed the run_ddrad function. Note again that
# # the re1 and re2 arguments are hard-coded here for consistency. This section may  take longer to run than others
# # due to the complexity of ddrad, but should still run in well under a minute. If it doesn't, check your code
# # for infinite loops
print('testing run_ddrad function')
ddrad_seq_sites = run_ddrad(dna_string, re1='AanI', re2='MroI')
print(f'located {len(ddrad_seq_sites)} ddRad sites for re1=AanI, re2=MroI')
print(f'number of sites matches expected number (848 sites): {len(ddrad_seq_sites) == 848}')
print(f'first sequence site determined to be {ddrad_seq_sites[0]}')
print(f'first sequence site matches expected site (91143, 91243): {ddrad_seq_sites[0] == (91143, 91243)}')
print('\n')
#
# # Task 6: Finding Sites with Variation. Uncomment this section when you have completed the find_variable_sites function.
# # Again, this section will be slightly slower to run, but should still run in less than a minute
print('testing find_variable_sites function')
variable_srad_seq_sites = find_variable_sites(vcf_file_path=args.VCFFile, sequenced_sites=srad_seq_sites)
print(f'found {len(variable_srad_seq_sites)} singleRad sites with variation')
print(f'number of sites found matches expectation (2483 sites): {len(variable_srad_seq_sites)==2483}')
print(f'{len(100*variable_srad_seq_sites)/len(srad_seq_sites):.3f}% of sites sequenced with single rad had variation')
variable_ddrad_seq_sites = find_variable_sites(vcf_file_path=args.VCFFile, sequenced_sites=ddrad_seq_sites)
print(f'found {len(variable_ddrad_seq_sites)} ddRad sites with variation')
print(f'number of sites found matches expectation (48 sites): {len(variable_ddrad_seq_sites)==48}')
print(f'{len(100*variable_ddrad_seq_sites)/len(ddrad_seq_sites):.3f}% of sites sequenced with ddRad had variation')
#
#
# ### MAIN CODE ###
# # Once you think your code is complete, comment out the above testing section, uncomment the below section,
# # and run the script.
#
args = my_parse_args()
target_chromosome = 'NC_036780.1'
dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)
#
if args.Mode == 'SingleRad':
    print(f'running SingleRad for chromosome {target_chromosome} and restriction enzyme {args.RE1}')
    seq_sites = run_single_rad(dna_string, args.RE1)
elif args.Mode == 'ddRad':
    print(f'running ddRad for chromosome {target_chromosome} and restriction enzymes {args.RE1} and {args.RE2}')
    seq_sites = run_ddrad(dna_string, args.RE1, args.RE2)

input_len = len(dna_string)
n_sites = len(seq_sites)
sequencing_length = seq_sites[0][1] - seq_sites[0][0]
sequenced_fraction = sequencing_length*n_sites/input_len
print(f'{args.Mode} sequencing complete.')
print(f'\tlength of input sequence: {input_len}')
print(f'\tnumber of sequencing sites located: {n_sites}')
print(f'\tpercentage of nucleotides sequenced: {100*sequenced_fraction:.3f}%')

print('checking for variation within sequenced sites')
variable_sites = find_variable_sites(args.VCFFile, sequenced_sites=seq_sites)
n_var_sites = len(variable_sites)
var_fraction = n_var_sites / n_sites
print(f'\tnumber of sites with variation: {n_var_sites}')
print(f'\tfraction of sites with variation: {var_fraction}')

print('\nAnalysis Complete')

print('Complete!')