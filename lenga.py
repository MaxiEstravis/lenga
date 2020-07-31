###### MIT License ######

###### Copyright (c) 2018 Maximiliano Estravis ######


######## MAKE SURE THAT THE SCRIPTS/ DIRECTORY IS ########
########## IN THE SAME DIRECTORY AS THIS SCRIPT ##########

# coding=utf-8

import os
import re
import time
import sys

def error(s):
    if e != 0:
        print('\n\n\t\t'+s+'\n\n')
        quit()

exper = sys.argv[1]
assembly = sys.argv[2]
blat = sys.argv[3]
exonerate = sys.argv[4]

print("\n\n\t\tBEGINNING WORK %s\n\n" % exper)


############## CUT CONTIGS' NAMES TO ONE WORD ##############

def names(o):
    s = open(o).read()
    names = re.findall(">.+\n",s)
    return names

def sep_seq(o):
    s = open(o).read()
    seqs = re.split(">.+\n",s)
    seqs.pop(0)
    return seqs

def join(n,s):
    tot = []
    for i in n:
        ind = n.index(i)
        tot.append(n[ind]+'\n'+s[ind])
    return tot

n=names(assembly)
s=sep_seq(assembly)

nn=[re.split(' ',i)[0] for i in n]
ss=[i[:-1] for i in s]

j=join(nn,ss)

with open(assembly+'_short_names.fasta','a') as file:
    for i in j:
        print >> file, i

assembly+='_short_names.fasta'

os.system('mkdir %s' % exper)
os.chdir(exper)

#####################################################
################ ASSEMBLY CORRECTION ################
#####################################################

os.system('mkdir correction_interm_files')      # creates a directory to keep intermediate files produced during assembly correction
os.system('cp ../%s correction_interm_files/assembly.fasta' % assembly)
os.chdir('correction_interm_files')

###### contig merging
###### generates 5 intermediate files: 
# 1) the reverse complement of all contigs (.fasta.rev)
# 2) all contigs, in both forward and reverse orientations (.fasta.fwd-rev)
# 3) a blat of 2) against itself (.fasta.fwd-rev.self.psl)
# 4) chain info file, containing which contigs were merged to produce each chain (.fasta.chains_info)
# 5) chain fasta file (.fasta.chains)
e= os.system('php ../../scripts/chains.php assembly.fasta %s %s' % (exonerate, blat))
error("ERROR IN CONTIG MERGING")

print("\n\n\t\tCONTIG MERGING DONE\n\n")

###### filtering
###### removes contigs used for chain generation, and keeps only chains and contigs longer than X nt
###### generates two intermediate files:
# 6) length-filtered chains (.fasta.chains.filter)
# 7) length-filtered contigs not used in chains (.fasta.filter)
X = "200"
e = os.system('php ../../scripts/filter_chains.php %s assembly.fasta.chains_info assembly.fasta.chains assembly.fasta' % X)
error("ERROR IN CHAIN AND CONTIG FILTERING")

print("\n\n\t\tCHAIN AND CONTIG FILTERING DONE\n\n")

##### first redundancy filtering
##### generates two intermediate files:
# 8) a blat of 6) against itself (.fasta.chains.filter.self.psl)
# 9) redundancy-filtered chains (.fasta.chains.filter.red1)
e = os.system('php ../../scripts/eliminate_revcomp.php assembly.fasta.chains.filter %s' % blat)
error("ERROR IN FIRST REDUNDANCY FILTERING")

print("\n\n\t\tFIRST REDUNDANCY FILTERING DONE\n\n")

##### reunite chains and contigs
##### generates one intermediate file:
# 10) length-filtered contigs not used in chains (file 7) and length- and redundancy-filtered chains (file 9) (.step01.fasta)
os.system('cat assembly.fasta.chains.filter.red1 > transcripts.step01.fasta')
os.system('cat assembly.fasta.filter >> transcripts.step01.fasta')

##### self-blat to eliminate further redundancy
##### generates one intermediate file:
# 11) a blat of 10) against itself (.step01.fasta.self.psl)
os.system('%s -maxIntron=0 -noHead -maxGap=0 -minScore=61 transcripts.step01.fasta transcripts.step01.fasta transcripts.step01.fasta.self.psl' % blat)

#### second redundancy filtering
#### OVL defines the minimum nt length evaluated to filter (in our case, the chosen k-mer size for assembly)
#### PERC defines the minimum overlap to filter 
#### generates one intermediate file:
# 12) redundancy filtered chains and contigs (.step02.fasta)
OVL = "67"
PERC = "0.80"
os.system('php ../../scripts/eliminate_redundant_fasta.php %s %s transcripts.step01.fasta transcripts.step01.fasta.self.psl > transcripts.step02.fasta' % (OVL, PERC))
error("ERROR IN SECOND REDUNDANCY FILTERING")

print("\n\n\t\tSECOND REDUNDANCY FILTERING DONE\n\n")

os.system('cp transcripts.step02.fasta ../%s.corrected' % assembly)

print("\n\n\t\tASSEMBLY CORRECTION DONE\n\n")
