# coding=utf-8

import os
import re
import sys
import time

####### FUNCTION DEFINITIONS #######

def open_file(x):           # opens a file given as an absolute or relative path
    a = open(x).read().splitlines()
    return a

def names(o):       # returns a list containing sequences' names
    s = open(o).read()
    names = re.findall(">.+\n",s)
    return names

def sep_seq(o):     # returns a list containing 6 grouped peptid sequences by contig (its translation for each frame)
    s = open(o).read()
    seqs = re.split(">.+\n",s)
    seqs.pop(0)
    sep = [seqs[x:x+6] for x in range(0,len(seqs),6)]
    return sep

def split_orfs(l):  # splits ORFs for each frame translation for each contig
    split = [[[] for x in range(6)] for x in range(len(l))]
    for i in l:
        ind = l.index(i)
        indf = int((float(ind)/float(len(l)))*100)        
        if ind%30000 == 0:
            print("\t\t\t\t"+time.ctime())
        if ind%3000 == 0:
            print "Splitting ORFs: ", indf, "%"
        for s in i:
            split[l.index(i)][i.index(s)] = s.split("*")
    return split

def longest_orf(l):     # finds the longest ORF per frame
    longest = []
    for h in l:
        for i in h:
            if i:
                longest.append(max(i,key = len))
            else:
                longest.append("ERROR")
        ind = l.index(h)
        indf = int((float(ind)/float(len(l)))*100)
        if ind%3000 == 0:
            print "Assigning longest ORF per frame: ", indf, "%"
        if ind%30000 == 0:
            print("\t\t\t\t"+time.ctime())
    return longest

def join(n,s):      # reunites each frame's name with its longest ORF
    tot = []
    for i in n:
        ind = n.index(i)
        indf = int((float(ind)/float(len(n)))*100)
        tot.append(n[ind]+s[ind])
        if ind%5000 == 0:
            print "Rejoining names and sequences: ", indf, "%"
        if ind%50000 == 0:
            print("\t\t\t\t"+time.ctime())
    return tot

def annotation(ss,nn):   # annotates the dataset ss with reference nn
    for i in ss:
        for j in nn:
            if i[1] in j[0]:
                i.append(j[1])        
        ind = ss.index(i)
        indf = int((float(ind)/float(len(ss)))*100)
        if ind%5000 == 0:
            print "Performing annotation:", indf,"%"
        if ind%50000 == 0:
            print("\t\t\t\t"+time.ctime())
    return ss


######## STATEMENTS ########

fasta = sys.argv[1]
ref_blat = sys.argv[2]
blat = sys.argv[3]
transeq = sys.argv[4]
merge = sys.argv[5]

print("STARTING ANNOTATION PIPELINE FOR %s" % fasta)

pep = fasta[:-6]+'_six_frames.pep'
query_blat = pep[:-4]+'_longest_clean.pep'
result_blat = query_blat[:-4]+'_vs_prot.psl'
def_file = fasta[:-6]+'_def'
sorted_file = fasta[:-6]+'_clean_sorted'

os.system(transeq+' '+fasta+' '+pep+' -frame=6')					######## the translation
print("\n\n\t\tTranslated data file: %s\n\n" % pep)

seq_names = names(pep)
sep = sep_seq(pep)
split = split_orfs(sep)
print "\n\n\t\tORF splitting done\n\n"

for i in split:
    if i == [[] for x in range(6)]:
        ll = [sep[split.index(i)],sep[split.index(i)+1]]
        split2 = orf.split(ll)
        split[split.index(i)] = split2[0]

longest = longest_orf(split)											############ longest ORF selection
print "\n\n\t\tLongest ORF selection done\n\n"

seq_clean = [re.sub('\n','',i) for i in longest]
tot = join(seq_names, seq_clean)
with open(query_blat,'a') as file:
    for i in tot:
        print >> file, i

print("\n\n\t\tLongest ORF per frame file: %s\n\n" % query_blat)

print "\n\n\t\tStarting BLAT\n\n"
os.system('%s -oneOff=1 -noHead -prot %s %s %s' % (blat, ref_blat, query_blat, result_blat))		######### blat
print("\n\n\t\tBLAT done: result %s\n\n" % result_blat)

os.system('php %s %s > %s' % (merge, result_blat, result_blat+'_merge'))            ######## blat correction

n1 = names(ref_blat)
n = [re.sub('\n','',i) for i in n1]
n2 = [re.sub('>','',i) for i in n]
nn = [re.split(' ',i,maxsplit = 1) for i in n2]         #### prepares the reference for annotation

q = [re.split('\t',i) for i in open_file(result_blat+'_merge')]

ssa = annotation(q,nn)                       ##### annotation
print "\n\n\t\tANNOTATION DONE\n\n"

### arranging and sorting of annotation

ss2 = [i[:7] if len(i)>7 else i for i in ssa]
sep = [i for i in ss2 if len(i) == 7]

query_id = [i[0] for i in sep]
target_id = [i[1] for i in sep]
query_len = [i[2] for i in sep]
target_len = [i[3] for i in sep]
align_len = [i[4] for i in sep]
score = [i[5] for i in sep]
target_anot = [i[6] for i in sep]

id_2 = [i[3:] for i in target_id]
target_id = [re.sub('\|.+$','',i) for i in id_2]
target_os = [re.split('=',i[-1])[1][:-3] for i in ssa]
target_short = [re.split('=',i[-1])[0][:-3] for i in ssa]

frames = [re.findall("[\d]$",i) for i in query_id]
frames_def = [''.join(i) for i in frames]

new = [query_id[i][:-2]+'\t'+frames_def[i]+'\t'+score[i]+'\t'+target_id[i]+'\t'+target_len[i]+'\t'+target_short[i]+'\t'+target_os[i] for i in range(len(score))]

with open(fasta[:-6]+'_clean_presorted','a') as file:
    for i in new:
        print >> file, i

os.system('sort -k1,1 -k2,2 -k3,3 %s > %s' % (fasta[:-6]+'_clean_presorted',sorted_file))
os.system('rm %s' % fasta[:-6]+'_clean_presorted')
print("\n\n\t\tAnnotated and sorted file: %s\n\n" % sorted_file)

###

### best hit selection

neww = [re.sub('\t','',i,count = 1) for i in open_file(sorted_file)]
sepp = [re.split('\t',i) for i in neww]
for i in sepp:
    i[1] = float(i[1])
    i[3] = int(i[3])

lil = [[sepp[0]]]
for li in sepp[1:]:
    if li[0] == lil[-1][-1][0]:
            lil[-1].append(li)
    else:
            lil.append([li])

single = [i[-1] for i in lil]

single_sep = [[a[0][:-1],a[0][-1],a[1],a[2],a[3],a[4],a[5]] for a in single]
single_sep_unique = [[single_sep[0]]]
for li in single_sep[1:]:
    if li[0] == single_sep_unique[-1][-1][0]:
            single_sep_unique[-1].append(li)
    else:
            single_sep_unique.append([li])

single_def = []
for i in single_sep_unique:
    if len(i) == 1:
        single_def.append(i)
    elif len(i) == 2:
        if i[0][2]>i[1][2]:
            single_def.append(i[0])
        else:
            single_def.append(i[1])            
    elif len(i) == 3:
        if i[0][2]>i[1][2] and i[0][2]>i[2][2]:
            single_def.append(i[0])
        elif i[1][2]>i[0][2] and i[1][2]>i[2][2]:
            single_def.append(i[1])
        elif i[2][2]>i[0][2] and i[2][2]>i[1][2]:
            single_def.append(i[2])
        else:
            print "\nTied contig: ", i, "\n"
    else:
        print "\nContig with more than 3 frames in blat: ", i, "\n"

single_deff = [i[0] if type(i[0]) == list else i for i in single_def]
for i in single_deff:
    i[2]=str(i[2])
    i[4]=str(i[4])

def_tab = ['\t'.join(i) for i in single_deff]

with open(def_file,'a') as file:
	for i in def_tab:
		print >> file, i

print("\n\n\t\tFINAL ANNOTATED FILE: %s\n\n" % def_file)
