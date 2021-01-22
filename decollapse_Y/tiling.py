#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

import pysam
from sys import argv
import csv

#bam
ref = argv[1]
#ref
bam = argv[2]
#bed
bed = argv[3]
#segdups
seg = argv[4]

def check_seq(seq, coordinates):
	return_seq = True
	for qstart, qend in coordinates.items():
		if seq.reference_start>qend or seq.reference_end<qstart:
			return_seq = True
		else:
			return_seq = False
			break
	if return_seq is True:
		coordinates.update({seq.reference_start:seq.reference_end})
	return return_seq
	
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

samfile = pysam.AlignmentFile(bam, "rb")

seqid = samfile.get_reference_name(0)

name_suffix = seqid+"_decollapsed"

print("Reference: "+samfile.get_reference_name(0))
print ("N of decollapsed sequences: "+str(samfile.count()))

seq = samfile.fetch()
seq_sorted_by_query_alignment_length = sorted(seq, key=lambda x: x.query_alignment_length, reverse=True)

samfile.close()

print("\nSummary:")
for seq in seq_sorted_by_query_alignment_length:
    print (str(seq.query_name)+"\t"+str(seq.query_alignment_length)+"\t"+str(seq.reference_start)+"\t"+str(seq.reference_end))

first_block_start = seq_sorted_by_query_alignment_length[0].reference_start
first_block_end = seq_sorted_by_query_alignment_length[0].reference_end

print("\nSeed:")
print(seq_sorted_by_query_alignment_length[0].query_name+"\t"+str(seq_sorted_by_query_alignment_length[0].query_alignment_length)+"\t"+str(first_block_start)+"\t"+str(first_block_end))

coordinates = {}
tiles = []

f0 = open(seqid+"_unplaced.fasta", "w")

for seq in seq_sorted_by_query_alignment_length:  
	new_seq = check_seq(seq, coordinates)
	if new_seq==True:
		tiles.append(seq)
	else:
		f0.write(">"+seqid+"_"+seq.query_name+" "+str(seq.reference_start)+"-"+str(seq.reference_end)+"\n")
		f0.write(seq.query_sequence+"\n")
		
f0.close()

sorted_tiles = sorted(tiles, key=lambda x: x.reference_start, reverse=False)

seq_bed = []
ref_bed = []

print("\nSegdup tiles (mapped bases):")		
for x in range(len(sorted_tiles)): 
    print(str(sorted_tiles[x].query_name)+"\t"+str(sorted_tiles[x].reference_start)+"\t"+str(sorted_tiles[x].reference_end))
    seq_bed.append([sorted_tiles[x].query_name,sorted_tiles[x].reference_start,sorted_tiles[x].reference_end]) 

print("\nReference tiles (1-based):")
with open(bed, newline = '') as bedfile: 
	bedfile = csv.reader(bedfile, delimiter='\t')
	bedfile = ([column[0], int(column[1]), int(column[2])] for column in bedfile)
	for line in bedfile:
		line[1]+=1
		line[2]+=1
		print(line[0]+"\t"+str(line[1])+"\t"+str(line[2]))
		ref_bed.append(line)

c1=1
c2=0
n=1
gap=1000

segdups = pysam.FastaFile(seg)

f1 = open(name_suffix+".agp", "w")
f2 = open(name_suffix+".fasta", "w")

f2.write(">"+name_suffix+"\n")

for x in range(len(seq_bed)-1):
	for y in range(len(ref_bed)):
		if ref_bed[y][2]<=seq_bed[0][1] and x==0:
			c2=ref_bed[y][2]-1
			length=c2-ref_bed[y][1]+1
			f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tW\t"+str(ref_bed[y][0])+"\t"+str(ref_bed[y][1])+"\t"+str(length)+"\t+\n")
			f2.write(segdups.fetch(str(ref_bed[y][0]), start=ref_bed[y][1]-1, end=ref_bed[y][2]-1))
			n=n+1
			c1=c2+1
			c2=c2+gap
			f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tN\t"+str(gap)+"\t"+"scaffold\tyes\tna\n")
			f2.write("N"*1000)
			n=n+1
			c1=c1+gap
	length=sorted_tiles[x].query_length
	c2=c2+length
	if sorted_tiles[x].is_reverse:
		sign="-"
		dna = revcomp(segdups.fetch(sorted_tiles[x].query_name))
	else:
		sign="+"
		dna = segdups.fetch(sorted_tiles[x].query_name)
	f1.write(str(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tW\t"+sorted_tiles[x].query_name)+"\t1\t"+str(length)+"\t"+sign+"\n")
	f2.write(dna)
	n=n+1
	c1=c2+1
	c2=c2+gap
	f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tN\t"+str(gap)+"\t"+"scaffold\tyes\tna\n")
	f2.write("N"*1000)
	n=n+1
	c1=c1+gap
	for y in range(len(ref_bed)):
		if ref_bed[y][1]>=seq_bed[x][2] and ref_bed[y][2]-1<=seq_bed[x+1][1]:
			length=ref_bed[y][2]-ref_bed[y][1]
			c2=c2+length
			f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tW\t"+str(ref_bed[y][0])+"\t"+str(ref_bed[y][1])+"\t"+str(ref_bed[y][2]-1)+"\t+\n")
			f2.write(segdups.fetch(str(ref_bed[y][0]), start=ref_bed[y][1]-1, end=ref_bed[y][2]-1))
			n=n+1
			c1=c2+1
			c2=c2+gap
			f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tN\t"+str(gap)+"\t"+"scaffold\tyes\tna\n")
			f2.write("N"*1000)
			n=n+1
			c1=c1+gap
length=sorted_tiles[len(seq_bed)-1].query_length
c2=c2+length
if sorted_tiles[len(seq_bed)-1].is_reverse:
	sign="-"
	dna = revcomp(segdups.fetch(sorted_tiles[len(seq_bed)-1].query_name))
else:
	sign="+"
	dna = segdups.fetch(sorted_tiles[len(seq_bed)-1].query_name)
f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tW\t"+str(sorted_tiles[len(seq_bed)-1].query_name)+"\t1\t"+str(length)+"\t"+sign+"\n")
f2.write(dna)

n=n+1
c1=c2+1

for y in range(len(ref_bed)):
	if ref_bed[y][1]>=seq_bed[len(seq_bed)-1][2]:

		c2=c2+gap		
		f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tN\t"+str(gap)+"\t"+"scaffold\tyes\tna\n")
		f2.write("N"*1000)
		n=n+1
		c1=c1+gap

		length=ref_bed[y][2]-ref_bed[y][1]
		c2=c2+length
		f1.write(name_suffix+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(n)+"\tW\t"+str(ref_bed[y][0])+"\t"+str(ref_bed[y][1])+"\t"+str(ref_bed[y][2]-1)+"\t+\n")
		f2.write(segdups.fetch(str(ref_bed[y][0]), start=ref_bed[y][1]-1, end=ref_bed[y][2]-1))
		c1=c2+1
		n=n+1

		
f1.close()
f2.close()