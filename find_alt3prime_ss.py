import sys, csv
import scipy.stats as sps

try:
	bed1 = open(sys.argv[1])
	bed2 = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	print('usage: script.py junctionswt junctionsmt outfilename type')  # all junctions (not just the ones unique to long reads)
	sys.exit(1)

try:
	sys.argv[4]
	fiveprimeon = True
except:
	fiveprimeon=False


alljuncs = {}


for line in bed1:
	line = line.rstrip().split('\t')
	if fiveprimeon:
		if line[5] == '+':
			line[1], line[2] = line[2], line[1] 			

	elif line[5] == '-':
		line[1], line[2] = line[2], line[1]  # reverse coordinates for junctions on - strand
	chrom, fiveprime, threeprime, name, count, strand = line
	chrom = strand+chrom
	if chrom not in alljuncs:  # chrom
		alljuncs[chrom] = {}
	if fiveprime not in alljuncs[chrom]:
		alljuncs[chrom][fiveprime] = {}  # 5' end anchor
	if threeprime not in alljuncs[chrom][fiveprime]:
		alljuncs[chrom][fiveprime][threeprime] = [0,0]
	alljuncs[chrom][fiveprime][threeprime][0] = int(count)

for line in bed2:
	line = line.rstrip().split('\t')
	if fiveprimeon:
		if line[5] == '+':
			line[1], line[2] = line[2], line[1] 			
	elif line[5] == '-':
		line[1], line[2] = line[2], line[1]  # reverse coordinates for junctions on - strand
	chrom, fiveprime, threeprime, name, count, strand = line
	chrom = strand+chrom
	if chrom not in alljuncs:  # chrom
		alljuncs[chrom] = {}
	if fiveprime not in alljuncs[chrom]:
		alljuncs[chrom][fiveprime] = {}  # 5' end anchor
	if threeprime not in alljuncs[chrom][fiveprime]:
		alljuncs[chrom][fiveprime][threeprime] = [0,0]
	alljuncs[chrom][fiveprime][threeprime][1] = int(count)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in alljuncs:
		for fiveprime in alljuncs[chrom]:
			if len(alljuncs[chrom][fiveprime]) == 1:
				continue
			for threeprime1 in alljuncs[chrom][fiveprime]:
				if sum(alljuncs[chrom][fiveprime][threeprime1]) == 1:
					continue
				allothercounts = [0,0]
				for threeprime2 in alljuncs[chrom][fiveprime]:
					if threeprime1 == threeprime2: #or sum(alljuncs[chrom][fiveprime][threeprime2]) == 0:
						continue
					allothercounts[0] += alljuncs[chrom][fiveprime][threeprime2][0]
					allothercounts[1] += alljuncs[chrom][fiveprime][threeprime2][1]
					# ctable = [alljuncs[chrom][fiveprime][threeprime1], alljuncs[chrom][fiveprime][threeprime2]]
				if sum(allothercounts) == 1:
					continue
				ctable = [alljuncs[chrom][fiveprime][threeprime1], allothercounts]
				writer.writerow([chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], chrom[0]] + ctable[0] + ctable[1])

