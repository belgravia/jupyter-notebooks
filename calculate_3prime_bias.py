import sys, csv

try:
	txt = open(sys.argv[1])
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	print('usage: script.py junctionfile gtf outfilename')
	sys.exit(1)

annotpos = {}
annotmin = {}
for line in gtf:   # i really should write a gtf parser..
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	# gene = gene.split(';')[2]
	# gene = gene[len('gene_name') + 3:-1]
	if strand == '+':
		if chrom not in annotpos:
			annotpos[chrom] = {}
		annotpos[chrom][end] = {}
	else:
		if chrom not in annotmin:
			annotmin[chrom] = {}
		annotmin[chrom][end] = {}
sys.stderr.write('finished reading gtf\n')

lost = 0
for line in txt:
	line = line.rstrip().split('\t')
	chrom, start, end, counts, strand = line[0], int(line[1]), int(line[2]), int(line[4]), line[5]
	if strand == '+':
		wiggle = 0
		while end + wiggle not in annotpos[chrom]:
			if wiggle == 100:
				lost += 1
				break
			if wiggle == 0:
				wiggle += 1
			elif wiggle >= 0:
				wiggle = wiggle * -1
			else:
				wiggle = (wiggle-1) * -1
		if wiggle == 100:
			continue
		if wiggle in annotpos[chrom][end+wiggle]:
			annotpos[chrom][end+wiggle][wiggle] += counts
		else:
			annotpos[chrom][end+wiggle][wiggle] = counts
	elif strand == '-':
		wiggle = 0
		while start - wiggle not in annotmin[chrom]:
			if wiggle == 100:
				break
			if wiggle == 0:
				wiggle += 1
			elif wiggle >= 0:
				wiggle = wiggle * -1
			else:
				wiggle = (wiggle-1) * -1
		if wiggle == 100:
			lost += 1
			continue
		if wiggle in annotmin[chrom][start-wiggle]:
			annotmin[chrom][start-wiggle][wiggle] += counts
		else:
			annotmin[chrom][start-wiggle][wiggle] = counts			

distribution = {}
for chrom in annotpos:
	endpos = sorted(annotpos[chrom].keys())
	for e in endpos:
		# find the mode of wiggles for that end
		wiggle = 0
		counts = 0
		for w in annotpos[chrom][e]:
			if annotpos[chrom][e][w] >= counts and w != 0 and abs(w) < abs(wiggle):
				wiggle = w
				counts = annotpos[chrom][e][w]
		if wiggle == 0:
			continue
		if wiggle in distribution:
			distribution[wiggle] += 1
		else:
			distribution[wiggle] = 1
	endmin = sorted(annotmin[chrom].keys())
	for e in endmin:
		wiggle = 0
		counts = 0
		for w in annotmin[chrom][e]:
			if annotmin[chrom][e][w] > counts and w != 0 and abs(w) < abs(wiggle):
				wiggle = w
				counts = annotmin[chrom][e][w]
		if wiggle == 0:
			continue
		if wiggle in distribution:
			distribution[wiggle] += 1
		else:
			distribution[wiggle] = 1

dkeys = sorted(distribution.keys())
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for d in dkeys:
		writer.writerow([d, distribution[d]])



sys.stderr.write('lost ' + str(lost) + ' ends\n')
