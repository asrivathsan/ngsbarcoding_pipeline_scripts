#usage is 
#python scriptname input_threshold_fastafile input_listfile withinformation_on_var outfilename
#listfile should be structured as a simple list of variables of interest eg (remove # is next lines): 
#Jan
#Feb
#Mar

import sys,os,itertools
fastafile=open(sys.argv[1])
inlines=fastafile.readlines()
parfile=open(sys.argv[2])
pars=[x.strip() for x in parfile.readlines()]
print pars
clustdict={}
for line in inlines:
	if ">" in line:
		try:
			clustdict[line.split(" ")[1]]=[[],{}]
		except KeyError:
			print line
for line in inlines:
	if ">" in line:
		clustdict[line.split(" ")[1]][0].append(" ".join(line.split(" ")[2:]).strip())
print clustdict
for k,v in clustdict.iteritems():
	for par in pars:
		parcount=0
		for seq in v[0]:
			if par in seq:
				parcount+=1
		clustdict[k][1][par]=parcount
		
outfile= open(sys.argv[3],'w')
outfile.write('\t')
for each in pars:
	outfile.write('\t'+each)
outfile.write('\n')
for k,v in clustdict.iteritems():
	outfile.write(k.replace(";",""))
	for par in pars:
		outfile.write('\t'+str(clustdict[k][1][par]))
	outfile.write('\n')	
