#usage script.py listfile pathtodirectory lengththresh

#listfile is a file with list of specimens. The script will find file for each specimen+".fa" and perform all merge. Final output is in the listfile+"_summarystats" file.
from collections import Counter
import sys,os

def allmerge(infilename,thresh):
	print infilename
	infile=open(infilename)
	def reverse_comp(inseq):
		comp_seq=''
		flag=True
		ambiguity_codes=['N','K','M','R','Y','S','W','B','V','H','D','N','X','n','k','m','r','y','s','w','b','v','h','d','x']
		for base in ambiguity_codes:
			if base in inseq:
				flag==False
		if flag==True:	
			for nucl in inseq:
				if nucl=='A':
					comp_seq=comp_seq+'T'
				elif nucl=='T':
					comp_seq=comp_seq+'A'
				elif nucl=='G':
					comp_seq=comp_seq+'C'
				elif nucl=='C':
					comp_seq=comp_seq+'G'
		elif flag==False:
			comp_seq=0
		return comp_seq[::-1]
		
	def findsubsequence(headseq,seqdict,dmplist,n):
		dmplist.append(headseq)
		subset=[headseq]
		totalcounts=seqdict[headseq]
		maxcounts=seqdict[headseq]
		dominantseq=headseq
		seqdict[headseq]=0
		for j in seqdict.keys():
			for ind in subset:
				if j!=ind:
					if len(j)>len(ind):	
						if ind in j:
							n=n+1
							totalcounts=totalcounts+seqdict[j]
							subset.append(j)
							if seqdict[j]>maxcounts:
								maxcounts=seqdict[j]
								dominantseq=j
							dmplist.append(j)
							seqdict[j]=0
							break
						if ind in reverse_comp(j):
							n=n+1
							totalcounts=totalcounts+seqdict[j]
							subset.append(reverse_comp(j))
							if seqdict[j]>maxcounts:
								maxcounts=seqdict[j]
								dominantseq=j
							dmplist.append(j)
							seqdict[j]=0
							break
					elif len(j)<len(ind):
						if j in ind:
							n=n+1
							totalcounts=totalcounts+seqdict[j]
							subset.append(j)
							if seqdict[j]>maxcounts:
								maxcounts=seqdict[j]
								dominantseq=j
							dmplist.append(j)
							seqdict[j]=0
							break
						elif j in reverse_comp(ind):
							n=n+1
							totalcounts=totalcounts+seqdict[j]
							subset.append(reverse_comp(j))
							if seqdict[j]>maxcounts:
								maxcounts=seqdict[j]
								dominantseq=j
							dmplist.append(j)
							seqdict[j]=0
							break
		return totalcounts,subset,maxcounts,dominantseq,n

	seqlist=[]
	for line in infile.readlines():
		if ">" not in line:
			if len(line.strip())>thresh:
				seqlist.append(line.strip())

	seqlistuniques={}
	for seq in seqlist:
		if seq not in seqlistuniques.keys():
			if reverse_comp(seq) not in seqlistuniques.keys():
				seqlistuniques[seq]=0

	for each in seqlist:
		try:
			seqlistuniques[each]=seqlistuniques[each]+1
		except KeyError:
			seqlistuniques[reverse_comp(each)]=seqlistuniques[reverse_comp(each)]+1

	o1=open(infilename+".uniq",'w')
	for i,each in enumerate(seqlistuniques.keys()):
		o1.write(">Unique_"+str(i)+"; "+str(seqlistuniques[each])+'\n'+each+'\n')
	o1.close()
	o2=open(infilename+".uniq.5",'w')
	for i,each in enumerate(seqlistuniques.keys()):
		if seqlistuniques[each]>=5:
			o2.write(">Unique_"+str(i)+"; "+str(seqlistuniques[each])+'\n'+each+'\n')
	o2.close()

	seqlistunique_sorted= sorted(seqlistuniques.items(), key=lambda x: x[1])
	countsvec_unique5 = []
	for i in [x[1] for x in seqlistunique_sorted]:
		if i>=5:
			countsvec_unique5.append(i)
	if len(countsvec_unique5)>0:
		countsvec_unique5.sort()
		countsvec_unique5=countsvec_unique5[::-1]

	params={}
	dmp=[]
	dmplist=[]
	n=0
	for each in seqlistuniques.keys():
		if seqlistuniques[each]>=5:
			if each not in dmplist:
				counts,seqset,maxcounts,dominantseq,n=findsubsequence(each,seqlistuniques,dmplist,n)
				params[each]=[counts,seqset,maxcounts,dominantseq,n]
				n=params[each][4]
				dmp.append([params[each][3],params[each][0]])
	print str(len(dmp)) + " sequences left"
	if len(dmp)>0:
		v=[]
		for each in dmp:
			v.append(each[1])

		print max(v)
		o=open(infilename+".merged",'w')
		for i,each in enumerate(dmp):
			o.write(">Merged_"+str(i)+"; "+str(each[1])+'\n'+each[0]+'\n')
		o.close()

		dmp= sorted(dmp, key=lambda x: x[1])

		if len(dmp)==1:
			seqoutset_fordisplay=[x[0] for x in dmp]
		if len(dmp)==2:
			seqoutset_fordisplay=[x[0] for x in dmp[-2:]]
		if len(dmp)>=3:
			seqoutset_fordisplay=[x[0] for x in dmp[-3:]]

		countsvec_merge5 = [x[1] for x in dmp]
		print countsvec_merge5
		if len(countsvec_merge5)>0:
			countsvec_merge5.sort()
			countsvec_merge5=countsvec_merge5[::-1]
		seqoutset_fordisplay=seqoutset_fordisplay[::-1]
	else:
		countsvec_merge5=[]
		seqoutset_fordisplay=[]
		
	return countsvec_unique5,countsvec_merge5,seqoutset_fordisplay
	
listfile=open(sys.argv[1])
l=listfile.readlines()
l2=[]
for x in l:
	x=x.strip()
	if len(x)>0:
		l2.append(x)
print l2
path=sys.argv[2]
thresh=int(sys.argv[3])
dirlist=os.listdir(path)
print thresh
outfile=open(sys.argv[1]+"_summarystats",'w')
outfile.write("Specimen ID\tCounts of unique sequences(>5)\tCounts of merged sequences (>5)\tFirst dominant sequence\t Second dominant sequence\t Third Dominant sequence\n")
for name in l2:
	countsunique,countsmerge,top3= allmerge(path+name+".fa",thresh)
	countsunique=[str(x) for x in countsunique]
	countsmerge=[str(x) for x in countsmerge]
	outfile.write(name+'\t'+','.join(countsunique)+'\t'+','.join(countsmerge))
	for each in top3:
		outfile.write('\t'+each)
	outfile.write('\n')
outfile.close()