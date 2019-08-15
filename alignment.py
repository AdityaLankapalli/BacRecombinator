#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import random

### read the alignment file and genome header list for analysis
### read the contents of alignment file into alignment dictionary with 'header names as key' and 'sequence as value'



def readfiles(func):
	def read_file(*args,**kwargs):
		with open(*args,**kwargs) as f:
			g=f.read().strip().split('\n')
		result=func(g)
		return result
	return read_file



def readfasta(func):
	def read_file(*args,**kwargs):
		with open(*args,**kwargs) as f:
			g=f.read().strip().split('>')
		result=func(g)
		return result
	return read_file



@readfasta
def fastafile(g):
	fasta_alignment={}
	for fasta in list(filter(None,g)):
		fasta_alignment[fasta.split('\n')[0]]=''.join(fasta.split('\n')[1:])
	return fasta_alignment


keywords=['Recipients','Donars','Targets','Sources','RecombiningGenomes','BackboneGenomes']


def addlist(wo):
	l=[]
	for i in wo:
		l+=i
	return l

@readfiles
def genomedetails(genome_list):
	keywords=['Recipients','Donars','Targets','Sources','RecombiningGenomes','BackboneGenomes']
	Recipients=[]
	Donars=[]
	for x in list(filter(None,genome_list)):
		if x in keywords:
			A=[]
			if keywords.index(x)%2==0:
				Recipients.append(A)
			elif keywords.index(x)%2==1:
				Donars.append(A)
		else:
			A.append(x)
	return addlist(Recipients),Donars

		

#raise genomesError("something is wrong no keywords found !! Please check again")


def header_modifinator(alignment):
	random_chars=''.join(random.sample('ABCDEFGHIJKLMNOPQRSTUVWXYZ',5))
	headeralias={}
	num=random.randint(0,100)
	for i in alignment:
		headeralias[i]=random_chars+str('%0'+str(5)+'d') % num
		num+=1
	return headeralias

def aliasfindinator(id,headeralias):
	if id in list(headeralias.keys()):
			return headeralias[str(id)]
	elif id in list(headeralias.values()):
		try:
			list(headeralias.values()).count(str(id))==1
			return list(headeralias.keys())[list(headeralias.values()).index(str(id))]
		except DuplicateError as err:
			raise ( "Duplicate error: {0} found multiple times".format(err))

'''#check error statement is it ###'''


## fun code ##
def headernames_modifier(i):
	'''Creates a temporary header name for given name based on the character size of the name. if less than
		10 a random number with leading zeros is added else the string is stripted to character size of 8
		and random number is added. Total character size of return header name is 10'''
	if len(i)<=10:
		return(str(i)+((str('%0'+str((10-len(i)))+'d') % (random.randint(0,10**(10-len(i)))))))
	else:
		return(i[:8]+str(random.randint(0,100)))
#####

def slidinator(alignment,start,end):
	'''Based on the start and end position, a block of alignment is extracted'''
	block={}
	for i in alignment:
		block[str(i)]=alignment[str(i)][int(start):int(end)]
	return(block)



def mismatchinator(alignment,w,target=None):
	'''Estimates the composition of bases per site and estimates the variants contributed by the specific strain'''
	mismatcher=[]
	for j in range(0,w):
		m=[]
		for i in alignment:
			if alignment[str(target)][j] in ['A','C','T','G'] and alignment[str(i)][j] in ['A','C','T','G']:
				if len(list(set([alignment[str(target)][j],alignment[str(i)][j]])))==2:
					m.append(1)
				else:
					m.append(0)
			else:
				m.append(0)
		mismatcher.append(sum(m)/(len(alignment)-1))
	return(mismatcher)


def make_files(align,ser,we,wr,threads):
	backbone_cmds={}
	Backbone_target_sequences={}
	for num,i in enumerate(wr):
		backbonename=str(num)+'_Backbone.fasta'
		backbonetreename=str(num)+'_Backbone.tree'
		backbone_cmds[backbonetreename]='/Users/lankapalli/Downloads/standard-RAxML-master/raxmlHPC -s '+str(backbonename)+' -n '+str(backbonetreename)+' -m GTRGAMMAI -f a -p 123456 -x 123456789 -# 100 -T '+str(threads-1)
		o=open(backbonename,'w')
		[o.write('>'+str(aliasfindinator(j,ser))+'\n'+align[str(j)]+'\n') for num1,j in enumerate(i)]
		o.close()
		Backbone_target=[]
		for n,k in enumerate(we):
			backbonewithgenome=str(num)+'_Backbone_'+str(k)+'.fasta'
			print(backbonewithgenome)
			Backbone_target_sequences[backbonewithgenome]=[k]+i
			t=open(backbonewithgenome,'w')
			[t.write('>'+str(aliasfindinator(m,ser))+'\n'+align[str(m)]+'\n') for m in Backbone_target_sequences[backbonewithgenome]]
			t.close()
	return Backbone_target_sequences,backbone_cmds


def random_St_En(alignment_length):
	fn=lambda x,y: (y,x) if x>y else (x,y)
	a,b=fn(random.randint(0,(alignment_length-1)),random.randint(0,(alignment_length-1)))
	if b-a>100:
		return a,b
	else:
		random_St_En(alignment_length)


align_length=lambda align:list(set([len(align[i]) for i in align]))[0] if len(set([len(align[i]) for i in align]))==1 else MainError("alignment is of unequal length")





def slide_loc(alignment_length,sliding_window,sliding_overlap):
	loc=[]
	for blocki in range(0,alignment_length-sliding_window,sliding_overlap):
		loc.append((blocki,blocki+sliding_window))
	loc.append((blocki+sliding_overlap,alignment_length))
	return loc

def phylip_format(align):
	seq2str=lambda align:[str(x)+'\t\t'+str(align[x])+'\n' for x in align]
	return str(len(align))+'\t'+str(align_length(align))+'\n\n'+''.join(seq2str(align))






