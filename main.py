#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mlgenerators import *
from alignment import *
from Tree import *
import argparse,sys,os,logging,time
import datetime
import subprocess
import multiprocessing as mp
from random import random,seed,randint
from collections import Counter,OrderedDict
from pandas import DataFrame, concat
import gc

Maincommands = argparse.ArgumentParser(prog='BACRECOMBINATOR',usage='./BACRECOMBINATOR.py [-h] [-i INPUT_FASTA] [-o OUTPUT_FASTA] [-l SEQ_LENGTH] [-g] \n python3 BACRECOMINATOR.py [-h] [-i INPUT_FASTA] [-o OUTPUT_FASTA] [-l SEQ_LENGTH] [-g]',description='''Our automated pipeline detects possible regions of recombination using maximum likelihood tree topology variation in bacterial species.\
	Generate 	''',
								 epilog="Please send your suggestions and comments to Aditya  < lankapalli@shh.mpg.de > ")
Maincommands.add_argument("-i", "--input",dest="input_fasta", help="Input multifasta")
Maincommands.add_argument("-t", "--trees",dest="backbone_tree", help="Backbone topology/topologies file (Not applicable if raxml/iqtree active)")
Maincommands.add_argument("-p", "--parameterfile",dest="params_file", help="Backbone topology/topologies file (Not applicable if raxml/iqtree active)")
Maincommands.add_argument("-pre", "--prefix",dest="prefix", help="prefix for output")
Maincommands.add_argument("-sw", "--windowlength", dest="sliding_window", help="Size of sliding window for each alginment block", type=int,default=1000)
Maincommands.add_argument("-ovl","--overlap", dest="sliding_overlap", help="Overlap size of sliding windows", type=int,default=100)
Maincommands.add_argument("-s","--seed", dest="randomseed", help="Seed value for RNG", type=int)
Maincommands.add_argument("-g", "--rep_genomes", dest="rep_genomes", help="List of represntative genomes for backbone topology and target genome")
Maincommands.add_argument("-b", "--bootstrap", dest="bootstrap", help="Bootstrap replicates for Backbone Topology", type=int,default=1000)
Maincommands.add_argument("-op", "--options", dest="options", help="Additional options from the user", type=int,default=1000)
Maincommands.add_argument("-pro", "--program", dest="options", help="Program for Backbone Topology and topology testing", type=int,default='RAxML,TREE-PUZZLE')
Maincommands.add_argument("-nt", "--threads", dest="threads", help="Number of threads", type=int,default=3)
Maincommands.add_argument("-r","--redo", dest="redo_run", help="Re run the analaysis", action='store_false', default=True)
Maincommands.add_argument("-incr","--rootinclude", dest="rootinc", help="include root tree", action='store_false', default=True)
args = Maincommands.parse_args()

#pool=mp.Pool(threads=3)

class MainError(Exception):
	def __init__(self,message):
		self.message=message

class AlignmentError(MainError):
	def __init__(self,message):
		self.message=message

class PuzzleError(Exception):
	def __init__(self,message):
		self.message=message

class RAxMLError(Exception):
	def __init__(self,message,code):
		self.message=message
		self.code=code
	def __str__(self):
		return {self.code:self.message}


class TreeError(Exception):
	def __init__(self,value):
		if value is None:
			value='\n Please check for following errors in newick tree.\n'
		Exception.__init__(self,value)



class elapsed_time(object):
	def __init__(self,elapsed):
		self.elapsed=elapsed.total_seconds()
	def __days__(self):
		return divmod(self.elapsed,86400)[0]
	def __hours__(self):
		return divmod(self.elapsed,3600)[0]
	def __min_secs__(self):
		return divmod(self.elapsed,60)[0],divmod(self.elapsed,60)[1]
	def timestring(self):
		if self.elapsed>86400:
			return("%d days, %d hours, %d minutes and %d seconds" % (self.__days__(), self.__hours__(), self.__min_secs__()))
		elif self.elapsed > 3600 and self.elapsed < 86400 :
			return("%d hours, %d minutes and %d seconds" % ( self.__hours__(), self.__min_secs__()[0],self.__min_secs__()[1]))
		elif self.elapsed > 60 and self.elapsed < 3600 :
			return("%d minutes and %d seconds" % (self.__min_secs__()[0],self.__min_secs__()[1]))
		else:
			return("%d seconds" % (self.__min_secs__()[1]))







###commands


if args.randomseed:
	print("Seed value for this run {0}".format(args.randomseed))
	seed(args.randomseed)
else:
	seednumber=randint(0,1234567890)
	print("Seed value for this run {0}".format(seednumber))
	seed(seednumber)


######## original
def main():
	manager=mp.Manager()
	tasksqueue=manager.Queue()
	results=[]
	results1={}
	results2=[]
	results3=[]
	P=mp.Pool(processes=threads)
	P.map(TREEPUZZLE1,[(puzzle_commands[i],i,tasksqueue) for i in puzzle_commands])
	nume=0
	while not tasksqueue.empty():
		try:
			cmdname,returncode=tasksqueue.get()
			print(cmdname)
			print(returncode)
			if returncode == 0:
				Logl,ELW=puzzleresults(cmdname,prefix)
				results2.append(DataFrame(Logl,index=[nume]))
				results3.append(DataFrame(ELW,index=['pval','sval']))
				cleanup(cmdname)
				nume=nume+1
			else:
				results1[cmdname]=None
		except Exception as e:
			raise(e)
	P.close()
	P.join()
	return results1,concat(results2),concat(results3)



def puzzleresults(cmdname,prefix):
	start,end=list(map(int,cmdname.split(prefix)[1].split('_')[1:]))
	Logl,Diff,SE,p1KH,SH,ELW=puzzle_table((cmdname)+'.out.puzzle',numTree)
	Loglres=OrderedDict([('Start',start),('End',end)]+Logl)
	ELWpval=OrderedDict([('Start',start),('End',end)]+ELW)
	return Loglres,ELWpval



## duplicate code present
def startend(cmdname,prefix):
	start,end=list(map(int,cmdname.split(prefix)[1].split('_')[1:]))
	return start,end

def cleanup(cmdname):
	files2remove=['.out.puzzle','.out.tree','.out.sitelh','.out.sitefreqs','.phy','.out.siterate','.out.dist']
	'''files2remove=['.phy','.out.siterate','.out.dist']'''
	for file in files2remove:
		os.remove(cmdname+file)

# raxml_path=<path to raxml>

alignmentfile=args.input_fasta
genome_list=args.rep_genomes
rootinc=args.rootinc
threads=args.threads
align=fastafile(alignmentfile)
ser=header_modifinator(align)
we,wr=genomedetails(genome_list)
Backbone_target_sequences,backbone_cmds=make_files(align,ser,we,wr,threads)




##### Parameter setup
for raxmlcmds in backbone_cmds:
	raxml_run(backbone_cmds[raxmlcmds])
	puzzledetails=raxml_parameters('RAxML_info.'+raxmlcmds)
	gtrmodelorder={'rate A <-> C': '1\n', 'rate A <-> G': '2\n', 'rate A <-> T': '3\n', 'rate C <-> G': '4\n', 'rate C <-> T': '5\n', 'rate G <-> T': '6\n'}
	writeparameters=open(raxmlcmds.replace('.tree','_parameters.txt'),'w')
	writeparameters.write('e\nx\nx\nm\nm\n')
	for i in gtrmodelorder:
		writeparameters.write(gtrmodelorder[i]+puzzledetails[i].strip()+'\n')
	writeparameters.write('w\nw\nw\nc\n8\ny\n')
	writeparameters.close()
	Alltrees=build_trees('RAxML_bestTree.'+raxmlcmds,Backbone_target_sequences,ser,rootinc)



## main
for jam in Backbone_target_sequences:
	treefilename=str(Backbone_target_sequences[jam][0])+'_'+'Treetopologynewick.nwk'
	puzzleparametersfilename=str(jam.replace(Backbone_target_sequences[jam][0]+'.fasta','parameters.txt'))
	numTree=len(Alltrees[jam])
	print(Alltrees[jam])
	print(numTree)
	fasta_align=fastafile(jam)
	print("{} Tree topology models found in {}".format(numTree,treefilename))
	location=slide_loc(align_length(fasta_align),args.sliding_window,args.sliding_overlap)
	if args.prefix:
		prefix=args.prefix
	else:
		prefix=str(Backbone_target_sequences[jam][0])+'_run'
	puzzle_commands=puzzle_run(fasta_align,location,prefix,treefilename,puzzleparametersfilename)
	resultsa,resultsb,resultsc=main()
	missing=open(str(Backbone_target_sequences[jam][0])+'_'+'missingregions','w')
	for i in resultsa:
		mstart,mend=startend(i,prefix)
		missing.write(str(mstart)+'\t'+str(mend)+'\n')
	missing.close()
	ELWpval=resultsc.sort_values('Start').loc['pval',]
	ELWsval=resultsc.sort_values('Start').loc['sval',]
	ELWpval.index=list(range(ELWpval.shape[0]))
	ELWsval.index=list(range(ELWsval.shape[0]))
	likelihoodvalues=resultsb.sort_values('Start')
	likelihoodvalues.to_csv(prefix+'.trees_loglh',sep='\t', encoding='utf-8', index=False)
	ELWpval.to_csv(prefix+'.trees_ELWpvalues',sep='\t', encoding='utf-8', index=False)
	ELWsval.to_csv(prefix+'.trees_ELWstats',sep='\t', encoding='utf-8', index=False)


'''
	## prune method
from Tree import *

ww,EE=build_topologies(args.backbone_tree)



align=fastafile(alignmentfile)
ser=header_modifinator(align)
we,wr=genomedetails(genome_list)
## modify
## align,ser,Backbone_target_sequences,backbone_cmds=make_files(alignmentfile,genome_list,threads,target_genome=None)


for node in EE:
	print(node.name)
	if len(node.descendents)==0:
		rf,rf_index=prune_node(node)
		rootinc=False
		EF=[]
		listgen(ww,EF)
		o=open(str(node.name)+'_tree.nwk','w')
		tree_topology_creator(ww,EF,node,o,rootinc)
		del EF
		rewire_nodes(rf,rf_index,node)
		del rf,rf_index
	print("finished the trees")




'''

