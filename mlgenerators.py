#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#############################################
##   lankapalli@shh.mpg.de
##   Program runs RAxML and TREE-PUZZLE
##   programs and evalutes the output of
##	TREE-PUZZLE output
##
#############################################
import os,sys
from alignment import *
import multiprocessing as mp
from collections import OrderedDict



import datetime

import subprocess
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

##Estimate the p-value or Confidence set
def __pval__(value):
	return 1 if value<=0.05 else 0

def __sval__(value):
	return 1 if value=='-' else 0


def RAXML(cmd):
	g=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
	retVal=0
	for line in g.stdout.readlines():
		g.wait()
		print(line.decode('utf-8'))
	retVal=g.returncode
	print("exit code {0}".format(retVal))

## check if previous files exist 
def RAxMLCheck(cmd):
	print("Checking for previous RAxML run")
	cmd_name=cmd.split(' ')[(cmd.split(' ').index('-n'))+1]
	print(cmd_name)
	if os.path.exists('RAxML_info.'+str(cmd_name)):
		print("Warning : {0} file already exists and will be renamed as {1}".format(cmd_name,'Old_'+cmd_name))
	return cmd_name

##Process a RAxML run

def raxml_run(backbone_cmds):
	cmd=backbone_cmds
	print(cmd)
	cmd_outname=RAxMLCheck(cmd)
	starttime=datetime.datetime.now()
	print("Starting backbone tree reconstruction at {0}".format(starttime.strftime('%H:%M:%S')))
	RA=mp.Process(target=RAXML,args=(cmd,))
	RA.start()
	RA.join()
	RA.terminate()
	endtime=datetime.datetime.now()
	process_elapsed=elapsed_time(endtime-starttime)
	print("Backbone tree reconstruction sucessfully completed at {0}".format(endtime.strftime('%H:%M:%S')))
	print("Time for {0} analysis took {1}".format(raxml_run.__name__,process_elapsed.timestring()))



def TREEPUZZLE(cmd,cmdname,q):
	try:
		g=subprocess.run(cmd,stderr=subprocess.PIPE,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		q.put((cmdname,g))
	except:
		raise("Something is wrong !!, Please check the ")

def create_table(g,numT,string='Tree   log L   difference    S.E.      p-1sKH     p-SH       c-ELW      2sKH'):
	''' Identifies TREE-PUZZLE table output for clock-like tree models and non clock-like tree models '''
	try:
		assert(g.count(string)==1)
		h=g[(g.index(string)+2):(g.index(string)+2+numT)]
		k=[list(filter(lambda x:x!='<----',filter(None,h[n].strip().split(' ')))) for n,i in enumerate(h)]
		return(k)
	except Exception:
		print("{0} was found {1} times".format(string,g.count(string)))

def puzzletable(k):
	'''For each TREE-PUZZLE table, the function returns highest Log-likelihood tree model , the Log-likelihood difference between average Log-likelihood of non-rejected models and average Log-likelihood of rejected models, ELW statistics for the all tree models,sorted order of tree models based on Log-likelihood scores,sorted Log-likelihood scores'''
	Logl,Diff,SE,p1KH,SH,ELW=([] for i in range(6))
	for i in k:
		Logl+=[('Tree_'+str(int(i[0])),float(i[1]))]
		Diff+=[('Tree_'+str(int(i[0])),float(i[2]))]
		p1KH+=[('Tree_'+str(int(i[0])),(__pval__(float(list(filter(None,i[4].split(' ')))[0])),__sval__((list(filter(None,i[5].split(' ')))[0]))))]
		SH+=[('Tree_'+str(int(i[0])),(__pval__(float(list(filter(None,i[6].split(' ')))[0])),__sval__((list(filter(None,i[7].split(' ')))[0]))))]
		ELW+=[('Tree_'+str(int(i[0])),(__pval__(float(list(filter(None,i[8].split(' ')))[0])),__sval__((list(filter(None,i[9].split(' ')))[0]))))]
		if i[3]=='best':
			BestTree=[('Best','Tree_'+str(int(i[0])))]
			SE+=[(int(i[0]),0.0)]
		else:
			SE+=[(int(i[0]),float(i[3]))]
	return(Logl,Diff,SE,BestTree+p1KH,BestTree+SH,BestTree+ELW)

def Logvaluesbased(Logl):
	Logl=OrderedDict(Logl)
	di=lambda x,y : x-y
	SS=sorted(Logl,key=lambda x: Logl[x],reverse=True)
	SA=list(map(lambda x: Logl[x],SS))
	bestllh_diff=[di(sorted(list(set(SA)),reverse=True)[0],sorted(list(set(SA)),reverse=True)[1]) if len(sorted(list(set(SA)),reverse=True))>1 else 0][0]
	return(bestllh_diff,SS,Logl)
## Wrapper function to read puzzle output files and return results from puzzletable function
def puzzle_table(g,numTree):
	@readfiles
	def puzzle_read(g):
		return g
	return puzzletable(create_table(puzzle_read(g),numTree))


@readfiles
def raxml_parameters(g):
	puzzledetails={i.split(':')[0]: i.split(':')[1] for i in list(filter(None,g[g.index('Model Information:')-2:g.index('Model Information:')+17]))}
	return puzzledetails


# Main process to run tree puzzle setup
def running_process(align,start,end,prefix,treefilename,puzzleparametersfilename):
	nametag=prefix+'_'+str(start)+'_'+str(end)
	outfile=open(nametag+'.phy','w')
	phylipalign=phylip_format(slidinator(align,start,end))
	outfile.write(phylipalign)
	outfile.close()
	return (nametag,'/usr/local/bin/puzzle '+str(nametag)+'.phy '+str(treefilename)+' -param='+str(puzzleparametersfilename)+' -wsl -wsr -wsitefreqs -prefix='+str(nametag)+'.out')




def puzzle_run(align,loc,prefix,treefilename,puzzleparametersfilename):
	puzzle_commands={}
	for start,end in loc:
		command_name,command=running_process(align,start,end,prefix,treefilename,puzzleparametersfilename)
		puzzle_commands[command_name]=command
	return puzzle_commands



def TREEPUZZLE1(commd):
	cmd,cmdname,q=commd
	print(cmd)
	print("Process id {0}".format(os.getpid()))
	g=subprocess.run(cmd.split(' '),stderr=subprocess.PIPE,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	q.put([cmdname,g.returncode])




