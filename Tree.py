#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from alignment import *
import multiprocessing as mp
from numpy import random
from xkcdrgb import xkcd_rgb

class Node():
	def __init__(self,name=None,root=False):
		self.name=name
		self.ancestor='ROOT'
		self.descendents=[]
		self.branch_length=None
		self.root=root


class TreeError(Exception):
	def __init__(self,value):
		if value is None:
			value='\n Please check for following errors in newick tree.\n'
		Exception.__init__(self,value)


def newick_check(ns,index=0):
	try:
		ns[index]=='(' and ns[len(ns)-1]==';'
		ns.count('(')==ns.count(')')
	except TreeError:
		raise TreeError('Improper Newick format or newick tree doesnot have trailing ;')
	except TreeError:
		raise TreeError('Opening and closing parentesis do not match')
	finally:
		print('Newick tree looks good.')


#########################################



def Nodelabel(string,index):
	label=''
	branchlen=None
	while string[index]!=';':
		if string[index]==':':
			index+=1
			branchlen,index=branch(string,index)
			break
		label+=string[index]
		index+=1
	return label,branchlen,index

def branch(string,index):
	branchlen=''
	a=index
	while string[index]!=';':
		if string[index]==',' or string[index]==')':
			break
		branchlen+=string[index]
		index+=1
	if branchlen is not '':
		branchlen=float(branchlen)
	else:
		branchlen=None
	return branchlen,index


## adapted from pyevolve
def read_tree(tree):
	current_node=None
	index=0
	count=(tree.count(':')-tree.count(')'))+1
	while tree[index]!=';':
		if tree[index]=='(' and index==0:
			count=count
			root_node=Node(index,root=True)
			current_node=root_node
			index+=1
		elif tree[index]=='(':
			index+=1
			Internal_node=Node(index)
			Internal_node.ancestor=current_node
			current_node.descendents.append(Internal_node)
			current_node=Internal_node
		elif tree[index]==',':
			index+=1
		elif tree[index]==')':
			index+=1
			name,branchlength,index=Nodelabel(tree,index)
			if name=='' or name is None:
				current_node.name="Internalnode_"+str(count)
				current_node.branch_length=branchlength
			else:
				current_node.name=name
				current_node.branch_length=branchlength
			current_node=current_node.ancestor
			count+=1
		else:
			name,branchlength,index=Nodelabel(tree,index)
			leaf_node=Node()
			leaf_node.name=name
			leaf_node.branch_length=branchlength
			leaf_node.ancestor=current_node
			current_node.descendents.append(leaf_node)
	return root_node

######

def write_topology(root_node,n=0,s=[]):
	if n>0:
		s.append(",")
	if len(root_node.descendents) is not 0:
		s.append("(")
		[write_topology(i,n,s) for n,i in enumerate(root_node.descendents)]
		s.append(")")
	else:
		s.append(root_node.name)
	return ''.join(s)+';'










## Lambda function to return branch_length for a given node
ert=lambda x,y:x.branch_length if x.name==y else None

## returns a list of node and node names
def write_description(root_node,node,s=[]):
	if len(root_node.descendents) is not 0:
		[write_description(i,s) for i in root_node.descendents]
		s.append(ert(root_node,node.name))
	else:
		s.append(ert(root_node,node.name))
	return list(filter(None,s))


def colourup(numTrees):
	from seaborn import xckd_rgb
	colormap={}
	return 



def add_sister(S,node):
	EDS=Node()
	EDS.ancestor=S.ancestor
	EDS.descendents.append(node)
	node.ancestor=EDS
	EDS.descendents.append(S)
	S.ancestor=EDS
	for n,i in enumerate(EDS.ancestor.descendents):
		if i==S:
			EDS.ancestor.descendents[n]=EDS
		else:
			pass
	return EDS


#EDS.ancestor.descendents.append(S)
#S.ancestor.descendents.remove(EDS)

def remove_sister(EDS,node):
	EDS.descendents.remove(node)
	S=EDS.descendents[0]
	S.ancestor=EDS.ancestor
	for n,i in enumerate(EDS.ancestor.descendents):
		if i==EDS:
			EDS.ancestor.descendents[n]=S
		else:
			pass
	del EDS


def addto_root(S,node):
	EDS=Node()
	EDS.descendents.append(node)
	if S.root is True:
		S.ancestor=EDS
		EDS.root=True
		S.root=False
	else:
		raise TreeError("This is not the root node !!")
	EDS.descendents.append(S)
	EDS.ancestor='ROOT'
	return EDS

def removefrom_root(EDS,node):
	EDS.descendents.remove(node)
	if EDS.root is True:
		EDS.descendents[0].ancestor=None
		EDS.descendents[0].root=True
		EDS.root=False
	else:
		raise TreeError("This is not the root node !!")
	del EDS


		



#WED=add_sister(root_node,node)
#	write_topology(ee,s=[])
#	remove_sister(WED,node)
#	write_topology(ee,s=[])
'''
def listgen(root_node,EE):
	if root_node.descendents is not 0:
		[listgen(i,EE) for i in root_node.descendents]
		print(root_node.name)
'''

def listgen(root_node,EE):
	if root_node.descendents is not 0:
		[listgen(i,EE) for i in root_node.descendents]
		EE.append(root_node)


def tiplabel_changer(root_node,p):
	if root_node.ancestor is not None:
		[tiplabel_changer(i,p) for i in root_node.descendents]
		if 'Internal' in root_node.name:
			pass
		else:
			root_node.name=p[root_node.name]
			print(str(root_node.name)+'\t'+"{:}".format(root_node.branch_length))

# test code
def write_tree(root_node,n=0,s=[]):
	if n>0:
		s.append(",")
	if len(root_node.descendents) is not 0:
		s.append("(")
		[write_tree(i,n,s) for n,i in enumerate(root_node.descendents)]
		s.append(")")
		if root_node.branch_length is not None:
			s.append(':'+"{:}".format(root_node.branch_length))
		else:
			pass
	else:
		s.append(root_node.name+':'+"{:}".format(root_node.branch_length))
	return ''.join(s)+';'




@readfiles
def build_topologies(treemodels):
	for k in treemodels:
		ee=read_tree(k)
		EE=[]
		listgen(ee,EE)
		return ee,EE

XKCD_RGBcolors=[xkcd_rgb[i] for n,i in enumerate(xkcd_rgb)]

def build_trees(treemodels,Backbone_target_sequences,ser,rootinc):
	ee,EE=build_topologies(treemodels)
	treecolors=random.choice(XKCD_RGBcolors,len(EE),replace=False)
	Alltrees={}
	for jam in Backbone_target_sequences:
		j=Node(str(aliasfindinator(Backbone_target_sequences[jam][0],ser)))
		print(j.name)
		o=open(str(Backbone_target_sequences[jam][0])+'_'+'Treetopologynewick.nwk','w')
		treedetails=open(str(Backbone_target_sequences[jam][0])+'_'+'Topologymodeldetails.txt','w')
		Alltrees[str(jam)]=[]
		for n,i in enumerate(EE):
			if i.root is False:
				Alltrees[str(jam)].append(i)
				treedetails.write('Tree_'+str(n+1)+'\t'+str(i.name)+'\t'+str(treecolors[n])+'\n')
				WED=add_sister(i,j)
				o.write(write_topology(ee,s=[])+'\n')
				remove_sister(WED,j)
			elif i.root and rootinc:
				Alltrees[str(jam)].append(i)
				treedetails.write('Tree_'+str(n+1)+'\t'+str(i.name)+'\t'+str(treecolors[n])+'\n')
				WED=addto_root(i,j)
				o.write(write_topology(WED,s=[])+'\n')
				removefrom_root(WED,j)
			else:
				pass
		o.close()
		treedetails.close()
	return Alltrees



def prunenbuild_trees(ser):
	Alltrees={}
	for jam in Backbone_target_sequences:
		j=Node(str(aliasfindinator(Backbone_target_sequences[jam][0],ser)))
		print(j.name)
		o=open(str(Backbone_target_sequences[jam][0])+'_'+'Treetopologynewick.nwk','w')
		treedetails=open(str(Backbone_target_sequences[jam][0])+'_'+'Topologymodeldetails.txt','w')
		Alltrees[str(jam)]=[]
		for n,i in enumerate(EE):
			if i.root is False:
				Alltrees[str(jam)].append(i)
				treedetails.write('Tree_'+str(n+1)+'\t'+str(i.name)+'\t'+str(treecolors[n])+'\n')
				WED=add_sister(i,j)
				o.write(write_topology(ee,s=[])+'\n')
				remove_sister(WED,j)
			elif i.root and rootinc:
				Alltrees[str(jam)].append(i)
				treedetails.write('Tree_'+str(n+1)+'\t'+str(i.name)+'\t'+str(treecolors[n])+'\n')
				WED=addto_root(i,j)
				o.write(write_topology(WED,s=[])+'\n')
				removefrom_root(WED,j)
			else:
				pass
		o.close()
		treedetails.close()
	return Alltrees



def tree_topology_creator(ee,EE,j,o,rootinc):
	for n,i in enumerate(EE):
		if i.root is False:
			WED=add_sister(i,j)
			o.write(write_topology(ee,s=[])+'\n')
			remove_sister(WED,j)
		elif i.root and rootinc:
			WED=addto_root(i,j)
			o.write(write_topology(WED,s=[])+'\n')
			removefrom_root(WED,j)
		else:
			pass
	o.close()
		




def tree_topology_creator(ee,EE,j,o,rootinc):
	for n,i in enumerate(EE):
		if i.root is False:
			WED=add_sister(i,j)
			o.write(write_topology(ee,s=[])+'\n')
			remove_sister(WED,j)
		elif i.root and rootinc:
			WED=addto_root(i,j)
			o.write(write_topology(WED,s=[])+'\n')
			removefrom_root(WED,j)
		else:
			pass
	o.close()

def prune_node(node):
	intnode,intnode_index=node.ancestor,node.ancestor.descendents.index(node)
	intnode.descendents.remove(node)
	if intnode.ancestor=='ROOT':
		pass
	else:
		for n,i in enumerate(intnode.descendents):
			i.ancestor=intnode.ancestor
			intnode.ancestor.descendents[intnode.ancestor.descendents.index(intnode)]=i
	return(intnode,intnode_index)

	


def rewire_nodes(intnode,intnode_index,node):
	if intnode.ancestor=='ROOT':
		pass
		node.ancestor=intnode
	else:
		for n,i in enumerate(intnode.descendents):
			if i==intnode.ancestor.descendents[intnode.ancestor.descendents.index(i)]:
				intnode.ancestor.descendents[intnode.ancestor.descendents.index(i)]=intnode
				i.ancestor=intnode
		node.ancestor=intnode
	intnode.descendents.insert(intnode_index,node)
	
