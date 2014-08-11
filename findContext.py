#-*- coding: <utf-16> -*- 
import codecs
import os
import sys
import string
import operator
import math
import collections
import numpy as np
import networkx as nx
import ctypes
import contextlib
import itertools
import time

import multiprocessing as mp

useMultiprocessing = True

cpu = mp.cpu_count()
if cpu == 1:
	howManyCoresToUse = 1
elif cpu == 2:
	howManyCoresToUse = 2
else:
	howManyCoresToUse = cpu-1

#howManyCoresToUse = 10

#--------------------------------------------------------------------------------------------------#

print time.strftime('\nlog_%Y-%m-%d_%H-%M-%S' ,time.localtime())

print '%d CPU cores available, multiprocessing using %d child processes' % (cpu, howManyCoresToUse)

def get_from_multiprocessing(f, args=None):
	pool = mp.Pool(processes = howManyCoresToUse)
	if args:
		result = pool.apply_async(f, args)
	else:
		result = pool.apply_async(f)
	pool.close()
	return result.get()

#useMultiprocessing_response = raw_input('use multiprocessing? (Y/n) ')
#if useMultiprocessing_response.lower() != 'n':
#	useMultiprocessing = True
#	print 'multiprocessing enabled'
#else:
#	useMultiprocessing = False
#	print 'multiprocessing disabled'
#if useMultiprocessing:
#	howManyCoresToUse_response = raw_input('how many child processes? (there are %d CPU cores on this machine) ' % (mp.cpu_count()))
#	howManyCoresToUse = int(howManyCoresToUse_response)

import re
def sep(s, thou=",", dec="."):
	integer, decimal = s.split(".")
	integer = re.sub(r"\B(?=(?:\d{3})+$)", thou, integer)
	return integer + dec + decimal

#-----------------------------------------------------------------------#
#									#
#	This program takes a trigram file and a word list		#	
#	and creates a file with lists of most similar words.		#
#	John Goldsmith and Wang Xiuli 2012.		Jackson Lee 2014		#
#									#
#-----------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#	helper functions
#---------------------------------------------------------------------------#
# 
#

# 1 May 2014: Changed contexts to a dict whose value is dict, not list. JG.



#---------------------------------------------------------------------------#
#	Variables to be changed by user
#---------------------------------------------------------------------------#
LatexFlag = True
PrintEigenvectorsFlag = True
unicodeFlag = False 
FileEncoding =  "ascii"

shortfilename 		= "english-brown"
outshortfilename 	= "english-brown"
languagename 		= "english"

datafolder			= "../../data/"

ngramfolder   		= datafolder + languagename + "/ngrams/"
outfolder	 		= datafolder + languagename + "/neighbors/"

 

NumberOfEigenvectors 		= 11

NumberOfNeighborsKeyboard = raw_input('\nhow many neighbors? (default=9) ')
if not NumberOfNeighborsKeyboard:
    NumberOfNeighbors = 9
else:
    NumberOfNeighbors = int(NumberOfNeighborsKeyboard)

NumberOfWordsForAnalysis_response = raw_input('number of words for analysis? (default=1000) ') #4000
if not NumberOfWordsForAnalysis_response.isdigit():
	NumberOfWordsForAnalysis = 1000 # default value
else:
	NumberOfWordsForAnalysis = int(NumberOfWordsForAnalysis_response)

punctuation 		= " $/+.,;:?!()\"[]"





#---------------------------------------------------------------------------#
#	File names
#---------------------------------------------------------------------------#

infileBigramsname 	= ngramfolder + shortfilename	+ "_bigrams.txt" 
infileTrigramsname 	= ngramfolder + shortfilename	+ "_trigrams.txt"
infileWordsname 	= ngramfolder   + shortfilename	+ "_words.txt" 
outfilenameEigenvectors	= outfolder	 + outshortfilename + "_contexts_eigenvectors" + ".txt"
outfilenameNeighbors	= outfolder	 + outshortfilename + "_" + str(NumberOfWordsForAnalysis) + "_" + str(NumberOfNeighbors) + "_nearest_neighbors_contexts.txt"
outfilenameNeighborsGraphml	= outfolder	 + outshortfilename + "_" + str(NumberOfWordsForAnalysis) + "_nearest_neighbors_contexts.graphml"
outfilenameLatex 	= outfolder	 + outshortfilename + "_latex_contexts.tex"
outfilenameContexts	 = outfolder	 + outshortfilename + "_contexts.txt"
outCSVname = outfilenameNeighbors + '.temp.csv'

if os.path.isfile(outfilenameNeighbors):
	continutKeyboard = raw_input('The file %s already exists.\nAre you sure you want to override it and re-do all the computation? [N/y] ' % (outfilenameNeighbors))
	if continutKeyboard.lower() != 'y':
		print '\nprogram terminated\n'
		sys.exit()

print "\n\nI am looking for: ", infileTrigramsname
print "Number of words that will be analyzed:", NumberOfWordsForAnalysis

#---------------------------------------------------------------------------#
#	Important data structures
#---------------------------------------------------------------------------#

mywords = collections.OrderedDict()

#analyzedwordlist 	= list() # these are the words that will be analyzed

analyzedworddict   	= collections.OrderedDict() # key is word, value is its index in analyzedwordlist

closestNeighbors 	= collections.OrderedDict() #a dict whose values are lists; the lists are the closest words to the key.

contexts 		= collections.defaultdict(set) # key is word, value is a dict of contexts (used to be a list of contexts)

coordinates 		= dict()

#CountOfSharedContexts   is a V by V matrix, where v is the index number of one of the words to be analyzed.

Diameter = dict()

from_word_to_context = collections.defaultdict(set) # this dict takes a word as key, and returns a set as value; values of set are the contexts in which the word appears.

HeavilyWeightedContexts = dict() # key is a word w1, value is a dict called WeightedContexts. Key of WeightedContexts is a context, value of WeightedContexts is the number of words within a ball around w1 that share that context with w1.


wordsdistance 		= dict() # key is a word, word1,  being analyzed, value is a pair of word-index-number and euclidean distance (word2, distance). This will be sorted to get the nearest neighbors to word1.


#---------------------------------------------------------------------------#
#	Variables
#---------------------------------------------------------------------------#

linecount 		= 0


 
#---------------------------------------------------------------------------#
#	Normalize function
#---------------------------------------------------------------------------#
def Normalize(NumberOfWordsForAnalysis, CountOfSharedContexts):
	diameterDict = dict()
	for w1 in range(NumberOfWordsForAnalysis):
		diameterDict[w1] = 0
		for w2 in range(NumberOfWordsForAnalysis):
			if w1 == w2:
				continue
			diameterDict[w1] += CountOfSharedContexts[w1,w2]
		if diameterDict[w1] == 0:
		 	diameterDict[w1] = 1
	return diameterDict

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
# this function calculates  contexts shared by two words 
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#def CountSharedContextsFunction(word1, word2, contexts,from_word_to_context):
#	return len(set(from_word_to_context[word1]).intersection(set(from_word_to_context[word2])))

#def SharedContextsFunction(word1, word2,  from_word_to_context):
#	return set(from_word_to_context[word1]).intersection(set(from_word_to_context[word2]))

def FindListOfSharedContexts(word1, word2,  from_word_to_context):
#def SharedContextsFunction(word1, word2,  from_word_to_context):
	returnedcontexts = ()
	for context in from_word_to_context[word1]:
		if context in from_word_to_context[word2]:
			returnedcontexts[context] = 1
	return returnedcontexts

def GetNumberOfSharedWords(context1, context2, from_context_to_wordDict):

#	word1contextSet = from_word_to_context[word1] # set of contexts 
#	word2contextSet = from_word_to_context[word2] # set of contexts

	return len(from_context_to_wordDict[context1] & from_context_to_wordDict[context2])


def GetNumberOfSharedContexts(word1, word2, from_word_to_context):

#	word1contextSet = from_word_to_context[word1] # set of contexts 
#	word2contextSet = from_word_to_context[word2] # set of contexts

	return len(from_word_to_context[word1] & from_word_to_context[word2])

#	count = 0
#	if len(from_word_to_context[word1]) < len(from_word_to_context[word2]):
#		for context in from_word_to_context[word1]:
#			if context in from_word_to_context[word2]:
#				count += 1
#		return count	
#	else:
#		for context in from_word_to_context[word2]:
#			if context in from_word_to_context[word1]:
#				count += 1
#		return count		



	

def WeightedSharedContextsFunction(word1, word2, from_word_to_context,HeavilyWeightedContexts, weight):
	count = 0
	for context in from_word_to_context[word1]:
		if context in from_word_to_context[word2]:
			if context in HeavilyWeightedContexts[word1]:
				count += weight
			else:
				count += 1
	return count			
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
#	Open files for reading and writing
#---------------------------------------------------------------------------#

if unicodeFlag:
	trigramfile 		=codecs.open(infileTrigramsname, encoding = FileEncoding)
	wordfile 		=codecs.open(infileWordsname, encoding = FileEncoding)
	if PrintEigenvectorsFlag:
		outfileEigenvectors = codecs.open (outfilename1, "w",encoding = FileEncoding)
	outfileNeighbors	= codecs.open (outfileneighborsname, "w",encoding = FileEncoding)

 
	
else:
	if PrintEigenvectorsFlag:
		outfileEigenvectors = open (outfilenameEigenvectors, "w")
	outfileNeighbors	= open (outfilenameNeighbors, "w")
	outfileLatex 		= open (outfilenameLatex, "w")
	outfileContexts 	= open (outfilenameContexts, "w")

	wordfile		= open(infileWordsname)
	trigramfile 		= open(infileTrigramsname)
	bigramfile 		= open(infileBigramsname)
print "Language is", languagename, ". File name:", shortfilename, ". Number of contexts", NumberOfWordsForAnalysis, "."

if PrintEigenvectorsFlag:
	print >>outfileEigenvectors,"#", \
			languagename, "\n#", \
			shortfilename,"\n#", \
			"Number of words analyzed", NumberOfWordsForAnalysis,"\n#", \
			"Number of neighbors identified", NumberOfNeighbors, "\n#","\n#" 

print >>outfileNeighbors, "# language:", \
		languagename, "\n# corpus:",\
		shortfilename, "\n#",\
		"Number of words analyzed", NumberOfWordsForAnalysis,"\n#", \
		"Number of neighbors identified", NumberOfNeighbors,"\n"

print >>outfileContexts, "#  The number with each context is the number of distinct words found in that context.\n#" 

 
#---------------------------------------------------------------------------#
#	Read trigram file
#---------------------------------------------------------------------------#


for line in wordfile:
	pieces = line.split()
	if pieces[0] == "#":
		continue
	mywords[pieces[0]] = int(pieces[1])		 
print "1. Word file is ", infileWordsname, '\t corpus has', len(mywords), 'words'
wordfile.close()

#analyzedwordlist = sorted(mywords,key=mywords.__getitem__,reverse=True)
#analyzedwordlist[NumberOfWordsForAnalysis:] = []

#if NumberOfWordsForAnalysis > len(analyzedwordlist) :
#	print "We will reduce NumberOfWordsForAnalysis to" , len(analyzedwordlist)
#	NumberOfWordsForAnalysis = len(analyzedwordlist)


if NumberOfWordsForAnalysis > len(mywords):
	NumberOfWordsForAnalysis = len(mywords)
	print 'number of words for analysis reduced to', NumberOfWordsForAnalysis

analyzedwordlist = mywords.keys()[ : NumberOfWordsForAnalysis]
analyzedwordset = set(analyzedwordlist)

for i in range(NumberOfWordsForAnalysis):
	analyzedworddict[analyzedwordlist[i]] = i


 
print "2. Reading in trigram file."
for line in trigramfile:
	if line.startswith('#'):
		continue
	thesewords = line.split()
	linecount += 1

	thisword = thesewords[1]
	if thisword in analyzedwordset:
		context = thesewords[0] + " __ " +  thesewords[2]
		wordno = analyzedwordlist.index(thisword)
		contexts[context].add(wordno)
#		from_word_to_context[wordno].add(context)

	#Left trigrams
	thisword = thesewords[0]
	if thisword in analyzedwordset:
		context = " __ " + thesewords[1] + " " + thesewords[2]
		wordno = analyzedwordlist.index(thisword)
		contexts[context].add(wordno)
#		from_word_to_context[wordno].add(context)
			
	#Right trigrams
	thisword = thesewords[2]
	if thisword   in analyzedwordset:	
		context = thesewords[0] + " " + thesewords[1] + " __ "
		wordno = analyzedwordlist.index(thisword)
		contexts[context].add(wordno)
#		from_word_to_context[wordno].add(context)

 
#---------------------------------------------------------------------------#
#	Read bigram file
#---------------------------------------------------------------------------#
if True: 
	print "...Reading in bigram file."
	for line in bigramfile:
		linecount += 1
		thesewords = line.split()
		if thesewords[0] == "#":
			continue
	  
		thisword = thesewords[1]
		if thisword in analyzedwordset:
			context = thesewords[0] + " __ " 
			wordno = analyzedwordlist.index(thisword)
			contexts[context].add(wordno)
#			from_word_to_context[thisword][context]=1

		 
		thisword = thesewords[0]
		if thisword in analyzedwordset:
			context = "__ " + thesewords[1]  
			wordno = analyzedwordlist.index(thisword)
			contexts[context].add(wordno)
#			from_word_to_context[thisword][context]=1
		 
 
NumberOfContextsForAnalysis = 1000

contextList = [context for (context, wordCount) in sorted(contexts.items(), key=lambda x: len(x[1]), reverse=True)][: NumberOfContextsForAnalysis]

 
print "3. End of words and counts."
 

#---------------------------------------------------------------------------#
#	Count context features shared by words
#---------------------------------------------------------------------------#

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "4. Counting context features shared by words...",



# ###### not using multiprocessing ########

#print 'not using multiprocessing...'
#CountOfSharedContexts1 = np.zeros( (NumberOfWordsForAnalysis,NumberOfWordsForAnalysis) )
#for wordno1 in range(NumberOfWordsForAnalysis):
#	word1 = analyzedwordlist[wordno1]
#	for wordno2 in  range(wordno1+1, len(analyzedwordlist)):		 
#		word2 = analyzedwordlist[wordno2]
#		CountOfSharedContexts1[wordno1,wordno2] =   GetNumberOfSharedContexts(word1, word2,  from_word_to_context) 
#		CountOfSharedContexts1[wordno2,wordno1] =   CountOfSharedContexts1[wordno1,wordno2] 
#CountOfSharedContexts = CountOfSharedContexts1

# ###### END OF not using multiprocessing ########

def counting_words_shared_by_contexts():
	CountOfSharedWords1 = np.zeros( (NumberOfContextsForAnalysis,NumberOfContextsForAnalysis) )
	for contextno1 in range(NumberOfContextsForAnalysis):
		for contextno2 in  range(contextno1+1, NumberOfContextsForAnalysis):		 
			CountOfSharedWords1[contextno1,contextno2] =   GetNumberOfSharedWords(contextList[contextno1], contextList[contextno2],  contexts)
	return CountOfSharedWords1

CountOfSharedWords = get_from_multiprocessing(counting_words_shared_by_contexts)
CountOfSharedWords = CountOfSharedWords + CountOfSharedWords.T


def counting_context_features():
	CountOfSharedContexts1 = np.zeros( (NumberOfWordsForAnalysis,NumberOfWordsForAnalysis) )
	for wordno1 in range(NumberOfWordsForAnalysis):
		for wordno2 in  range(wordno1+1, NumberOfWordsForAnalysis):		 
			CountOfSharedContexts1[wordno1,wordno2] =   GetNumberOfSharedContexts(wordno1, wordno2,  from_word_to_context)
	return CountOfSharedContexts1

#CountOfSharedContexts = get_from_multiprocessing(counting_context_features)
#CountOfSharedContexts = CountOfSharedContexts + CountOfSharedContexts.T

del from_word_to_context

print 'done.'

#---------------------------------------------------------------------------#
#	Normalize function call
#---------------------------------------------------------------------------#


print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "5. Normalizing nearness measurements....",


Diameter = get_from_multiprocessing(Normalize, [NumberOfContextsForAnalysis, CountOfSharedWords])

print "\t Done." 


#---------------------------------------------------------------------------#
#	Incidence graph
#---------------------------------------------------------------------------#

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "6. We compute the incidence graph....",

def compute_incidence_graph():
	incidencegraph= np.zeros( (NumberOfContextsForAnalysis,NumberOfContextsForAnalysis) )

	for (c1, c2) in itertools.product(range(NumberOfContextsForAnalysis), repeat=2):
		if c1 == c2:
			incidencegraph[c1,c1] = Diameter[c1]
		else:
			incidencegraph[c1,c2] = CountOfSharedWords[c1,c2]	
	return incidencegraph

incidencegraph = get_from_multiprocessing(compute_incidence_graph)

del CountOfSharedWords

print "Done."

 		
#---------------------------------------------------------------------------#
#	Compute raw/naive neighbors.
#---------------------------------------------------------------------------#

 
 
#---------------------------------------------------------------------------#
#	Normalize the laplacian.
#---------------------------------------------------------------------------#

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "7. We normalize the laplacian....",

def compute_laplacian():
	mylaplacian = np.zeros((NumberOfContextsForAnalysis,NumberOfContextsForAnalysis) )

	for (i, j) in itertools.product(range(NumberOfContextsForAnalysis), repeat=2):
		if i == j:
			mylaplacian[i,j] = 1
		else:
			if incidencegraph[i,j] == 0:
				mylaplacian[i,j]=0
			else:
				mylaplacian[i,j] = -1 * incidencegraph[i,j]/ math.sqrt ( Diameter[i] * Diameter[j] )
	return mylaplacian

mylaplacian = get_from_multiprocessing(compute_laplacian)

del incidencegraph
del Diameter

print "Done." 



#---------------------------------------------------------------------------#
#	Compute eigenvectors.
#---------------------------------------------------------------------------#

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "8. Compute eigenvectors...",

myeigenvalues, myeigenvectors = get_from_multiprocessing(np.linalg.eigh, [mylaplacian])

print "done."

del mylaplacian
  
#---------------------------------------------------------------------------#
#	Generate latex output.
#---------------------------------------------------------------------------#
if LatexFlag:
	#Latex output
	print >>outfileLatex, "%",  infileWordsname
	print >>outfileLatex, "\\documentclass{article}" 
	print >>outfileLatex, "\\usepackage{booktabs}" 
	print >>outfileLatex, "\\begin{document}" 

data = dict() # key is eigennumber, value is list of triples: (index, word, eigen^{th} coordinate) sorted by increasing coordinate
print ("9. Printing contexts to latex file.")
formatstr = '%20s   %15s %10.3f'
headerformatstr = '%20s  %15s %10.3f %10s'
NumberOfWordsToDisplayForEachEigenvector = 20
	

	
		 
if PrintEigenvectorsFlag:

	for eigenno in range(NumberOfEigenvectors):
		print >>outfileEigenvectors
		print >>outfileEigenvectors,headerformatstr %("Eigenvector number", eigenno, myeigenvalues[eigenno], "word" )
		print >>outfileEigenvectors,"_"*50 

		eigenlist=list()		
		for wordno in range (NumberOfWordsForAnalysis):		 
			eigenlist.append( (wordno,myeigenvectors[wordno, eigenno]) )			
		eigenlist.sort(key=lambda x:x[1])			

		for wordno in range(NumberOfWordsForAnalysis):	
			word = analyzedwordlist[eigenlist[wordno][0]]
			coord =  eigenlist[wordno][1]		
			print >>outfileEigenvectors, formatstr %(eigenno, word, eigenlist[wordno][1])


 

if LatexFlag:
	for eigenno in range(NumberOfEigenvectors):
		eigenlist=list()	
		data = list()
		for wordno in range (NumberOfWordsForAnalysis):		 
			eigenlist.append( (wordno,myeigenvectors[wordno, eigenno]) )			
		eigenlist.sort(key=lambda x:x[1])			
		print >>outfileLatex			 
		print >>outfileLatex, "Eigenvector number", eigenno, "\n" 
		print >>outfileLatex, "\\begin{tabular}{lll}\\toprule"
		print >>outfileLatex, " & word & coordinate \\\\ \\midrule "

		for i in range(NumberOfWordsForAnalysis):			 
			word = analyzedwordlist[eigenlist[i][0]]
			coord =  eigenlist[i][1]
			if i < NumberOfWordsToDisplayForEachEigenvector or i > NumberOfWordsForAnalysis - NumberOfWordsToDisplayForEachEigenvector:
				data.append((i, word , coord ))
		for (i, word, coord) in data:
			if word == "&":
				word = "\&" 
			print >>outfileLatex,  "%5d & %10s &  %10.3f \\\\" % (i, word, coord) 

		print >>outfileLatex, "\\bottomrule \n \\end{tabular}", "\n\n"
		print >>outfileLatex, "\\newpage" 
print >>outfileLatex, "\\end{document}" 
#---------------------------------------------------------------------------#
#	Finding coordinates in space of low dimensionality
#---------------------------------------------------------------------------#

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "10. Finding coordinates in space of low dimensionality."

def compute_coordinates():
	coordinates = dict()
	for contextno in range(NumberOfContextsForAnalysis):
		coordinates[contextno]= list() 
		for eigenno in range(NumberOfEigenvectors):
			coordinates[contextno].append( myeigenvectors[ contextno, eigenno ] )
	return coordinates

coordinates = get_from_multiprocessing(compute_coordinates)

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print '       coordinates computed. Now computing distances between words...',

del myeigenvectors


def compute_contexts_distance():
	contextsdistance = dict()

	for contextno in range(NumberOfContextsForAnalysis):
		contextsdistance[contextno] = dict()

	for (contextno1, contextno2) in itertools.combinations(range(NumberOfContextsForAnalysis), 2):
		distance = 0
		for coordno in range(NumberOfEigenvectors):
			x = coordinates[contextno1][coordno] - coordinates[contextno2][coordno]
			distance += abs(x ** 3)
#		wordsdistance[(wordno1, wordno2)] = distance
		contextsdistance[contextno1][contextno2] = distance
		contextsdistance[contextno2][contextno1] = distance

	return contextsdistance

contextsdistance = get_from_multiprocessing(compute_contexts_distance)

#print '\nworddistance:\n'
#for (i,j) in itertools.product(range(10),repeat=2):
#	print '{0:10s} {1:10s} {2:5.10f}'.format(analyzedwordlist[i], analyzedwordlist[j], wordsdistance[i][j])

print 'Done.'

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print '       Writing the word distance dict to CSV...',

outCSVfile = open(outCSVname, 'w')

for (contextno1, contextno2) in itertools.combinations(range(NumberOfContextsForAnalysis), 2):
	outCSVfile.write('%d,%d,%f\n' % (contextno1, contextno2, contextsdistance[contextno1][contextno2]))

outCSVfile.close()

print 'Done.'

#---------------------------------------------------------------------------#
#	 Finding closest neighbors on the manifold's approximation
#---------------------------------------------------------------------------#

print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime())
print "11. Finding closest neighbors on the manifold('s approximation)."


def compute_graph():
	mygraph = nx.Graph()
	for (contextno, context1)  in enumerate(contextList):
		mygraph.add_node(context1, contextindex=contextno)
	return mygraph

mygraph = get_from_multiprocessing(compute_graph)


print '      graph initialized. Computing nearest neighbors now...'
print time.strftime('        log_%Y-%m-%d_%H-%M-%S' ,time.localtime()),


def compute_closest_neighbors():
	closestNeighbors = dict()

#	wordsdistanceSorted = sorted(wordsdistance.items(), key=lambda x:x[1])

	for (contextno1, context1) in enumerate(contextList):
#		closestNeighbors[wordno1] = list() # list of word indices
#		neighborWordNumberList = [wordno for (wordno, distance) in sorted(wordsdistance[wordno1].items(), key=lambda x:x[1])]

		neighborContextNumberList = [contextno2 for (contextno2, distance) in sorted(contextsdistance[contextno1].items(), key=lambda x:x[1])]

#################

#		neighborWordNumberList = list()
#		neighborDistanceList = list()

#		neighborCount = 0
#		currentMaxDistance = 0

#		for ((wordno1_inPair, wordno2_inPair), distance) in wordsdistanceSorted:
#			if wordno1 not in (wordno1_inPair, wordno2_inPair):
#				continue

#			if currentMaxDistance and (neighborCount > NumberOfNeighbors) and \
#				(distance > currentMaxDistance):
#				break

#			neighborCount += 1
#			if wordno1 == wordno1_inPair:
#				wordno2 = wordno2_inPair
#			else:
#				wordno2 = wordno1_inPair

#			if distance > currentMaxDistance:
#				currentMaxDistance = distance

#			neighborDistanceList.append(distance)
#			neighborWordNumberList.append(wordno2)

#################

		neighborContextNumberList = neighborContextNumberList[: NumberOfNeighbors]

		closestNeighbors[contextno1] = neighborContextNumberList

		for (idx, contextno2) in enumerate(neighborContextNumberList):
			mygraph.add_edge(context1, contextList[contextno2], rank=idx+1)


#		for (idx, wordno2) in enumerate(neighborWordNumberList):		
#			if wordno2 ==  wordno1:			 
#				continue			
#			if idx > NumberOfNeighbors:
#				break
#			word2 = analyzedwordlist[wordno2]

#			closestNeighbors[wordno1].append(wordno2)
#			mygraph.add_edge(word1, word2, rank=idx)		

	return closestNeighbors

#closestNeighbors = get_from_multiprocessing(compute_closest_neighbors)
closestNeighbors = compute_closest_neighbors()

for (contextno, context) in enumerate(contextList):
	print >>outfileNeighbors, context+'\t'+ '\t'.join([contextList[idx] for idx in closestNeighbors[contextno]])

#for (word, neighborList) in sorted(closestNeighbors.items(), key=lambda x: analyzedworddict[x[0]]):
#	print >>outfileNeighbors, word, ' '.join(neighborList)

print 'done.'

nx.write_graphml(mygraph,outfilenameNeighborsGraphml)
 


#---------------------------------------------------------------------------#
#	 Reweight contexts based on word neighborhoods: change, March 2014. JG. 
#---------------------------------------------------------------------------#
NewContexts = dict()
if (False):
	print "12. Reweighting contexts based on word neighborhoods",

	print >>outfileContexts, "For each word, we print how often it shares a context with one of its nearest neighbors (with a threshold set on how many words shared that context).\n\n"
	for word1 in analyzedwordlist:	 
		print >>outfileContexts, "\n\n---------------\n", word1
		HeavilyWeightedContexts[word1] = dict()
		for word2 in closestNeighbors[word1]:
			for context in FindListOfSharedContexts(word1, word2,  from_word_to_context):
				if context not in HeavilyWeightedContexts[word1]:
					HeavilyWeightedContexts[word1][context] = list()
				HeavilyWeightedContexts[word1][context].append(word2)
		for context in HeavilyWeightedContexts[word1]:
			if  len(HeavilyWeightedContexts[word1][context]) > 5:	
				print >>outfileContexts, "\n\n\t", context		
				for word in HeavilyWeightedContexts[word1][context]:
					print >>outfileContexts, "%12s" % word, 
				print >>outfileContexts				
				 
	print "...Done" 	


if False:
	print "13. Counting context features shared by words...",

	CountOfSharedContexts = zeros( (NumberOfWordsForAnalysis,NumberOfWordsForAnalysis) )
	 
	weight = 10
	for word1 in analyzedwordlist:	 
		#print  word1,
		wordno1 = analyzedworddict[word1]
		for word2 in analyzedwordlist:
			wordno2 = analyzedworddict[word2]
			CountOfSharedContexts[wordno1,wordno2] =   WeightedSharedContextsFunction(word1, word2, from_word_to_context ,HeavilyWeightedContexts, weight)
		
	print "... Done." 

 
		
outfileNeighbors.close()
#---------------------------------------------------------------------------#
#	 Print contexts shared by nearby words: not finished
#---------------------------------------------------------------------------#

if False:
	numberperrow= 5
	for word in analyzedwordlist:
		print >>outfileContexts,"\n\n", word
		number = 0
		if len(from_word_to_context[word]) < 200:
			continue
		for context in from_word_to_context[word]:
			if len(contexts[context]) < 20:
				continue
			if number == 0:
				print >>outfileContexts, "\n\t",
			thesecontexts = set(from_word_to_context[word])
			print >>outfileContexts, "%3d  %-20s " %(len(contexts[context]), context, ),
			number += 1
			if number == numberperrow:
				number = 0
	 

 
print "Exiting successfully."

#os.popen("latex " + outfilenameLatex ) 

if PrintEigenvectorsFlag:
	outfileEigenvectors.close()
outfileNeighbors.close()
 
print time.strftime('log_%Y-%m-%d_%H-%M-%S' ,time.localtime())



