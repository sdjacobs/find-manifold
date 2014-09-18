#-----------------------------------------------------------------------#
#
#	This program takes n-gram files and a word list	
#	and creates a file with lists of most similar words.
#	John Goldsmith and Wang Xiuli 2012.
#	Jackson Lee 2014
#
#-----------------------------------------------------------------------#

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
import itertools
import time
import subprocess
from contextlib import closing

import multiprocessing as mp

#--------------------------------------------------------------------------------------------------#

timeFormat = '      current time: %Y-%m-%d %H:%M:%S'

beginTime = time.localtime()
print
print time.strftime(timeFormat , beginTime)

#--------------------------------------------------------------------------------------------------#

# multiprocessing settings and associated functions  
   
cpu = mp.cpu_count()  
if cpu > 2:  
    howManyCoresToUse = cpu-1  
else:  
   howManyCoresToUse = cpu  
   
print 'multiprocessing spawning %d child processes' % (howManyCoresToUse)  


#def get_from_multiprocessing(f, args=None):
#	global howManyCoresToUse
#	pool = mp.Pool(processes = howManyCoresToUse)
#	if args:
#		result = pool.apply_async(f, args)
#	else:
#		result = pool.apply_async(f)
#	return result.get()

def compute_chunksize(N, howManyCoresToUse):
	floatNum = float(N/howManyCoresToUse) / 10
	return int(math.ceil(floatNum) * 10)

#---------------------------------------------------------------------------#
#	Variables to be changed by user
#---------------------------------------------------------------------------#
LatexFlag = False
PrintEigenvectorsFlag = False
unicodeFlag = False 
FileEncoding =  "ascii"

shortfilename 		= "english-brown"
outshortfilename 	= "english-brown"
languagename 		= "english"

datafolder			= "../../data/"

ngramfolder   		= datafolder + languagename + "/ngrams/"
outfolder	 		= datafolder + languagename + "/neighbors/"
wordcontextfolder	= datafolder + languagename + "/word_contexts/"

NumberOfEigenvectors 		= 11

punctuation 		= " $/+.,;:?!()\"[]"

if '-m' in sys.argv:
	runningOnMidway = True
	sys.argv.remove('-m')
else:
	runningOnMidway = False


#def get_NumberOfNeighbors():
#	NumberOfNeighborsKeyboard = raw_input('\nhow many neighbors? (default=9) ')
#	if not NumberOfNeighborsKeyboard:
#		NumberOfNeighbors = 9
#	else:
#		NumberOfNeighbors = int(NumberOfNeighborsKeyboard)

#def get_NumberOfWordsForAnalysis():
#	NumberOfWordsForAnalysis_response = raw_input('number of words for analysis? (default=1000) ') #4000
#	if not NumberOfWordsForAnalysis_response.isdigit():
#		NumberOfWordsForAnalysis = 1000 # default value
#	else:
#		NumberOfWordsForAnalysis = int(NumberOfWordsForAnalysis_response)

try:
	if len(sys.argv) < 2:
		NumberOfNeighbors = 9
		NumberOfWordsForAnalysis = 1000
		continueKeyboard = raw_input('default values: NumberOfWordsForAnalysis = 1000, NumberOfNeighbors = 9\nContinue? [N/y] ')
		if continueKeyboard.lower() != 'y':
			print '\nprogram terminated\n'
			sys.exit()
	elif len(sys.argv) == 2:
		NumberOfNeighbors = 9
		NumberOfWordsForAnalysis = int(sys.argv[1])
	else:
		NumberOfNeighbors = int(sys.argv[2])
		NumberOfWordsForAnalysis = int(sys.argv[1])
except:
	print 'usage: python findManifold.py [NumberOfWordsForAnalysis] [NumberOfNeighbors]'
	sys.exit()

#---------------------------------------------------------------------------#
#	File names
#---------------------------------------------------------------------------#

infileBigramsname 	= ngramfolder + shortfilename	+ "_bigrams.txt" 
infileTrigramsname 	= ngramfolder + shortfilename	+ "_trigrams.txt"
infileWordsname 	= ngramfolder   + shortfilename	+ "_words.txt" 
outfilenameEigenvectors	= outfolder	 + outshortfilename + "_words_eigenvectors" + ".txt"
outfilenameNeighbors	= outfolder	 + outshortfilename + "_" + str(NumberOfWordsForAnalysis) + "_" + str(NumberOfNeighbors) + "_nearest_neighbors.txt"
#outfilenameNeighborsGraphml	= outfolder	 + outshortfilename + "_" + str(NumberOfWordsForAnalysis) + "_nearest_neighbors.graphml"
outfilenameLatex 	= outfolder	 + outshortfilename + "_latex.tex"
outfilenameContexts	 = outfolder	 + outshortfilename + "_contexts.txt"
outfilenameFromWordToContexts	 = outfolder	 + outshortfilename + "_" + str(NumberOfWordsForAnalysis)  + "_from-word-to-contexts.txt"
#outCSVname = outfilenameNeighbors + '.temp.csv'

if (not runningOnMidway) and os.path.isfile(outfilenameNeighbors):
	continutKeyboard = raw_input('The file %s already exists.\nAre you sure you want to override it and re-do all the computation? [N/y] ' % (outfilenameNeighbors))
	if continutKeyboard.lower() != 'y':
		print '\nprogram terminated\n'
		sys.exit()

print "\nI am looking for: ", infileTrigramsname
print "Number of words that will be analyzed:", NumberOfWordsForAnalysis
print "Number of neighbors:", NumberOfNeighbors

#---------------------------------------------------------------------------#
#	Important data structures
#---------------------------------------------------------------------------#

mywords = collections.OrderedDict()

#analyzedwordlist 	= list() # these are the words that will be analyzed

#analyzedworddict   	= collections.OrderedDict() # key is word, value is its index in analyzedwordlist

closestNeighbors 	= collections.OrderedDict() #a dict whose values are lists; the lists are the closest words to the key.

contexts 		= dict() # key is word, value is a dict of contexts (used to be a list of contexts)

coordinates 		= dict()

#CountOfSharedContexts   is a V by V matrix, where v is the index number of one of the words to be analyzed.

Diameter = dict()

from_word_to_context = collections.defaultdict(collections.Counter) # this dict takes a word as key, and returns a collections.Counter dict as the value; the value is a dict with (context, frequency count) pairs.
#word_context_graph = collections.defaultdict(nx.DiGraph) # this dict takes a word as key, and returns a networkx Graph as value; the Graph is a word context graph

#HeavilyWeightedContexts = dict() # key is a word w1, value is a dict called WeightedContexts. Key of WeightedContexts is a context, value of WeightedContexts is the number of words within a ball around w1 that share that context with w1.


wordsdistance 		= dict() # key is a word, word1,  being analyzed, value is a pair of word-index-number and euclidean distance (word2, distance). This will be sorted to get the nearest neighbors to word1.



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
# this function calculates  contexts shared by two words 
#---------------------------------------------------------------------------#

#def FindListOfSharedContexts(word1, word2,  from_word_to_context): ## function not used currently
##def SharedContextsFunction(word1, word2,  from_word_to_context):
#	returnedcontexts = ()
#	for context in from_word_to_context[word1]:
#		if context in from_word_to_context[word2]:
#			returnedcontexts[context] = 1
#	return returnedcontexts

def GetNumberOfSharedContexts(word1, word2):
	return len(set(from_word_to_context[word1]) & set(from_word_to_context[word2]))

#def WeightedSharedContextsFunction(word1, word2, from_word_to_context,HeavilyWeightedContexts, weight): # function not used currently
#	count = 0
#	for context in from_word_to_context[word1]:
#		if context in from_word_to_context[word2]:
#			if context in HeavilyWeightedContexts[word1]:
#				count += weight
#			else:
#				count += 1
#	return count			

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
	outfileFromWordToContexts = open(outfilenameFromWordToContexts, 'w')

	wordfile		= open(infileWordsname)
	trigramfile 		= open(infileTrigramsname)
	bigramfile 		= open(infileBigramsname)
print "Language is", languagename, ". File name:", shortfilename, ". Number of words", NumberOfWordsForAnalysis, "."

if PrintEigenvectorsFlag:
	print >>outfileEigenvectors,"#", \
			languagename, "\n#", \
			shortfilename,"\n#", \
			"Number of words analyzed", NumberOfWordsForAnalysis,"\n#", \
			"Number of neighbors identified", NumberOfNeighbors, "\n#","\n#" 

for outfile in [outfileNeighbors, outfileFromWordToContexts]:
	print >>outfile, "# language:", \
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

if NumberOfWordsForAnalysis > len(mywords):
	NumberOfWordsForAnalysis = len(mywords)
	print 'number of words for analysis reduced to', NumberOfWordsForAnalysis

analyzedwordlist = mywords.keys()[ : NumberOfWordsForAnalysis]
analyzedwordset = set(analyzedwordlist)

del mywords

#for i in range(NumberOfWordsForAnalysis):
#	analyzedworddict[analyzedwordlist[i]] = i

print "2. Reading in trigram file."
for line in trigramfile:
	if line.startswith('#'):
		continue
	thesewords = line.split()

	thisword = thesewords[1]
	if thisword in analyzedwordset:
		context = thesewords[0] + " __ " +  thesewords[2]
		wordno = analyzedwordlist.index(thisword)
#		if not context in contexts:
#			contexts[context] = dict()
#		contexts[context][thisword]  = trigram_count
#		from_word_to_context[wordno].add(context)
		from_word_to_context[wordno][context] += 1
#		word_context_graph[wordno].add_edge(thesewords[0], '_')
#		word_context_graph[wordno].add_edge('_', thesewords[2])

	#Left trigrams
	thisword = thesewords[0]
	if thisword in analyzedwordset:
		context = " __ " + thesewords[1] + " " + thesewords[2]
		wordno = analyzedwordlist.index(thisword)
#		if not context in contexts:
#			contexts[context] = dict()
#		contexts[context][thisword] = trigram_count
#		from_word_to_context[wordno].add(context)
		from_word_to_context[wordno][context] += 1
#		word_context_graph[wordno].add_edge(thesewords[1], thesewords[2])
#		word_context_graph[wordno].add_edge('_', thesewords[1])

	#Right trigrams
	thisword = thesewords[2]
	if thisword   in analyzedwordset:	
		context = thesewords[0] + " " + thesewords[1] + " __ "
		wordno = analyzedwordlist.index(thisword)
#		if not context in contexts:
#			contexts[context] = dict()
#		contexts[context][thisword] = trigram_count
#		from_word_to_context[wordno].add(context)
		from_word_to_context[wordno][context] += 1
#		word_context_graph[wordno].add_edge(thesewords[0], thesewords[1])
#		word_context_graph[wordno].add_edge(thesewords[1], '_')

 
#---------------------------------------------------------------------------#
#	Read bigram file
#---------------------------------------------------------------------------#
print "...Reading in bigram file."
for line in bigramfile:
	thesewords = line.split()
	if thesewords[0] == "#":
		continue
  
	thisword = thesewords[1]
	if thisword in analyzedwordset:
		context = thesewords[0] + " __ " 
#			if not context in contexts:
#				contexts[context] = dict()
#			contexts[context][thisword] =1
#			from_word_to_context[wordno].add(context)
		from_word_to_context[wordno][context] += 1
#		word_context_graph[wordno].add_edge(thesewords[0], '_')

	 
	thisword = thesewords[0]
	if thisword in analyzedwordset:
		context = "__ " + thesewords[1]  
#			if not context in contexts:
#				contexts[context] = dict()
#			contexts[context][thisword] = 1
#			from_word_to_context[wordno].add(context)
		from_word_to_context[wordno][context] += 1
#		word_context_graph[wordno].add_edge('_', thesewords[1])

#---------------------------------------------------------------------------#

if False:
	print '...writing in from-word-to-contexts file.',

	for (wordno, word) in enumerate(analyzedwordlist):
		contextOrderedList = [context for (context, freq) in from_word_to_context[wordno].most_common()]
		print>> outfileFromWordToContexts, '%s\t%s' % (word, '\t'.join(contextOrderedList))
	outfileFromWordToContexts.close()

	print '   done.'

#---------------------------------------------------------------------------#

if False:
	print '...writing in word context graphs',

	for (wordno, word) in enumerate(analyzedwordlist):
		if '/' in word:
			continue
		numOfDigits = len(str(NumberOfWordsForAnalysis))
		g = word_context_graph[wordno]
		outGexf = wordcontextfolder + str(wordno).zfill(numOfDigits) + '_' + word + '_' + str(g.number_of_nodes()) + '.gexf'
		nx.write_gexf(g, outGexf)

	print '   done.'

	del word_context_graph

#---------------------------------------------------------------------------#
 
print "3. End of words and counts."


#---------------------------------------------------------------------------#
#	Count context features shared by words
#---------------------------------------------------------------------------#

print time.strftime(timeFormat ,time.localtime())
print "4. Counting context features shared by words...",

datatype = ctypes.c_int


#################################################################
# multiprocessing method 1

def counting_context_features_1core():
    arr = np.zeros((NumberOfWordsForAnalysis, NumberOfWordsForAnalysis))
    for word1 in range(0, NumberOfWordsForAnalysis):
        for word2 in range(word1+1, NumberOfWordsForAnalysis):
            x = GetNumberOfSharedContexts(word1, word2)
            arr[word1, word2] = x
            arr[word2, word1] = x
    return arr

CountOfSharedContexts = counting_context_features_1core()

del from_word_to_context

print 'done.'

#---------------------------------------------------------------------------#
#	Normalize function call
#---------------------------------------------------------------------------#


print time.strftime(timeFormat ,time.localtime())
print "5. Normalizing nearness measurements....",

Diameter = Normalize(NumberOfWordsForAnalysis, CountOfSharedContexts)
#Diameter = get_from_multiprocessing(Normalize, [NumberOfWordsForAnalysis, CountOfSharedContexts])

print "\t Done." 


#---------------------------------------------------------------------------#
#	Incidence graph
#---------------------------------------------------------------------------#

print time.strftime(timeFormat ,time.localtime())
print "6. We compute the incidence graph....",

def compute_incidence_graph():
	incidencegraph= np.zeros( (NumberOfWordsForAnalysis,NumberOfWordsForAnalysis), dtype=np.int32)

	for (w1, w2) in itertools.product(range(NumberOfWordsForAnalysis), repeat=2):
		if w1 == w2:
			incidencegraph[w1,w1] = Diameter[w1]
		else:
			incidencegraph[w1,w2] = CountOfSharedContexts[w1,w2]	
	return incidencegraph

incidencegraph = compute_incidence_graph()
#incidencegraph = get_from_multiprocessing(compute_incidence_graph)

del CountOfSharedContexts

print "Done."

 		
#---------------------------------------------------------------------------#
#	Compute raw/naive neighbors.
#---------------------------------------------------------------------------#

 
 
#---------------------------------------------------------------------------#
#	Normalize the laplacian.
#---------------------------------------------------------------------------#

print time.strftime(timeFormat ,time.localtime())
print "7. We normalize the laplacian....",

def compute_laplacian():
	mylaplacian = np.zeros((NumberOfWordsForAnalysis,NumberOfWordsForAnalysis), dtype=np.float32 )

	for (i, j) in itertools.product(range(NumberOfWordsForAnalysis), repeat=2):
		if i == j:
			mylaplacian[i,j] = 1
		else:
			if incidencegraph[i,j] == 0:
				mylaplacian[i,j]=0
			else:
				mylaplacian[i,j] = -1 * incidencegraph[i,j]/ math.sqrt ( Diameter[i] * Diameter[j] )
	return mylaplacian

mylaplacian = compute_laplacian()
#mylaplacian = get_from_multiprocessing(compute_laplacian)

del incidencegraph
del Diameter

print "Done." 



#---------------------------------------------------------------------------#
#	Compute eigenvectors.
#---------------------------------------------------------------------------#

print time.strftime(timeFormat ,time.localtime())
print "8. Compute eigenvectors...",

myeigenvalues, myeigenvectors = np.linalg.eigh(mylaplacian)
#myeigenvalues, myeigenvectors = get_from_multiprocessing(np.linalg.eigh, [mylaplacian])

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

print time.strftime(timeFormat ,time.localtime())
print "10. Finding coordinates in space of low dimensionality."

def compute_coordinates():
    coordinates = dict()
    for wordno in range(NumberOfWordsForAnalysis):
        coordinates[wordno]= list() 
        for eigenno in range(NumberOfEigenvectors):
            coordinates[wordno].append( myeigenvectors[ wordno, eigenno ] )
    return coordinates

#coordinates = get_from_multiprocessing(compute_coordinates)
coordinates = compute_coordinates()

print time.strftime(timeFormat ,time.localtime())
print '       coordinates computed. Now computing distances between words...',

del myeigenvectors

datatype = ctypes.c_float

def compute_words_distance_1core():
    arr = np.zeros((NumberOfWordsForAnalysis, NumberOfWordsForAnalysis))

    for wordno1 in range(NumberOfWordsForAnalysis):
        for wordno2 in range(wordno1+1, NumberOfWordsForAnalysis):
            distance = 0
            for coordno in range(NumberOfEigenvectors):
                x = coordinates[wordno1][coordno] - coordinates[wordno2][coordno]
                distance += abs(x ** 3)
            arr[wordno1, wordno2] = distance
            arr[wordno2, wordno1] = distance
    return arr

wordsdistance = compute_words_distance_1core()

#wordsdistance_shared = mp.Array(datatype, NumberOfWordsForAnalysis ** 2)

def compute_words_distance(wordno1):
	arr = np.frombuffer(wordsdistance_shared.get_obj(), dtype=datatype)
	arr.shape = (NumberOfWordsForAnalysis, NumberOfWordsForAnalysis)

	for wordno2 in range(wordno1+1, NumberOfWordsForAnalysis):
		distance = 0
		for coordno in range(NumberOfEigenvectors):
			x = coordinates[wordno1][coordno] - coordinates[wordno2][coordno]
			distance += abs(x ** 3)
		arr[wordno1, wordno2] = distance
#		arr[wordno2, wordno1] = distance

	arr.shape = NumberOfWordsForAnalysis ** 2

#with closing(mp.Pool(processes=howManyCoresToUse, initializer=init, initargs=(wordsdistance_shared, ))) as p:
#	p.map_async(compute_words_distance, range(NumberOfWordsForAnalysis), chunksize=compute_chunksize(NumberOfWordsForAnalysis, howManyCoresToUse))

#p.join()

#wordsdistance = np.frombuffer(wordsdistance_shared.get_obj(), dtype=datatype)

print wordsdistance.shape

#wordsdistance.shape = (NumberOfWordsForAnalysis, NumberOfWordsForAnalysis)

#wordsdistance = wordsdistance + wordsdistance.T

print 'Done.'

#print time.strftime(timeFormat ,time.localtime())
#print '       Writing the word distance dict to CSV...',

#outCSVfile = open(outCSVname, 'w')

#for (wordno1, wordno2) in itertools.combinations(range(NumberOfWordsForAnalysis), 2):
#	outCSVfile.write('%d,%d,%f\n' % (wordno1, wordno2, wordsdistance[wordno1, wordno2]))

#outCSVfile.close()

#print 'Done.'

#---------------------------------------------------------------------------#
#	 Finding closest neighbors on the manifold's approximation
#---------------------------------------------------------------------------#

#print time.strftime(timeFormat ,time.localtime())
#print "11. Finding closest neighbors on the manifold('s approximation)."


#def compute_graph():
#	mygraph = nx.Graph()
#	for (wordno, word1)  in enumerate(analyzedwordlist):
#		mygraph.add_node(word1, wordindex=wordno)
#	return mygraph

#mygraph = get_from_multiprocessing(compute_graph)


print time.strftime(timeFormat ,time.localtime())
#print '      graph initialized. Computing nearest neighbors now... ',
print '      computing nearest neighbors now... ',



def compute_closest_neighbors():
	closestNeighbors = dict()

#	wordsdistanceSorted = sorted(wordsdistance.items(), key=lambda x:x[1])

	for (wordno1, word1) in enumerate(analyzedwordlist):

#		neighborWordNumberList = [wordno2 for (wordno2, distance) in sorted(wordsdistance[wordno1].items(), key=lambda x:x[1])]

		neighborWordNumberList = [wordno2 for (wordno2, distance) in sorted(enumerate(list(wordsdistance[wordno1])), key=lambda x:x[1])][1:]
		# ignore the first element in list because the first element is wordno1 itself

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

		neighborWordNumberList = neighborWordNumberList[: NumberOfNeighbors]

		closestNeighbors[wordno1] = neighborWordNumberList

#		for (idx, wordno2) in enumerate(neighborWordNumberList):
#			mygraph.add_edge(word1, analyzedwordlist[wordno2], rank=idx+1)


#		for (idx, wordno2) in enumerate(neighborWordNumberList):		
#			if wordno2 ==  wordno1:			 
#				continue			
#			if idx > NumberOfNeighbors:
#				break
#			word2 = analyzedwordlist[wordno2]

#			closestNeighbors[wordno1].append(wordno2)
#			mygraph.add_edge(word1, word2, rank=idx)		

	return closestNeighbors

closestNeighbors = compute_closest_neighbors()
#closestNeighbors = get_from_multiprocessing(compute_closest_neighbors)

for (wordno, word) in enumerate(analyzedwordlist):
	print >>outfileNeighbors, word, ' '.join([analyzedwordlist[idx] for idx in closestNeighbors[wordno]])

#for (word, neighborList) in sorted(closestNeighbors.items(), key=lambda x: analyzedworddict[x[0]]):
#	print >>outfileNeighbors, word, ' '.join(neighborList)

outfileNeighbors.close()

print 'done.'
print time.strftime(timeFormat ,time.localtime())

#nx.write_graphml(mygraph,outfilenameNeighborsGraphml)
 


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

endTime = time.localtime()

timeDifference = (time.mktime(endTime) - time.mktime(beginTime)) / 60

print time.strftime('log_%Y-%m-%d_%H-%M-%S' ,endTime)
print 'amount of time taken:', timeDifference, 'minutes'

#subprocess.call(('cp', outfilenameFromWordToContexts, '.'))
subprocess.call(('cp', outfilenameNeighbors, '.'))
