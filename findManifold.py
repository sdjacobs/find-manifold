#-----------------------------------------------------------------------#
#
#    This program takes n-gram files and a word list    
#    and creates a file with lists of most similar words.
#    John Goldsmith and Wang Xiuli 2012.
#    Jackson Lee 2014
#
#-----------------------------------------------------------------------#

import codecs
#import os
import sys
#import string
import operator
import math
import collections
import numpy as np
#import networkx as nx
import ctypes
import itertools
import time
#import subprocess
#from contextlib import closing

import multiprocessing as mp

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

def GetMyWords(wordfile):
    mywords = collections.OrderedDict()
    for line in wordfile:
        pieces = line.split()
        if pieces[0] == "#":
            continue
        mywords[pieces[0]] = int(pieces[1])         
    return mywords

def ReadInTrigrams(trigramfile, analyzedwordlist, analyzedwordset, from_word_to_context):
    for line in trigramfile:
        if line.startswith('#'):
            continue
        thesewords = line.split()

        thisword = thesewords[1]
        if thisword in analyzedwordset:
            context = thesewords[0] + " __ " +  thesewords[2]
            wordno = analyzedwordlist.index(thisword)
            from_word_to_context[wordno][context] += 1

        #Left trigrams
        thisword = thesewords[0]
        if thisword in analyzedwordset:
            context = " __ " + thesewords[1] + " " + thesewords[2]
            wordno = analyzedwordlist.index(thisword)
            from_word_to_context[wordno][context] += 1

        #Right trigrams
        thisword = thesewords[2]
        if thisword   in analyzedwordset:    
            context = thesewords[0] + " " + thesewords[1] + " __ "
            wordno = analyzedwordlist.index(thisword)
            from_word_to_context[wordno][context] += 1

def ReadInBigrams(bigramfile, analyzedwordlist, analyzedwordset, from_word_to_context):
    print "...Reading in bigram file."
    for line in bigramfile:
        thesewords = line.split()
        if thesewords[0] == "#":
            continue 
        thisword = thesewords[1]
        if thisword in analyzedwordset:
            context = thesewords[0] + " __ " 
            wordno = analyzedwordlist.index(thisword) # TODO: is this a bug in the older version?
            from_word_to_context[wordno][context] += 1
        thisword = thesewords[0]
        if thisword in analyzedwordset:
            context = "__ " + thesewords[1]
            wordno = analyzedwordlist.index(thisword)
            from_word_to_context[wordno][context] += 1

def MakeContextArray(NumberOfWordsForAnalysis, from_word_to_context):
    context_list = list(set(context for i in from_word_to_context for context in from_word_to_context[i]))
    context_array = np.zeros((NumberOfWordsForAnalysis, len(context_list)))
    for wordno in range(NumberOfWordsForAnalysis):
        for contextno in range(len(context_list)):
            if (from_word_to_context[wordno][(context_list[contextno])] > 0):
                context_array[wordno, contextno] = 1
    return context_array

def QuickGetNumberOfSharedContexts(word1, word2):
    return np.dot(context_array[word1], context_array[word2])

def GetNumberOfSharedContexts(word1, word2, from_word_to_context):
    return len(set(from_word_to_context[word1]) & set(from_word_to_context[word2]))

def counting_context_features_old(NumberOfWordsForAnlysis, from_word_to_context):
    arr = np.zeros((NumberOfWordsForAnalysis, NumberOfWordsForAnalysis))
    for word1 in range(0, NumberOfWordsForAnalysis):
        for word2 in range(word1+1, NumberOfWordsForAnalysis):
            x = GetNumberOfSharedContexts(word1, word2, from_word_to_context)
            arr[word1, word2] = x
            arr[word2, word1] = x
    return arr

def counting_context_features(context_array):
    return np.dot(context_array, context_array.T) 


def compute_incidence_graph(NumberOfWordsForAnalysis, Diameter, CountOfSharedContexts):
    incidencegraph= np.zeros( (NumberOfWordsForAnalysis, NumberOfWordsForAnalysis), dtype=np.int32)

    for (w1, w2) in itertools.product(range(NumberOfWordsForAnalysis), repeat=2):
        if w1 == w2:
            incidencegraph[w1,w1] = Diameter[w1]
        else:
            incidencegraph[w1,w2] = CountOfSharedContexts[w1,w2]    
    return incidencegraph


def compute_laplacian(NumberOfWordsForAnalysis, Diameter, incidencegraph):
    mylaplacian = np.zeros((NumberOfWordsForAnalysis, NumberOfWordsForAnalysis), dtype=np.float32 )

    for (i, j) in itertools.product(range(NumberOfWordsForAnalysis), repeat=2):
        if i == j:
            mylaplacian[i,j] = 1
        else:
            if incidencegraph[i,j] == 0:
                mylaplacian[i,j]=0
            else:
                mylaplacian[i,j] = -1 * incidencegraph[i,j]/ math.sqrt ( Diameter[i] * Diameter[j] )
    return mylaplacian

def LatexAndEigenvectorOutput(LatexFlag, PrintEigenvectorsFlag, infileWordsname, outfileLatex, outfileEigenvectors, NumberOfEigenvectors, myeigenvalues, NumberOfWordsForAnalysis):
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


def compute_coordinates(NumberOfWordsForAnalysis, NumberOfEigenvectors, myeigenvectors):
    Coordinates = dict()
    for wordno in range(NumberOfWordsForAnalysis):
        Coordinates[wordno]= list() 
        for eigenno in range(NumberOfEigenvectors):
            Coordinates[wordno].append( myeigenvectors[ wordno, eigenno ] )
    return Coordinates

def compute_words_distance(nwords, neigs, coordinates):
    arr = np.zeros((nwords, nwords))

    for wordno1 in range(nwords):
        for wordno2 in range(wordno1+1, nwords):
            distance = 0
            for coordno in range(neigs):
                x = coordinates[wordno1][coordno] - coordinates[wordno2][coordno]
                distance += abs(x ** 3)
            arr[wordno1, wordno2] = distance
            arr[wordno2, wordno1] = distance
    return arr

def compute_closest_neighbors(analyzedwordlist, wordsdistance, NumberOfNeighbors):
    closestNeighbors = dict()
    for (wordno1, word1) in enumerate(analyzedwordlist):
        neighborWordNumberList = [wordno2 for (wordno2, distance) in sorted(enumerate(list(wordsdistance[wordno1])), key=lambda x:x[1])][1:]
        neighborWordNumberList = neighborWordNumberList[: NumberOfNeighbors]
        closestNeighbors[wordno1] = neighborWordNumberList
    return closestNeighbors

def main(argv):
    timeFormat = '      current time: %Y-%m-%d %H:%M:%S'

    beginTime = time.localtime()
    print
    print time.strftime(timeFormat , beginTime)
   
#---------------------------------------------------------------------------#
#    Variables to be changed by user
#---------------------------------------------------------------------------#
    LatexFlag = False
    PrintEigenvectorsFlag = False
    unicodeFlag = False 
    FileEncoding =  "ascii"

    shortfilename         = "english-brown"
    outshortfilename     = "english-brown"
    languagename         = "english"

    datafolder            = "../../data/"

    ngramfolder           = datafolder + languagename + "/ngrams/"
    outfolder             = datafolder + languagename + "/neighbors/"
    wordcontextfolder    = datafolder + languagename + "/word_contexts/"

    NumberOfEigenvectors         = 11

    punctuation         = " $/+.,;:?!()\"[]"

    try:
        if len(argv) < 2:
            NumberOfNeighbors = 9
            NumberOfWordsForAnalysis = 1000
            continueKeyboard = raw_input('default values: NumberOfWordsForAnalysis = 1000, NumberOfNeighbors = 9\nContinue? [N/y] ')
            if continueKeyboard.lower() != 'y':
                print '\nprogram terminated\n'
                exit()
        elif len(argv) == 2:
            NumberOfNeighbors = 9
            NumberOfWordsForAnalysis = int(argv[1])
        else:
            NumberOfNeighbors = int(argv[2])
            NumberOfWordsForAnalysis = int(argv[1])
    except: 
        print 'usage: python findManifold.py [NumberOfWordsForAnalysis] [NumberOfNeighbors]'
        exit()

#---------------------------------------------------------------------------#
#    File names
#---------------------------------------------------------------------------#

    infileBigramsname     = ngramfolder + shortfilename    + "_bigrams.txt" 
    infileTrigramsname     = ngramfolder + shortfilename    + "_trigrams.txt"
    infileWordsname     = ngramfolder   + shortfilename    + "_words.txt" 
    outfilenameEigenvectors    = outfolder     + outshortfilename + "_words_eigenvectors" + ".txt"
    outfilenameNeighbors    = outfolder     + outshortfilename + "_" + str(NumberOfWordsForAnalysis) + "_" + str(NumberOfNeighbors) + "_nearest_neighbors.txt"
    outfilenameLatex     = outfolder     + outshortfilename + "_latex.tex"
    outfilenameContexts     = outfolder     + outshortfilename + "_contexts.txt"
    outfilenameFromWordToContexts     = outfolder     + outshortfilename + "_" + str(NumberOfWordsForAnalysis)  + "_from-word-to-contexts.txt"

    print "\nI am looking for: ", infileTrigramsname
    print "Number of words that will be analyzed:", NumberOfWordsForAnalysis
    print "Number of neighbors:", NumberOfNeighbors

#---------------------------------------------------------------------------#
#    Important data structures
#---------------------------------------------------------------------------#
    closestNeighbors     = collections.OrderedDict() #a dict whose values are lists; the lists are the closest words to the key.
    contexts         = dict() # key is word, value is a dict of contexts (used to be a list of contexts)
    coordinates         = dict()
    Diameter = dict()

    from_word_to_context = collections.defaultdict(collections.Counter) # this dict takes a word as key, and returns a collections.Counter dict as the value; the value is a dict with (context, frequency count) pairs.

    wordsdistance         = dict() # key is a word, word1,  being analyzed, value is a pair of word-index-number and euclidean distance (word2, distance). This will be sorted to get the nearest neighbors to word1.

#---------------------------------------------------------------------------#
#    Open files for reading and writing
#---------------------------------------------------------------------------#

    if unicodeFlag:
        trigramfile         =codecs.open(infileTrigramsname, encoding = FileEncoding)
        wordfile         =codecs.open(infileWordsname, encoding = FileEncoding)
        if PrintEigenvectorsFlag:
            outfileEigenvectors = codecs.open (outfilename1, "w",encoding = FileEncoding)
        outfileNeighbors    = codecs.open (outfileneighborsname, "w",encoding = FileEncoding)
    else:
        if PrintEigenvectorsFlag:
            outfileEigenvectors = open (outfilenameEigenvectors, "w")
        outfileNeighbors    = open (outfilenameNeighbors, "w")
        outfileLatex         = open (outfilenameLatex, "w")
        outfileContexts     = open (outfilenameContexts, "w")
        outfileFromWordToContexts = open(outfilenameFromWordToContexts, 'w')

        wordfile        = open(infileWordsname)
        trigramfile         = open(infileTrigramsname)
        bigramfile         = open(infileBigramsname)
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

    mywords = GetMyWords(wordfile)
    print "1. Word file is ", infileWordsname, '\t corpus has', len(mywords), 'words'
    wordfile.close()

    if NumberOfWordsForAnalysis > len(mywords):
        NumberOfWordsForAnalysis = len(mywords)
        print 'number of words for analysis reduced to', NumberOfWordsForAnalysis

    analyzedwordlist = mywords.keys()[ : NumberOfWordsForAnalysis]
    analyzedwordset = set(analyzedwordlist)
    print "2. Reading in trigram file."
    ReadInTrigrams(trigramfile, analyzedwordlist, analyzedwordset, from_word_to_context)
    print "... reading in bigrams" 
    ReadInBigrams(bigramfile, analyzedwordlist, analyzedwordset, from_word_to_context)
    print "3. End of words and counts."
   

#---------------------------------------------------------------------------#
#    Count context features shared by words
#---------------------------------------------------------------------------#

    print time.strftime(timeFormat ,time.localtime())
    print "4. Counting context features shared by words...",
    datatype = ctypes.c_int
    context_array = MakeContextArray(NumberOfWordsForAnalysis, from_word_to_context)
    CountOfSharedContexts = counting_context_features(context_array) # TODO: which one?
    print 'done.'

#---------------------------------------------------------------------------#
#    Normalize function call
#---------------------------------------------------------------------------#


    print time.strftime(timeFormat ,time.localtime())
    print "5. Normalizing nearness measurements....",
    Diameter = Normalize(NumberOfWordsForAnalysis, CountOfSharedContexts)
    print "\t Done." 


#---------------------------------------------------------------------------#
#    Incidence graph
#---------------------------------------------------------------------------#

    print time.strftime(timeFormat ,time.localtime())
    print "6. We compute the incidence graph....",
    incidencegraph = compute_incidence_graph(NumberOfWordsForAnalysis, Diameter, CountOfSharedContexts)
    print "done." 

         
#---------------------------------------------------------------------------#
#    normalize the laplacian.
#---------------------------------------------------------------------------#

    print time.strftime(timeFormat ,time.localtime())
    print "7. we normalize the laplacian....",
    mylaplacian = compute_laplacian(NumberOfWordsForAnalysis, Diameter, incidencegraph)
    print "done." 

#---------------------------------------------------------------------------#
#    compute eigenvectors.
#---------------------------------------------------------------------------#

    print time.strftime(timeFormat ,time.localtime())
    print "8. compute eigenvectors...",
    myeigenvalues, myeigenvectors = np.linalg.eigh(mylaplacian)
    print "done."

  
#---------------------------------------------------------------------------#
#    generate latex and eigenvector output.
#---------------------------------------------------------------------------#
    if LatexFlag and PrintEigenvectorsFlag:
        LatexAndEigenvectorOutput(LatexFlag, PrintEigenvectorsFlag, infileWordsname, outfileLatex, outfileEigenvectors, NumberOfEigenvectors, myeigenvalues, NumberOfWordsForAnalysis)

#---------------------------------------------------------------------------#
#    finding coordinates in space of low dimensionality
#---------------------------------------------------------------------------#

    print time.strftime(timeFormat ,time.localtime())
    print "10. finding coordinates in space of low dimensionality."
    Coordinates =  compute_coordinates(NumberOfWordsForAnalysis, NumberOfEigenvectors, myeigenvectors)
    print time.strftime(timeFormat ,time.localtime())
    
    print '       coordinates computed. now computing distances between words...',
    wordsdistance = compute_words_distance(NumberOfWordsForAnalysis, NumberOfEigenvectors, Coordinates)
    print 'Done.'

    print time.strftime(timeFormat ,time.localtime())
    print '      computing nearest neighbors now... ',

    closestNeighbors = compute_closest_neighbors(analyzedwordlist, wordsdistance, NumberOfNeighbors)


#---------------------------------------------------------------------------#
# Output to files 
#---------------------------------------------------------------------------#
    for (wordno, word) in enumerate(analyzedwordlist):
        print >>outfileNeighbors, word, ' '.join([analyzedwordlist[idx] for idx in closestNeighbors[wordno]])

    outfileNeighbors.close()

    print 'done.'
    print time.strftime(timeFormat, time.localtime())

    if PrintEigenvectorsFlag:
        outfileEigenvectors.close()

    endTime = time.localtime()
    timeDifference = (time.mktime(endTime) - time.mktime(beginTime)) / 60
    print time.strftime('log_%Y-%m-%d_%H-%M-%S' ,endTime)
    print 'amount of time taken:', timeDifference, 'minutes'

    # I would recommend that you do file management stuff in the sbatch script, and leave this script for the numerical processing.
    # subprocess.call(('cp', outfilenameNeighbors, '.'))

# Don't execute any code if we are loading this file as a module.
if __name__ == "__main__":
    main(sys.argv)
