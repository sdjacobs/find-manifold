from findManifold import *
import argparse

def makeArgParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("nWords", help="Number of words for analysis", type=int, default=9)
    parser.add_argument("nNeighbors", help="Number of neighbors", type=int, default=1000)
    parser.add_argument("nEigenvectors", help="Number of eigenvectors", type=int, default=11)
    parser.add_argument("--bigrams", help="Bigrams file to use", type=str,
            default="../../data/english/ngrams/english-brown_bigrams.txt")
    parser.add_argument("--trigrams", help="Trigrams file to use", type=str,
            default="../../data/english/ngrams/english-brown_trigrams.txt")
    parser.add_argument("--words", help="Words file to use", type=str,
            default="../../data/english/ngrams/english-brown_words.txt")
    parser.add_argument("--output", help="Output folder to use", type=str,
            default="output")
    parser.add_argument("--name", help="Corpus name", type=str, default="english-brown")
    parser.add_argument("--languagename", help="Language name", type=str, default="english")


    return parser


def main(argv):

    punctuation         = " $/+.,;:?!()\"[]"

    args = makeArgParser().parse_args()
    NumberOfWordsForAnalysis = args.nWords
    NumberOfNeighbors = args.nNeighbors
    NumberOfEigenvectors = args.nEigenvectors
    infileBigramsname = args.bigrams
    infileTrigramsname = args.trigrams
    infileWordsname = args.words
    outfilenameNeighbors = args.output + "/" + args.name + "_" + str(NumberOfWordsForAnalysis) + "_" + str(NumberOfNeighbors) + "_nearest_neighbors.txt"
    outfilenameContexts     = args.output + "/" + args.name + "_contexts.txt"
    outfilenameFromWordToContexts     = args.output + "/" + args.name + "_" + str(NumberOfWordsForAnalysis)  + "_from-word-to-contexts.txt"

    print "\nI am looking for: ", infileTrigramsname
    print "Number of words that will be analyzed:", NumberOfWordsForAnalysis
    print "Number of neighbors:", NumberOfNeighbors

    # initialize data structures
    closestNeighbors     = collections.OrderedDict() #a dict whose values are lists; the lists are the closest words to the key.
    contexts         = dict() # key is word, value is a dict of contexts (used to be a list of contexts)
    coordinates         = dict()
    Diameter = dict()

    from_word_to_context = collections.defaultdict(collections.Counter) # this dict takes a word as key, and returns a collections.Counter dict as the value; the value is a dict with (context, frequency count) pairs.

    wordsdistance         = dict() # key is a word, word1,  being analyzed, value is a pair of word-index-number and euclidean distance (word2, distance). This will be sorted to get the nearest neighbors to word1.


    # open files
    outfileFromWordToContexts = open(outfilenameFromWordToContexts, 'w')
    outfileNeighbors = open(outfilenameNeighbors, 'w')
    outfileContexts = open(outfilenameContexts, 'w')
    wordfile        = open(infileWordsname)
    trigramfile         = open(infileTrigramsname)
    bigramfile         = open(infileBigramsname)
    


    for outfile in [outfileNeighbors, outfileFromWordToContexts]:
        print >>outfile, "# language:", \
                args.languagename, "\n# corpus:",\
                args.name, "\n#",\
                "Number of words analyzed", NumberOfWordsForAnalysis,"\n#", \
                "Number of neighbors identified", NumberOfNeighbors,"\n"

    print >>outfileContexts, "#  The number with each context is the number of distinct words found in that context.\n#" 

    mywords = GetMyWords(wordfile)
    wordfile.close()
    print "1. Word file is ", infileWordsname, '\t corpus has', len(mywords), 'words'


    if NumberOfWordsForAnalysis > len(mywords):
        NumberOfWordsForAnalysis = len(mywords)
        print 'number of words for analysis reduced to', NumberOfWordsForAnalysis

    analyzedwordlist = mywords.keys()[ : NumberOfWordsForAnalysis]
    analyzedwordset = set(analyzedwordlist)
    
    print "Reading in trigrams and bigrams"
    ReadInTrigrams(trigramfile, analyzedwordlist, analyzedwordset, from_word_to_context)
    ReadInBigrams(bigramfile, analyzedwordlist, analyzedwordset, from_word_to_context)
   
    print "Counting shared contexts"
    context_array = MakeContextArray(NumberOfWordsForAnalysis, from_word_to_context)
    CountOfSharedContexts = counting_context_features(context_array) # TODO: which one?

    print "Computing diameter"
    Diameter = Normalize(NumberOfWordsForAnalysis, CountOfSharedContexts)

    print "Computing incidence graph"
    incidencegraph = compute_incidence_graph(NumberOfWordsForAnalysis, Diameter, CountOfSharedContexts)
    
    print "Computing mylaplacian"
    mylaplacian = compute_laplacian(NumberOfWordsForAnalysis, Diameter, incidencegraph)

    print "8. compute eigenvectors...",
    myeigenvalues, myeigenvectors = np.linalg.eigh(mylaplacian)
    
    print "10. finding coordinates in space of low dimensionality."
    Coordinates =  compute_coordinates(NumberOfWordsForAnalysis, NumberOfEigenvectors, myeigenvectors)
    
    print '       coordinates computed. now computing distances between words...',
    wordsdistance = compute_words_distance(NumberOfWordsForAnalysis, NumberOfEigenvectors, Coordinates)

    print '      computing nearest neighbors now... ',
    closestNeighbors = compute_closest_neighbors(analyzedwordlist, wordsdistance, NumberOfNeighbors)

    print "Output to files"
    for (wordno, word) in enumerate(analyzedwordlist):
        print >>outfileNeighbors, word, ' '.join([analyzedwordlist[idx] for idx in closestNeighbors[wordno]])

    outfileNeighbors.close()



# Don't execute any code if we are loading this file as a module.
if __name__ == "__main__":
    main(sys.argv)
