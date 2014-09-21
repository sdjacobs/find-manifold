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
    print "Word file is ", infileWordsname, '\t corpus has', len(mywords), 'words'


    if NumberOfWordsForAnalysis > len(mywords):
        NumberOfWordsForAnalysis = len(mywords)
        print 'number of words for analysis reduced to', NumberOfWordsForAnalysis

    analyzedwordlist = mywords.keys()[ : NumberOfWordsForAnalysis] 
   
    print "Reading bigrams/trigrams"

    context_array = GetContextArray(NumberOfWordsForAnalysis, mywords, bigramfile, trigramfile)
    
    print "Computing shared contexts"
    CountOfSharedContexts = context_array.dot(context_array.T).todense()

    print "Computing diameter"
    Diameter = Normalize(NumberOfWordsForAnalysis, CountOfSharedContexts)

    print "Computing incidence graph"
    incidencegraph = compute_incidence_graph(NumberOfWordsForAnalysis, Diameter, CountOfSharedContexts)
    
    print "Computing mylaplacian"
    mylaplacian = compute_laplacian(NumberOfWordsForAnalysis, Diameter, incidencegraph)

    print "Compute eigenvectors..."
    myeigenvalues, myeigenvectors = np.linalg.eigh(mylaplacian)
    
    print 'Coordinates computed. now computing distances between words...'
    coordinates = myeigenvectors[:,:NumberOfEigenvectors] # take first N columns of eigenvector matrix
    wordsdistance = compute_words_distance(NumberOfWordsForAnalysis, coordinates)

    print 'Computing nearest neighbors now... '
    closestNeighbors = compute_closest_neighbors(analyzedwordlist, wordsdistance, NumberOfNeighbors)

    print "Output to files"
    for (wordno, word) in enumerate(analyzedwordlist):
        print >>outfileNeighbors, word, ' '.join([analyzedwordlist[idx] for idx in closestNeighbors[wordno]])

    outfileNeighbors.close()



# Don't execute any code if we are loading this file as a module.
if __name__ == "__main__":
    main(sys.argv)
