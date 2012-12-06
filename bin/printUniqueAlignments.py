import importSomaData
import MatchResult
pickleFile = 'fusarium_SigMatches.refinedUniqueMatchList.pickle'
contigToMatch = importSomaData.getContigToUniqueMatch(pickleFile)
ml = sorted(contigToMatch.values(), key = lambda mr: mr.size)[::-1]

for mr in ml:
    print '*'*50
    MatchResult.printAlignment(mr)
