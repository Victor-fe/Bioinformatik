seq = 'CTTCTCACGTACAACAAAATC'

symbol2number = {"A":0,"C":1,"G":2,"T":3}
    

def PatternToNumber(Pattern):
    if not Pattern:
        return 0
    symbol = Pattern[-1]
    prefix = Pattern[:-1]
    return ((4*PatternToNumber(prefix))+symbol2number[symbol])

def NumberToPattern(index, k):
	bases = ['A', 'C', 'G', 'T']
	pattern = ''
	for i in range(k):
		pattern += bases[index % 4]
		index = index // 4
	return pattern[::-1]

def ComputingFrequencies(text,k):
    FrequencyArray =[]
    for i in range(0,((4**k))):
        FrequencyArray.append(0)
    for i in range(0,(len(text)-1)):
        pattern = text[i:(i+k)]
        j = PatternToNumber(pattern)
        FrequencyArray[j] = FrequencyArray[j]+1
    return FrequencyArray

def FasterFrequentWords(text,k):
    FrequentPatterns = []
    FrequencyArray = ComputingFrequencies(text,k)
    maxCount = max(FrequencyArray)
    for i in range(0,(4**k)):
        if FrequencyArray[i] == maxCount:
            pattern = NumberToPattern(i,k)
            FrequentPatterns.append(pattern)
    return FrequentPatterns

print(FasterFrequentWords("ACGCGGCTCTGAAA",2))