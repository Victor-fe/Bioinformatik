seq = 'CTTCTCACGTACAACAAAATC'

symbol2number = {"A":0,"C":1,"G":2,"T":3}
    

def PatternToNumber(Pattern):
    if not Pattern:
        return 0
    symbol = Pattern[-1]
    prefix = Pattern[:-1]
    return ((4*PatternToNumber(prefix))+symbol2number[symbol])

print(PatternToNumber(seq))  