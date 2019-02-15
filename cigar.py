import re


def cigar_parse(cigar,start):
    insertions = []
    alignments = []
    aligns=[]

    pattern = re.compile('([MIDNSHPX=])')
    values = pattern.split(cigar)[:-1] ## turn cigar into tuple of values
    paired = (values[n:n+2] for n in range(0,len(values),2)) ## pair values by twos
    i = 0 ## alignment coordinate index
    g = start ## genomic coordinate index
    gstop=0
    for pair in paired:
        l = int(pair[0]) ## length of CIGAR event
        t = pair[1] ## type of CIGAR event
        if t == 'M': ## if match, return consecutive coordinates
            alignments.append((g, g+i+l-i,(i, i + l))) ## (genomic offset, (alignment.start, alignment.end))
            aligns.append((g, g+i+l-i, i, i+l)) ## (genomic offset, (alignment.start, alignment.end))
            i += l
            g += l
            gstop=g+i+l-i
        elif t == 'D': ## skip 'l' number of coordinates in reference
            g += l
        elif t == 'I': ## insertion of 'l' length
            insertions.append((i, i + l))
            i += l
        elif t == 'N': ## skipped region from the reference
            g += l
        elif t == 'S': ## soft clipping (clipped sequences present in SEQ)
            i += l
        elif t == 'H': ## hard clipping (clipped sequences NOT present in SEQ)
            pass
        elif t == 'P': ## padding (silent deletion from padded reference)
            pass
        elif t == '=': ## sequence match
            pass
        elif t == 'X': ## sequence mismatch
            pass
    return gstop
