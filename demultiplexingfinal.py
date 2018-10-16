#Index hops here require having no N's and an average quality score for each index of > 25, and are othewise considered bad reads
import argparse
import numpy as np
import gzip
def filelist():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--read1', help = "Assigns first read file",
        required = True, type=str)
    parser.add_argument('-R', '--read2', help = "Assigns second read file",
        required = True, type=str)
    parser.add_argument('-i', '--index1', help = "Assigns first index file",
        required = True, type=str)
    parser.add_argument('-I', '--index2', help = "Assigns second index file",
        required = True, type=str)
    return parser.parse_args()
args = filelist()
read1 = args.read1
read2 = args.read2
index1 = args.index1
index2 = args.index2
linecount = 0
mismatch = 0
hopped = 0
indexlist = {
    "GTAGCGTA":"B1",
    "CGATCGAT":"A5",
    "GATCAAGG":"C1",
    "AACAGCGA":"B9",
    "TAGCCATG":"C9",
    "CGGTAATC":"C3",
    "CTCTGGAT":"B3",
    "TACCGGAT":"C4",
    "CTAGCTCA":"A11",
    "CACTTCAC":"C7",
    "GCTACTCT":"B2",
    "ACGATCAG":"A1",
    "TATGGCAC":"B7",
    "TGTTCCGT":"A3",
    "GTCCTAAG":"B4",
    "TCGACAAG":"A12",
    "TCTTCGAC":"C10",
    "ATCATGCG":"A2",
    "ATCGTGGT":"C2",
    "TCGAGAGT":"A10",
    "TCGGATTC":"B8",
    "GATCTTGC":"A7",
    "AGAGTCCA":"B10",
    "AGGATAGC":"A8"
} #Dictionary of index IDs against actual barcodes
indexmismatch = {
    "B1":0,
    "A5":0,
    "C1":0,
    "B9":0,
    "C9":0,
    "C3":0,
    "B3":0,
    "C4":0,
    "A11":0,
    "C7":0,
    "B2":0,
    "A1":0,
    "B7":0,
    "A3":0,
    "B4":0,
    "A12":0,
    "C10":0,
    "A2":0,
    "C2":0,
    "A10":0,
    "B8":0,
    "A7":0,
    "B10":0,
    "A8":0
} #Dictionary for total number of incorrect reads per index
indexreads = {
    "B1":0,
    "A5":0,
    "C1":0,
    "B9":0,
    "C9":0,
    "C3":0,
    "B3":0,
    "C4":0,
    "A11":0,
    "C7":0,
    "B2":0,
    "A1":0,
    "B7":0,
    "A3":0,
    "B4":0,
    "A12":0,
    "C10":0,
    "A2":0,
    "C2":0,
    "A10":0,
    "B8":0,
    "A7":0,
    "B10":0,
    "A8":0
} #Dictionary for total number of reads for each index

def filemaker(indexdictionary):
    """Generates output files from index dictionary"""
    filelist = {}
    for key in indexlist.keys():
        filelist[indexlist[key]] = [str(indexlist[key] + ".1out.fastq"), str(indexlist[key] + ".2out.fastq")]
    filelist["badindex"] = ["badindex.1out.fastq", "badindex.2out.fastq"]
    return(filelist)

filelist = filemaker(indexlist)


def convert_phred(letter):
    """Converts a single character into a phred score"""
    letter = ord(letter)-33
    return letter

def rev_comp(barcode1, barcode2):
    """Tests two DNA strings to see if they are the reverse-complements of one-another"""
    for x in range(0,len(barcode1)):
        verify = 0
        if barcode1[x] == "A":
            if barcode2[x] == "T":
                continue
        elif barcode1[x] == "C":
            if barcode2[x] == "G":
                continue
        elif barcode1[x] == "G":
            if barcode2[x] == "C":
                continue
        elif barcode1[x] == "T":
            if barcode2[x] == "A":
                continue
        else:
            verify = 1
    if verify == 0:
        return(True)
    else:
        return(False)

def n_test(barcode1, barcode2):
    """Finds Ns in barcodes or not"""
    if "N" in str(barcode1) or "N" in str(barcode2):
        return(False)
    else:
        return(True)

def index_quality(barcode1, barcode2):
    """Finds if indices are high-enough quality to verify, using 25 as a cutoff"""
    b1score = 0
    b2score = 0
    for x in barcode1:
        val1 = convert_phred(x)
        b1score += val1
    for y in barcode2:
        val2 = convert_phred(y)
        b2score += val2
    b1avg = b1score/len(barcode1)
    b2avg = b2score/len(barcode2)
    if b1avg < 25 or b2avg < 25:
        return(False)
    else:
        return(True)

def id_index(barcode1, barcode2):
    """Figures out which files to output to"""

    if barcode1 in indexlist.keys():
        output_index = filelist[indexlist[barcode1]]
        actual = barcode1
    elif barcode2 in indexlist.keys():
        output_index = filelist[indexlist[barcode2]]
        actual = barcode2
    else:
        output_index = filelist["badindex"]
        actual = "badindex"
    return(output_index, actual)

def hop_test(barcode1, barcode2):
    """Tests if mismatch is due to index hopping or a bad read"""
    if barcode1 in indexlist.keys() and barcode2 in indexlist.keys():
        return(True)
    else:
        return(False)

def writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence_, r2_sequence, r2_space, r2_score, barcode):
    """Writes out the header, sequence, spacer, and quality scores to their specified files"""
    with open(filelist[barcode][0], "a") as o1, open(filelist[barcode][1], "a") as o2:
        o1.write(str(str(r1_header) + ":" + str(i1_sequence) + "\n" + str(r1_sequence) + "\n" + str(r1_space) + "\n" + str(r1_score) + "\n"))
        o2.write(str(str(r2_header) + ":" + str(i2_sequence) + "\n" + str(r2_sequence) + "\n" + str(r2_space) + "\n" + str(r2_score) + "\n"))
        o1.close()
        o2.close()
    return(None)

#for x in range(0,2): # Loop finds quality distribution for the input files
#    with gzip.open(fil[x], "rt") as f:
#        with open(distfiles[x], "w") as out:
#            LN = 0
#            for line in f:
#                if LN == 1:
#                    qscores = np.zeros(len(line)-1, dtype=int)
#                if LN % 4 == 3:
#                    for q in range(len(line)-1):
#                        qscores[q] += convert_phred(str(line[q]))
#                LN += 1
#            qscores = qscores/(LN/4)
#            np.savetxt(out, qscores)
#    f.close()

with gzip.open(read1, 'rt') as r1, gzip.open(read2, 'rt') as r2, gzip.open(index1, 'rt') as i1, gzip.open(index2, 'rt') as i2: # Beginning of demultiplexing loop
    while True:

        #read1 info
        r1_header = r1.readline().strip()
        r1_sequence = r1.readline().strip()
        r1_space = r1.readline().strip()
        r1_score = r1.readline().strip()

        #read2 info
        r2_header = r2.readline().strip()
        r2_sequence = r2.readline().strip()
        r2_space = r2.readline().strip()
        r2_score = r2.readline().strip()


        #index1 info
        i1_header = i1.readline().strip()
        i1_sequence = i1.readline().strip()
        i1_space = i1.readline().strip()
        i1_score = i1.readline().strip()

        #index2 info
        i2_header = i2.readline().strip()
        i2_sequence = i2.readline().strip()
        i2_space = i2.readline().strip()
        i2_score = i2.readline().strip()

        if len(r1_header) == 0 and len(r1_sequence) == 0 and len(r1_score) == 0:
            break
        linecount += 1

        if n_test(i1_sequence, i2_sequence) == True:
            if index_quality(i1_score, i2_score):
                if hop_test(i1_sequence, i2_sequence) == True:
                    hopped += 1
                    print("hopped")
                    print(linecount)
                    writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence, r2_sequence, r2_space, r2_score, "badindex")
                else:
                    if rev_comp(i1_sequence, i2_sequence) == True:
                        barcode, real = id_index(i1_sequence, i2_sequence)
                        if real == "badindex":
                            writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence, r2_sequence, r2_space, r2_score, "badindex")
                        else:
                            id = indexlist[real]
                            indexreads[id] +=  1
                            print(id)
                            writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence, r2_sequence, r2_space, r2_score, id)
                    else:
                        mismatch += 1
                        writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence, r2_sequence, r2_space, r2_score, "badindex")
            else:
                mismatch += 1
                writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence, r2_sequence, r2_space, r2_score, "badindex")
        else:
            mismatch += 1
            writeout(r1_header, i1_sequence, r1_sequence, r1_space, r1_score, r2_header, i2_sequence, r2_sequence, r2_space, r2_score, "badindex")
for key in indexreads:
    print("Index ", key, " had ", indexreads[key], " reads. This was ", round((indexreads[key]/linecount)*100, 2), "% of the sample.", sep = "")
print("There were", linecount, "total reads.")
print("There were", mismatch, "total bad reads not due to index hopping.")
print("There were ", hopped, " reads that had index hopping, which is ", round((hopped/linecount)*100,2), "% of total reads.", sep = "")
