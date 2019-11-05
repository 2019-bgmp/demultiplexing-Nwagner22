import gzip
import numpy as np
import matplotlib.pyplot as plt
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="A program to k-merize our data set at different sizes for k and plot the result")
    parser.add_argument("-f", "--file", help="use to specify file name", required=True, type = str)
    parser.add_argument("-l", "--length", help="use to specify length of line (101 for sequence and 8 for index)", required=True, type = int)
    parser.add_argument("-o", "--output", help="use to specify output file name", type = str)
    return parser.parse_args()

args = get_args()               # calls get_args method from above assigns the arguments to args
GIVEN_FILE = args.file          # assigning file name as string to global varible
LINE_LENGTH = args.length       # assigning line length as int to global variable
OUTPUT_FILE = args.output       # assigning output file name as string to global variable

def convert_phred(letter):                                # method to convert a single character into a phred score
    """Converts a single character into a phred score""" 
    return(ord(letter)-33)                                # ord returns the unicode code for any given unicode character. Then 33 is subracted off because that is the standard in order to avoid the first 33 non-visable characters


R1 = np.zeros((LINE_LENGTH),dtype=float)
LN = 0
len_line1 = 0
with gzip.open(GIVEN_FILE, 'rt') as fastqFile:
    for line in fastqFile:
        if(LN%4==3):
            line = line.strip()
            for num in range(len(line)):
                len_line1 = len(line)
                R1[num] += convert_phred(line[num])    # going through every element in the line and keeping a running sum of each phred value at each position
        LN += 1
for i in range(len(R1)):
    R1[i] = R1[i] / (LN / 4)     # calculating the mean phred score at each position

plt.figure(0)
plt.plot(range(len_line1),R1,'ro')   # bar plot to compare the base pair number and the mean quality score
plt.xlabel("#_Base_Pair")
plt.ylabel("Mean_Quality_Score")
plt.title("#_Base_Pair vs. Mean Q Score")

plt.savefig(OUTPUT_FILE)
