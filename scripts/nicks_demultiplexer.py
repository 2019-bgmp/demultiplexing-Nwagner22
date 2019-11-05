#!/usr/bin/env python

# Define the problem:
# We have recieved 4 raw fastq output files from an illumina sequencer. It is possible
# for index hopping to take place during the sequencing process. Out goal is to get rid
# of the low coverage reads, the index hopped reads, and any read that has an index that is
# not one of the 24. Storing all of these full records in their respective files. There should be
# 52 files in total. 48 for the forward and reverse index files, 2 for undetermined, 2 for index hopped
# We would only care about the 48 files afterward because those are the reads that pass.

import gzip
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="A program to demultiplex raw illumina output")
    parser.add_argument("-r1", "--read1", help="use to specify the forward read file name as a string", required=True, type = str)
    parser.add_argument("-i1", "--index1", help="use to specify the forward index file name as a string", required=True, type = str)
    parser.add_argument("-i2", "--index2", help="use to specify reverse index file name as a string", required=True ,type = str)
    parser.add_argument("-r2", "--read2", help="use to specify the reverse read file name as a string", required=True, type = str)
    parser.add_argument("-c", "--covcut", help="use to specify the coverage cutoff value", default=30, type=int)
    return parser.parse_args()

args = get_args()               # calls get_args method from above assigns the arguments to args
FORWARD_READ = args.read1          # assigning forward read file name as string to global varible
FORWARD_INDEX = args.index1       # assigning forward index file name as string to global variable
REVERSE_INDEX = args.index2       # assigning reverse index file name as string to global varible
REVERSE_READ = args.read2          # assigning reverse read file name as string to global variable
COVERAGE_CUTOFF = args.covcut       #assigning coverage cutoff value as int to global varibale (defaults to 30)


ALL_INDEXES = [] # I use this to store the 24 indexes in from the file
DUAL_COUNTS_DICT = {} # each index as the key and the value INTEGER value that holds that amount that each one of these indexes is seen in a dual matched
HOPPED_COUNT = 0   # Initialize a counter variable to keep track of the number of index hopped records
UNDETERMINED_COUNT = 0    # Initialize a counter variable to keep track of the number of Undetermined records

REVERSE_COMP_DICT = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"} # Key: each of the 4 nucleotides   Value: the compliment to that nucleotide

INDEX_DICT_FW = {} # I would store each index as the key and the value as a filepointer to the respective forward file
INDEX_DICT_RV = {} # I would store each index as the key and the value as a filepointer to the respective reverse file
HOPPED_DICT = {'forward':open('../demultiplexed_output/hopped_R1.fastq', 'w'), 'reverse':open('../demultiplexed_output/hopped_R2.fastq', 'w')} # key: forward or reverse    # filepointer to respective file
UNDETERMINED_DICT = {'forward':open('../demultiplexed_output/undetermined_R1.fastq', 'w'), 'reverse':open('../demultiplexed_output/undetermined_R2.fastq', 'w')} # Key: forward or reverse    value: filepointer to respective file

# Read in all 24 of the indexes
with open("indexes.txt", 'r') as indexFile:
    titles = indexFile.readline().strip().split()
    for line in indexFile:    # loops through all of the indexes in the file indexes.txt
        line = line.strip().split()
        # Here i am creating filepointers for each index, forward and reverse, so that I can store them in their respective dictionaries so I can use them later on when looping
        current_fw_filepointer = open('../demultiplexed_output/R1_' + line[4] + '.fastq', 'w')
        current_rv_filepointer = open('../demultiplexed_output/R2_' + line[4] + '.fastq', 'w')

        # Here is where I am storing the keys as each of the 24 indexes
        # and values as the filepointers in their respective dictionaries
        INDEX_DICT_FW[line[4]] = current_fw_filepointer
        INDEX_DICT_RV[line[4]] = current_rv_filepointer

        DUAL_COUNTS_DICT[line[4]] = 0  # stores all of the indexes into the dictionary and sets their value to zero

        ALL_INDEXES.append(line[4])  # appending each of the indexes to my ALL_INDEXES list that holds them
# unit test:
assert len(ALL_INDEXES) == 24



# FUNCTION that converts a single character into a numerical phred score
def convert_phred(letter):                                # method to convert a single character into a phred score
    """Converts a single character into a phred score""" 
    return(ord(letter)-33)                                # ord returns the unicode code for any given unicode character. Then 33 is subracted off because that is the standard in order to avoid the first 33 non-visable characters
# unit test:
assert convert_phred("A") == 32


# FUNCTION to check if the quality of the given phred score meets my quality cutoff
def coverage_check(phred_letter):
    """ calls convert_phred and compares the numerical phred value to my chosen cutoff"""
    if (convert_phred(phred_letter) > COVERAGE_CUTOFF):
        return True         # then it passes the cutoff check
    else:
        return False        # the read does not pass and will be put into the R1 or R2 undetermined file
# unit test:
assert coverage_check('?') == False
assert coverage_check('G') == True


# FUNCTION to calculate the reverse compliment
def reverse_compliment(index):
    """Converts a given string to the reverse compliment string"""
    temp_converted = ""
    for letter in index:
        temp_converted = REVERSE_COMP_DICT[letter] + temp_converted
    return temp_converted
# unit test:
assert reverse_compliment("TATG") == "CATA"


# FUNCTION to check if an index is one of the 24 given
def check_index(current_index):
    """checks to see if the current_index is in my list ALL_INDEXES. Returns true if index is one of the valid 24"""
    if(current_index in ALL_INDEXES):       #current_index passes check
        return True
    else:
        return False    #does not pass and read will be put into the R1 or R2 undetermined file
# unit test:
assert check_index("GTAGCGTA") == True
assert check_index("GTACCGTA") == False

TOTAL_COUNT = 0
with gzip.open(FORWARD_READ,'rt') as R1, gzip.open(FORWARD_INDEX,'rt') as R2, gzip.open(REVERSE_INDEX,'rt') as R3, gzip.open(REVERSE_READ,'rt') as R4:
    while(True):
        R1_RECORD = []
        R2_RECORD = []
        R3_RECORD = []
        R4_RECORD = []


        TOTAL_COUNT += 1
        # We want to start off by grabbing the first record from each of the files
        # grab current R1 record and store as seperate values in the list
        R1_RECORD.append(R1.readline().strip())
        R1_RECORD.append(R1.readline().strip())
        R1_RECORD.append(R1.readline().strip())
        R1_RECORD.append(R1.readline().strip())

        # grab current R2 record and store as seperate values in the list
        R2_RECORD.append(R2.readline().strip())
        R2_RECORD.append(R2.readline().strip())
        R2_RECORD.append(R2.readline().strip())
        R2_RECORD.append(R2.readline().strip())

        # grab current R3 record and store as seperate values in the list
        R3_RECORD.append(R3.readline().strip())
        R3_RECORD.append(R3.readline().strip())
        R3_RECORD.append(R3.readline().strip())
        R3_RECORD.append(R3.readline().strip())

        # grab current R4 record and store as seperate values in the list
        R4_RECORD.append(R4.readline().strip())
        R4_RECORD.append(R4.readline().strip())
        R4_RECORD.append(R4.readline().strip())
        R4_RECORD.append(R4.readline().strip())

        # condition to break out of the loop
        # only will break out when the end of the file has been read and blank values are stored in the lists
        if(R1_RECORD[0] == ""):
            break

        # After we have access to the current record in each file, we then want to test for each of the four possibilities before we can determine if it is dual matched
        # 1. if either of the indexes have Ns, put the whole record in the respective undetermined files
        # 2. if the phred score of any of the index bases is below the coverage cutoff, put the whole record in the respective undetermined files
        # 3. if the forward index or the reverse-compliment of the reverse index do not match the 24 given indexes, put the whole record in the respective undetermined files
        # 4. if both are valid indexes, but are not the reverse compliment of eachother, then the whole record is put in the respective index hopped files
        
        unreversed_uncomplimented = reverse_compliment(R3_RECORD[1])

        # check 1 and 3:
        if(check_index(R2_RECORD[1]) != True or check_index(unreversed_uncomplimented) != True):
            # One of the indexes had either an N or it just did not match one of the 24 given indexes.
            
            # The forward sequence record will be appended to the undetermined_R1.fastq file with the header modified to contain both indices
            UNDETERMINED_DICT['forward'].write(R1_RECORD[0] + ' index:' + R2_RECORD[1] + ' RC:' + R3_RECORD[1] + '\n')
            UNDETERMINED_DICT['forward'].write(R1_RECORD[1] + '\n')
            UNDETERMINED_DICT['forward'].write(R1_RECORD[2] + '\n')
            UNDETERMINED_DICT['forward'].write(R1_RECORD[3] + '\n')
            #     The reverse sequence record will be appended to the undetermined_R2.fastq file with the header modified to contain both indices
            UNDETERMINED_DICT['reverse'].write(R4_RECORD[0] + ' index:' + R3_RECORD[1] + ' RC:' + R2_RECORD[1] + '\n')
            UNDETERMINED_DICT['reverse'].write(R4_RECORD[1] + '\n')
            UNDETERMINED_DICT['reverse'].write(R4_RECORD[2] + '\n')
            UNDETERMINED_DICT['reverse'].write(R4_RECORD[3] + '\n')

            UNDETERMINED_COUNT += 1  # iterate undetermined count

            continue # SKIP TO NEXT ITERATION OF WHILE LOOP
        else:
            #     Both indexes do not have Ns and match one of the 24 given indexes
            pass # DO NOTHING

        # check 2:
        for i in range(len(R2_RECORD[1])):
            if(coverage_check(R2_RECORD[1][i]) != True):
                # The forward sequence record will be appended to the undetermined_R1.fastq file with the header modified to contain both indices
                UNDETERMINED_DICT['forward'].write(R1_RECORD[0] + ' index:' + R2_RECORD[1] + ' RC:' + R3_RECORD[1] + '\n')
                UNDETERMINED_DICT['forward'].write(R1_RECORD[1] + '\n')
                UNDETERMINED_DICT['forward'].write(R1_RECORD[2] + '\n')
                UNDETERMINED_DICT['forward'].write(R1_RECORD[3] + '\n')
                # The reverse sequence record will be appended to the undetermined_R2.fastq file with the header modified to contain both indices
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[0] + ' index:' + R3_RECORD[1] + ' RC:' + R2_RECORD[1] + '\n')
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[1] + '\n')
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[2] + '\n')
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[3] + '\n')

                UNDETERMINED_COUNT += 1
                
                continue # SKIP TO NEXT ITERATION OF WHILE LOOP
            else:
                # The current base in this index passes the coverage check
                pass # DO NOTHING
            if(coverage_check(R3_RECORD[1][i]) != True):
                # The forward sequence record will be appended to the undetermined_R1.fastq file with the header modified to contain both indices
                UNDETERMINED_DICT['forward'].write(R1_RECORD[0] + ' index:' + R2_RECORD[1] + ' RC:' + R3_RECORD[1] + '\n')
                UNDETERMINED_DICT['forward'].write(R1_RECORD[1] + '\n')
                UNDETERMINED_DICT['forward'].write(R1_RECORD[2] + '\n')
                UNDETERMINED_DICT['forward'].write(R1_RECORD[3] + '\n')
                # The reverse sequence record will be appended to the undetermined_R2.fastq file with the header modified to contain both indices
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[0] + ' index:' + R3_RECORD[1] + ' RC:' + R2_RECORD[1] + '\n')
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[1] + '\n')
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[2] + '\n')
                UNDETERMINED_DICT['reverse'].write(R4_RECORD[3] + '\n')

                UNDETERMINED_COUNT += 1

                continue # SKIP TO NEXT ITERATION OF WHILE LOOP
            else:
                # The current base in this index passes the coverage check
                pass # DO NOTHING

        # check 4:
        if(R2_RECORD[1] != unreversed_uncomplimented):
            # this read has been index hopped

            # The forward sequence record will be appended to the hopped_R1.fastq file with the header modified to contain both indices
            HOPPED_DICT['forward'].write(R1_RECORD[0] + ' index:' + R2_RECORD[1] + ' RC:' + R3_RECORD[1] + '\n')
            HOPPED_DICT['forward'].write(R1_RECORD[1] + '\n')
            HOPPED_DICT['forward'].write(R1_RECORD[2] + '\n')
            HOPPED_DICT['forward'].write(R1_RECORD[3] + '\n') 
            # The reverse sequence record will be appended to the hopped_R2.fastq file with the header modified to contain both indices
            HOPPED_DICT['reverse'].write(R4_RECORD[0] + ' index:' + R3_RECORD[1] + ' RC:' + R2_RECORD[1] + '\n')
            HOPPED_DICT['reverse'].write(R4_RECORD[1] + '\n')
            HOPPED_DICT['reverse'].write(R4_RECORD[2] + '\n')
            HOPPED_DICT['reverse'].write(R4_RECORD[3] + '\n')

            HOPPED_COUNT += 1

            continue # SKIP TO NEXT ITERATION OF WHILE LOOP
        else:
            pass
            # There is no index hopping and therefore the two reads are dual matched

            # The forward sequence record will be appended to the dual matched file that goes with this index (ex. indx1_R1.fastq) with the header modified to contain both indices
            INDEX_DICT_FW[R2_RECORD[1]].write(R1_RECORD[0] + ' index:' + R2_RECORD[1] + ' RC:' + R3_RECORD[1] + '\n')
            INDEX_DICT_FW[R2_RECORD[1]].write(R1_RECORD[1] + '\n')
            INDEX_DICT_FW[R2_RECORD[1]].write(R1_RECORD[2] + '\n')
            INDEX_DICT_FW[R2_RECORD[1]].write(R1_RECORD[3] + '\n') 
            # The reverse sequence record will be appended to the ual matched file that goes with this index (ex. indx1_R2.fastq) with the header modified to contain both indices
            INDEX_DICT_RV[R2_RECORD[1]].write(R4_RECORD[0] + ' index:' + R3_RECORD[1] + ' RC:' + R2_RECORD[1] + '\n')
            INDEX_DICT_RV[R2_RECORD[1]].write(R4_RECORD[1] + '\n')
            INDEX_DICT_RV[R2_RECORD[1]].write(R4_RECORD[2] + '\n')
            INDEX_DICT_RV[R2_RECORD[1]].write(R4_RECORD[3] + '\n')

            DUAL_COUNTS_DICT[R2_RECORD[1]] += 1


# Close all filepointers
for key in INDEX_DICT_FW:
    INDEX_DICT_FW[key].close()
for key in INDEX_DICT_RV:
    INDEX_DICT_RV[key].close()
for key in HOPPED_DICT:
    HOPPED_DICT[key].close()
for key in UNDETERMINED_DICT:
    UNDETERMINED_DICT[key].close()

sum = 0
sum+=(HOPPED_COUNT/TOTAL_COUNT) * 100
sum+=(UNDETERMINED_COUNT/TOTAL_COUNT) * 100

percent_hopped = str((HOPPED_COUNT/TOTAL_COUNT) * 100)[0:4]
percent_undeterminded = str((UNDETERMINED_COUNT/TOTAL_COUNT) * 100)[0:4]

# Report values at the end
print("")
print("Total index hopped records:", HOPPED_COUNT)
print("Percentage of index hopped:", percent_hopped, '%')
print("")
print("Total undetermined records:", UNDETERMINED_COUNT)
print("Percentage of undetermined:", percent_undeterminded, '%')
print("")

for i in ALL_INDEXES:
    sum+=(DUAL_COUNTS_DICT[i]/TOTAL_COUNT) * 100
    print("Percentage from index '" + i + "':", str((DUAL_COUNTS_DICT[i]/TOTAL_COUNT) * 100)[0:5], '%')

print("")
print("Total dual matched records:", TOTAL_COUNT-HOPPED_COUNT-UNDETERMINED_COUNT)
print("Total percentage Dual matched:", str(sum - (float(percent_undeterminded) + float(percent_hopped)))[0:5], '%')
print("")
print("Overall total percentage:", str(sum)[0:5], '%')
print("")