# Demultiplexing

**Final output from demultiplexing the full files can be seen in the FINAL_OUTPUT.txt file**

### Usage:
To run the demultiplexer you must run the command: **python nicks_demultiplexer.py**

Note: All files used with this program must be g-zipped (.gz)

##### Required flags:
* -r1: The file path to the forward biological read .fastq file
* -r2: The file path to the reverse biological read .fastq file
* -i1: The file path to the forward index .fastq file
* -i2: The file path to the reverse index .fastq file
* -c: The value you would like to set as the coverage cutoff value when throwing out low quality reads

#### Required packages:
* gzip
* argparse

**Example:**
```
python nicks_demultiplexer.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
```

#### OR:

Use the **run_demultiplexing.srun** script on talapas which should work for everyone since abosolute paths are used.

```
forward_read='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
reverse_read='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'
forward_index='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
reverse_index='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'


/usr/bin/time -v python nicks_demultiplexer.py -r1 $forward_read -r2 $reverse_read -i1 $forward_index -i2 $reverse_index
```

### Runtime:

To demultiplex the full files the script took a little bit longer than an hour to run. Here are the specs:

```
Command being timed: "python part2_demultiplexing.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
	User time (seconds): 5784.10
	System time (seconds): 74.86
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:38:32
	Maximum resident set size (kbytes): 126268
	Exit status: 0
```
