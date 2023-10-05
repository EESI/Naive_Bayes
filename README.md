# NBC++ / Naive Bayes C++
A Naive Bayes classifier for metagenomic data.  
Development version NOT FOR RELEASE, © Drexel University EESI Lab, 2017.  
Maintainer: Alexandru Cristian, alex dot cristian at drexel dot edu  
Owner: Gail Rosen, gailr at ece dot drexel dot edu  

## Usage
Pick one of the execution [modes]: train to create or add new sequences to a savefile, or classify to use existing savefiles and classify reads. In both cases, the [source dir] will contain the DNA sequences / reads. The benchmark mode is currently unsupported.

When training, the tags will be expressed by placing each sequence into a subfolder of the source directory, by the name of its respective class. The same applies when classifying, since NBC++ currently checks at the end of its execution how many picks were correct, providing the user with its accuracy.

If an unsupervised run is desired, for now the solution is to simply group all reads into one folder with an arbitrary tag.
```
./NB.run [mode: train/classify/benchmark] [source dir] [options]

Generic options:
  -h [ --help ]         Print help message
  -v [ --version ]      Print version information

Allowed options:
  -s [ --savedir ] arg (=./NB_save)     Path to save folder
  -k [ --kmersize ] arg (=6)            Kmer size used in count files
  -m [ --memlimit ] arg (=0)            Cap memory use to a predefined value (MBs).
  -t [ --nthreads ] arg (=1)            Number of threads to spawn, 1 by default
  -e [ --ext ] arg (=.kmr)              Extension of kmer count files, ".kmr" by default
  -n [ --nbatch ] arg (=0)              Number of genomes to load at one time in memory for training, all at once by default
  -f [ --full_result ]                  Print log likelihoods for class for every classified read.
  -o [ --output_prefix ] arg (=log_likelihood)  Output prefix.
  -d [ --temp-dir ] arg (=/tmp)         Temporary (working) directory path
  -r [ --row ] arg (=150000)            Maximum number of rows (classify mode)
  -c [ --col ] arg (=100000)            Maximum number of columns (classify mode)


```

### -m : Using the --memlimit memory cap option
In the current build of NBC++, the "-m" option refers strictly to the memory dedicated to loading the DNA reads/sequences into memory (depending on the mode used, train/classify).

Additional considerations must be taken when this option is combined with NBC++'s multithreading capability. For example, if using 15-mers, the maximum theoretical savefile size would be of 16GBs. In our experimental testing on our NCBI dataset, the maximum size approaches 4GBs for 15-mers. NBC++ will store in memory a number of classes equivalent to that of running threads, and their size is not accounted for by "--memlimit".

As a result, one must manually calculate the space available to load DNA reads, while also leaving some space to load the actual class savefiles. NBC++ will automatically adjust the number of DNA reads loaded to best fill the space allocated by the "--memlimit" option.
The default value is -1, which disables this option.

### -n : Using the --nbatch read cap option
This option fulfills the same function that --memlimit provides, the exception being that, in this case, a static cap is specified on the number of DNA reads loaded at one time into memory.

Note that while the number of reads stays the same, their size and that of their kmer count may vary. If a cap on memory usage is desired, use --memlimit instead.

The default value is -1, which disables this option.

### -t : Using the --nthreads option
This option specifies the number of concurrent threads to run. Each thread will compute data for a separate class, so each thread will load one savefile. DNA reads/sequences will be shared across all threads.

The default value is 1, which will run the program in single threaded mode.
Note: For classification, the number of threads should be at least 2.

### -e : Using the --ext option
This option specifies the extension of the files containing the kmer counts for each DNA read/sequence. NBC++ does not include a kmer counter by design; this part of the process should be taken care of by a separate module of the user's choosing. However, we do supply a kmer generating script that computes counts using Jellyfish 2 for all FASTA files in a given directory.

The default value is ".kmr", which will make the program look for files like "MyCount.kmr".

### -k : Using the --kmersize option
A different kmer size than the default may be selected by modifying this option.

The default value is 6, which will make the program output and expect 6-mers.

### -s : Using the --savedir option
This option specifies a target directory in which to save or load training data, depending on the mode. Each savefile in this directory corresponds to one class.

The default value is "./NB_save", which will direct the program to search for a folder named NB_save in its current directory.

## Project dependencies
This release requires the following dependencies to be available on the host machine:
- Grinder (https://sourceforge.net/projects/biogrinder/), tested on v. 0.5.4
- Jellyfish (http://www.genome.umd.edu/jellyfish.html), tested on v. 2.1.3

## Bug reports, feature requests
Please post an issue on our repo.

## Contributions and pull requests
...are welcome!

## Building on a local machine
First, install the dependencies. On Ubuntu, it suffices to run:
```
sudo apt install make g++ libboost-all-dev
```

Then, to build on your local machine, simply run:
```
make
```
Should you wish to debug a crash, a debug-friendly binary can be produced by building the "debug" Make target.
The experiment instructions below are written for Proteus runs, but they can be used for local runs as well - just make sure to compile in the way outlined above. Also, we'd recommend at least 16GBs for a relatively smooth experience - that should be able to allow for two cores to work concurrently. SSDs are highly recommended due to NBC++'s heavy use of disk I/O.


## Experiment Instruction (How to train and test the model):
1. Compile the NB classifier on Proteus (boost/openmpi/gcc is the dependency):
```
#! /bin/bash

module load boost/openmpi/gcc/64/1.57.0


make proteus
```
2. Train NB classifier:
change the variable `EXEC_ROOT` in nb-train.bash to the path to `NB.run` (the compiled NB classifier binary). Then run command:
```
nb-train.bash [path-to-data] [fold index] [total_num_of_folds] [K-mer size]
```
for example, let's assume we train NB classifier with 5-fold cross validation, the path to training data is `~/genomes/` (there are 5 directory in it, `fold1/` to `fold5`) and k-mer size is 15, then there will be 5 independent models. The command to train NB classifier and get the first model is:
```
nb-train.bash ~/genomes/ 1 5 15
```
Similarly, the third model out of five in 5-fold cross validation can be obtained by running:
```
nb-train.bash ~/genomes/ 3 5 15
```
3. Testing NB classifier:
change the variable `EXEC_ROOT` in nb-classify.bash to the path to `NB.run` (the compiled NB classifier binary). Then run command:
```
nb-classify.bash [path-to-model] [path-to-testing_data] [fold index] [K-mer size]
```
for example, let's assume we have trained NB classifier based on the configurations we described above. Then the 5 trained models are stored in the same directory as your data (`~/genomes/`). Hence, the path to model is `~/genomes/` too and k-mer size is 15, then there will be 5 independent models (`save_1/` to `save_5`). The command to classify reads in testing dataset (path: `~/ncbi_reads_noerr`) using the first model is:
```
nb-classify.bash ~/genomes/ ~/ncbi_reads_noerr 1 15
```
Similarly, classification results using the third model out of five in 5-fold cross validation can be obtained by running:
```
nb-classify.bash ~/genomes/ ~/ncbi_reads_noerr 3 15
```
The program will print out the classification progress as well as the results in stdout.

4. Experimental dataset (can be found in Proteus group folder):

Train dataset - ready on Dec. 15th 2018 (`~/genomes` directory in the example above):
```
/mnt/HA/groups/rosenGrp/zz374/INBC_latest/genomes
```
Test dataset - ready on Dec. 15th 2018  (`~/ncbi_reads_noerr` directory in the example above):
```
/mnt/HA/groups/rosenGrp/zz374/INBC_latest/ncbi_reads_noerr
```

