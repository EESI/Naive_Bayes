# NBC++ / Naive Bayes C++
A Naive Bayes classifier for metagenomic data.  
Development version NOT FOR RELEASE, © Drexel University EESI Lab, 2017, last update 2024.  
Maintainer: Haozhe Duan and Alexandru Cristian, alex dot cristian at drexel dot edu  
Owner: Gail Rosen, gailr at ece dot drexel dot edu  

## Usage (For the more detailed instruction, please read "NBC_Manual.pdf"
Pick one of the execution [modes]: train to create or add new sequences to a savefile (the .dat files), or classify to use existing savefiles and classify reads. In both cases, the [source dir] will contain the DNA sequences / reads. The benchmark mode is currently unsupported.

When training, the tags will be expressed by placing each sequence into a subfolder of the source directory, by the name of its respective class. When classifying, make sure the folder only contains ONE file to be classified, and no other files. If there are multiple files, please concatenate them into one file.

For training, NBC takes the kmer count file as input. If the file is sequence file, please use kmer count tools to generate the kmer count file.

For classifying, by default, the program will output the max likelihood class for each read. If the -f option is specified, the program will output the log likelihoods for each class for each read, and the max likelihood class. If the sequence read contains no valid kmers (kmer = 0), the program will output "sequence contains no valid kmers" for the read.

```
./NB.run [mode: train/classify] [source dir] [options]

Generic options:
  -h [ --help ]         Print help message
  -v [ --version ]      Print version information

Allowed options:
  -s [ --savedir ] arg (=./NB_save)     Path to training data "save" folder (where the training data is to be or was saved)
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
In the current build of NBC++, the "-m" option refers to the max memory (in MBs) the program can use.

For classify mode, please allocate at least 2000 MBs.

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

The save files are named "[class name]-save.dat", where [class name] is the name of the class the savefile corresponds to.

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
