//
//  NB.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#ifndef NB_hpp
#define NB_hpp

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include <thread>
#include <mutex>
#include "Genome.hpp"
#include "Class.hpp"
#include "Diskutil.hpp"
#include <condition_variable>
#include <sys/resource.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <dirent.h>
#include <regex>
using namespace std;

template <class T>
class Class;

class NB {
private:

    /**
     * @brief Searches a specified directory for a particular class header file and a set of temporary files that need to be appended, 
     *        and sorts the temporary files based on a specified naming convention.
     *
     * This function scans a given directory for a class header file and a set of temporary files. The class header 
     * file should have a name that matches the pattern "<output_class_index>.clshd". The temporary files should 
     * have names that follow the format "<filter_number>_<sequence_number>.tmp", where <filter_number> matches the 
     * specified filter (class group) number, and <sequence_number> is a numeric value indicating the sequence order. 
     * The function collects and sorts the temporary files in ascending order based on their sequence numbers, and 
     * returns the path of the class header file corresponding to the list of temporary files.
     *
     * @param filter_number The class group number used to identify relevant temporary files.
     * @param file_names A vector that will be populated with the paths of the sorted temporary files.
     * @param directory The path of the directory to be scanned.
     * @return The path of the class header file, or an empty string if no such file is found.
     */
    string getAppendFiles(int filter_number, std::vector<std::string>& file_names, const std::string& directory);


    /**
     * @brief Scans a specified directory and groups the files that need to be concatenated based on their naming patterns.
     *
     * This function traverses a given directory searching for files that adhere to specific naming conventions:
     * 1. Files named "concat_[sequence_number].tmp".
     * 2. Files named "[class_number]_[sequence_number]_f.[extension]".
     *
     * These files are categorized into distinct groups based on their group numbers and saved in the provided `file_groups` map.
     * The key of this map signifies the group (sequence) number, and its value is a vector containing the file names from that particular group.
     * Within each group, the files are sorted as follows:
     * - "concat_" files are placed before other temp files.
     * - The other files are organized by the numbers present in their names.
     *
     * @param file_groups An unordered map where the key is the group number and the value is a vector holding file names from that group.
     * @param directory The directory path where the function will look for the concatenated files.
     */
    void getConcatFiles(std::unordered_map<size_t, std::vector<std::string>>& file_groups, const std::string& directory);


    /**
     * @brief Loads the max from a list of files and stores the values into a provided list.
     *
     * This function iterates over a given list of files (`max_files`), where each file contains a series of `double` values.
     * It reads these files by memory-mapping them, extracts their `double` values, and appends them into the `max_list` vector.
     *
     * Each file corresponds to a separate vector<double> within `max_list`.
     *
     * In case of any issues while opening or mapping a file, an error message is printed to `std::cerr`, and the program exits.
     *
     * @param  max_files A vector of pairs: <sequence_number, file_path>. This function processes files based on these paths.
     * @param  max_list A vector of vectors where each inner vector stores `double` values read from the corresponding file in `max_files`.
     */
    void loadMaxResults(vector<pair<int,string>>& max_files, vector<vector<double>>& max_list);


    /**
     * @brief Compares values from a `max_vector` with a `compare_vector`, updates the max values then write the max into a binary file.
     *
     * The function works in the following way:
     * - It iterates through the odd-indexed values of the `max_vector` (<index, max_value> pair).
     * - For each odd-indexed value (max_value) , it compares the value with the corresponding value in `compare_vector`.
     * - If the value from `max_vector` is greater and also non-positive (i.e., `<= 0`), it updates the `compare_vector` at that index and also updates the preceding even-indexed value (considered as the 'index').
     * - After the comparison and modifications are complete, the function writes the modified `compare_vector` to a binary file specified by the `filepath`.
     *
     * If there's any issue in opening the output file for writing, an error message is printed to `std::cerr` and program exits.
     *
     * @param  max_vector      A vector containing values to be compared. It has a structure where even-indexed elements represent 'index' and odd-indexed elements represent the value.
     * @param  compare_vector A vector containing values that will be compared against and possibly modified. Has the same structure as `max_vector`.
     * @param  filepath        The path to the output file where the modified `compare_vector` will be written.
     */
    void compareAndWriteMax(vector<double> max_vector, vector<double>& compare_vector, const string &filepath);


    /**
     * @brief Reads sequences of strings from a binary file and populates a provided vector with those strings.
     *
     * This function processes a given binary file that contains sequences of strings in a specific format. Each string 
     * in the file is stored with a preceding size value (of `size_t` type) which indicates the length of the string. 
     * The function reads the size value first and then reads the subsequent string data. Each extracted string is 
     * appended to the provided `seq_headers` vector. 
     *
     * The function continues reading strings from the file until the end of the file is reached. If there's any 
     * issue during the reading process (like reading the size or the string data), an error message is printed to 
     * `std::cerr` and then program exits.
     *
     * @param seq_headers   A vector that will be populated with the strings extracted from the binary file.
     * @param file_name     The path to the binary file to be read.
     */
    void getStringSeqHeaders(vector<string>& seq_headers, const string& file_name);


    /**
     * @brief Extracts class headers from a series of memory-mapped regions, updating header vectors and offsets accordingly.
     *
     * This function iterates over a collection of memory-mapped regions, each expected to begin with a size_t value
     * indicating the total number of headers present. Immediately following this value, each header's length (size_t)
     * and actual string data are sequentially stored.
     * 
     * The function processes each mapped region to retrieve all the class headers. These headers are appended to the 
     * `cls_headers` vector. Furthermore, it calculates the offset in each mapped region after reading all headers and 
     * stores these offsets in the `offsets` vector. The total number of headers for each mapped region is also saved 
     * in the `cols` vector.
     *
     * @param cls_headers  A vector that will be populated with class headers extracted from the memory-mapped regions.
     * @param maps         A vector containing pairs where the first element is a pointer to the start of a memory-mapped 
     *                          region and the second element is the size of that region.
     * @param offsets      A vector that will be populated with offsets indicating the start position in each mapped region after 
     *                          reading all headers.
     * @param cols         A vector that will be populated with the total number of headers for each mapped region.
     */
    void getStringClassHeaders(vector<string>& cls_headers, vector<pair<void*, size_t>>& maps, vector<size_t>& offsets, vector<size_t>& cols);


    /**
     * @brief Scans a predefined temporary directory for sequence header files, extracts sequence numbers, and populates a vector with file details.
     *
     * This function iterates over the files in the temp (working) directory searching for files with names 
     * that match the pattern `<class_number>_<sequence_number>_f.hd`. The function extracts the sequence number from 
     * each matched file name and pairs it with the file's full path. The resulting vector of pairs (`seq_header_files`) 
     * is sorted based on the sequence numbers in ascending order.
     * 
     * For example, if the directory contains the files `1_10_f.hd`, `1_5_f.hd`, and `2_1_f.hd`, the resulting vector will 
     * contain pairs (5, <full_path_to_1_5_f.hd>), (10, <full_path_to_1_10_f.hd>), and (1, <full_path_to_2_1_f.hd>).
     *
     * @param[out] seq_header_files  A vector that will be populated with pairs. Each pair consists of a sequence number 
     *                              and the full path to the corresponding file.
     */
    void getSeqHeaderFiles(vector<pair<int,string>>& seq_header_files);


    /**
     * @brief Scans the temp (working) directory for `.max` files, extracts their sequence numbers, and populates a vector with file details.
     *
     * This function iterates over the files in a predefined directory (`temp_dir`) searching for files with names 
     * that match the pattern `<sequence_number>.max`. The function extracts the number 
     * from each matched file name and pairs it with the file's full path. The resulting vector of pairs (`max_files`) 
     * is sorted based on the sequence numbers in ascending order.
     * 
     * For example, if the directory contains the files `5.max`, `10.max`, and `1.max`, the resulting vector will 
     * contain pairs (1, <full_path_to_1.max>), (5, <full_path_to_5.max>), and (10, <full_path_to_10.max>).
     *
     * @param[out] max_files  A vector that will be populated with pairs. Each pair consists of a sequence number 
     *                        and the full path to the corresponding `.max` file.
     */
    void getMaxFiles(vector<pair<int,string>>& max_files);


    /**
     * @brief Processes and outputs maximum results from given header files and associated max files to a CSV format.
     *
     * This function pairs sequence header files with their corresponding max files. For each pairing, the function 
     * processes the max file and pairs its entries with the appropriate sequence headers. The sequence headers are 
     * used to extract metadata and results are written to a CSV file. The resulting CSV file will contain each 
     * sequence header, the associated class string (extracted from the filename of the training genome), and the 
     * data value from the max file.
     * 
     * The CSV format is as follows:
     * - Sequence Header,Class String,Data Value
     * 
     * If a data value in the max file is negative, the entry "invalid kmer" is written instead of the class string.
     * 
     * For example, if the max file contains the values `[1, 0.5, -1, 0.3]` and the sequence headers are `["header1", "header2"]`, 
     * the resulting CSV might look like:
     * ```
     * header1,class_string1,0.5
     * header2,invalid kmer,
     * ```
     *
     * @param seq_header_files A vector of pairs containing sequence numbers and the full path to the corresponding sequence header file.
     * @param max_files        A vector of pairs containing sequence numbers and the full path to the corresponding max file.
     */
    void outputMaxResults(vector<pair<int,string>>& seq_header_files, vector<pair<int,string>>& max_files);


    /**
     * @brief Removes files from the temp (working) directory that have one of the given extensions.
     *
     * This function searches through the `temp_dir` directory and removes all files that have 
     * a file extension that matches one of the provided extensions in the `extensions` vector.
     * If a file fails to be removed, an error message is printed to the standard error stream.
     * 
     * Usage example:
     * ```
     * NB nb_instance;
     * nb_instance.removeFilesWithExtensions([".tmp", ".bak"]);
     * ```
     * This will remove all files in `temp_dir` that have the `.tmp` or `.bak` extensions.
     *
     * @param[in] extensions A vector containing the list of file extensions (including the leading dot) to search for and remove.
     */
    void removeFilesWithExtensions(const std::vector<std::string>& extensions);


    /**
     * @brief Finalizes processing of a batch of outputs, manages synchronization and appends results.
     *
     * This function ensures that all the queued up outputs are written and processed correctly. 
     * It synchronizes with the writing tasks, appends results to the relevant output files, and 
     * manages internal counters and flags. Depending on the value of `NB::OUTPUT_FULL_LOG_LIKELIHOOD`, 
     * this function either appends full log likelihoods or simply appends maximum likelihoods.
     *
     */
    void wrapUpBatch();
public:
    /**
     * Instantiates a NB instance.
     * @param _kmer_size The kmer size to use.
     * @param _save_dir  The root directory containing savefiles for all
     *                   classes.
     * @param _nthreads  The number of threads with which to simultaneously
     *                   process classes.
     */
    NB(int _kmer_size, path _save_dir, int _nthreads = 1);

    ~NB();

    /**
     * @brief Loads training data from saved files in `save_dir`.
     *
     * Scans the `save_dir` for files with the extension `Diskutil::SAVE_FILE_EXT`. 
     * For each valid file, it extracts class information and populates the class list.
     */
    void loadTrain();



    /**
     * @brief Loads training data from saved files in `save_dir`, sort them and record their file size.
     *
     * Scans the `save_dir` for files with `Diskutil::SAVE_FILE_EXT`. For each valid 
     * file, it adds the file to the training genomes and calculates the estimated 
     * size of training data for it. 
     *
     */
    void loadClassify();

    /**
     * @brief Directs all classes to save their data to disk.
     * Classes that haven't been loaded will be skipped.
     */
    void save();

    /**
     * @brief Getter for the path of the save directory.
     * @return A boost::path object pointing to the save directory.
     */
    path getSavedir();

    /**
     * @brief Getter for the kmer size to be used within NB.
     * @return An integer - the current kmer size.
     */
    int getKmerSize();

    /**
     * @brief Getter for the total number of threads to be used within NB.
     * @return An integer - the total number of threads.
     */
    int getThreadNumber();

    /** 
     * @brief Adds a class to the classifier.
     * @param cl The class to be added.
     */
    void addClass(Class<int>* cl);

    /**
     * @brief Looks up a class in the classifier by its ID.
     * @param cl_id Unique class ID.
     * @return A pointer to the class, or nullptr if the class is not found.
     */
    Class<int>* getClass(string cl_id);

    /**
     * @brief Adds a class to the update queue, instead of processing updates right
     * away.
     * @param cl Class to add to the queue.
     */
    void addClassToUpdateQueue(Class<int>* cl);

    /**
     * @brief Triggers all queued classes to update.
     */
    void processClassUpdates();


    /**
     * @brief Sets the buffer height for writing.
     * 
     * @param size The desired buffer height.
     */
    void setWriteBufferSize(size_t size);
    

    /**
     * @brief Writes the output data to CSV files. If the `OUTPUT_FULL_LOG_LIKELIHOOD` 
     *        flag is enabled and it's the first sequence index, it creates a header 
     *        file containing the names of the training genomes (.clshd). Otherwise, it writes 
     *        the output values row by row to temporary files.

     * 
     * @note Errors in file operations (e.g., file creation failure) lead to an 
     *       error message and program termination.
     */
    void writeToCSV();


    /**
     * @brief Processes the loading of classes in threads, ensuring that memory usage doesn't exceed the provided limit.
     * 
     * @param memoryLimit The maximum memory limit allowed for loading classes.
     * @param used_memory The current amount of memory used by the program.
     */
    void loadClassesThreadProcess(uint64_t &memoryLimit, uint64_t &used_memory);


    /**
     * @brief Loads the classes in a multithreaded manner while ensuring the memory usage doesn't exceed the provided limit.
     * 
     * @param memoryLimit The maximum memory limit allowed for loading classes.
     * 
     */
    void loadClasses(uint64_t &memoryLimit);


    /**
     * @brief Processes the unloading of classes in threads. It removes classes from memory and cleans up resources.
     * 
     */
    void unloadClassesThreadProcess();


    /**
     * @brief Unloads classes that have been previously loaded into memory. Clean up and print the result if all training genomes have been processed.
     * 
     * This function pushes all the loaded classes to a queue for processing, then uses multiple threads (`nthreads`)
     * to unload each class. After all classes have been unloaded, it wraps up the batch and 
     * updates the output class index. Additionally, if all training genomes have been processed,
     * it retrieves sequence header files and max files, and processes the results. Finally, it
     * removes files with specific extensions.
     *
     */
    void unloadClasses();


    /**
     * @brief Checks if all training genomes have been loaded.
     * 
     * This function checks whether the current index of genomes being loaded (`load_start_index`)
     * has reached or surpassed the total number of training genomes available. 
     * 
     * @return bool Returns `true` if all training genomes have been loaded; otherwise, returns `false`.
     *
     */
    bool hasLoadedAll();


    /**
     * @brief Queues a job for classification for 1 se  uence and handles buffer-related operations.
     * 
     * This function is responsible for:
     * - Queuing a sequence classification job.
     * - Managing memory buffers for outputs.
     * - Writing output to CSV when a specified buffer size threshold is reached.
     * 
     * @param buffer A reference to the buffer containing sequence data to be classified.
     * @param sequence_data A tuple containing data about the sequence:
     *     - The first value represents the start position of the sequence.
     *     - The second value represents the end position of the sequence.
     *     - The third value represents the total number of sequences processed.
     *
     */
    void queueClassifyJob(vector<char>& buffer, tuple<uint64_t, uint64_t, uint64_t>& sequence_data);


    /**
     * @brief Initializes and starts the threads responsible for classification.
     * 
     * This function:
     * - Initializes the output buffer.
     * - Sets the `job_done` flag to false, indicating that jobs are still pending.
     * - Starts multiple classification threads based on the number specified in `nthreads`.
     *
     */
    void startClassifyThreads();


    /**
     * @brief Waits for all classification threads to complete and clears the threads.
     * 
     * This function:
     * - Sets the `job_done` flag to true, indicating that no more jobs are pending.
     * - Notifies all waiting threads via the condition variable `jobUpdateStatus`.
     * - Waits for all classification threads to complete using the `join` method.
     * - Clears the list of threads.
     *
     */
    void joinClassifyThreads();


    /**
     * @brief Sets the temporary directory path for the instance.
     * 
     * @param _temp_dir The desired path to the temporary directory.
     *
     */
    void setTempDir(const string& temp_dir);


    /**
     * @brief Concatenates CSV files from the temporary directory and writes the concatenated result to the specified output name.
     * 
     * The concatenate the csv files of the same sequence number from the temporary directory, 
     * and writes the concatenated result for each group. 
     * After concatenation, the temporary CSV files used in the concatenation are deleted.
     * 
     * @param output_name The prefix for the output file where the concatenated result will be saved.
     *
     */
    void concatenateCSVs(string& output_name);


    /**
     * @brief Waits for all k-mers in sequences to be processed.
     * 
     * This function blocks until all queued k-mers in sequences have been processed. 
     * It's useful to ensure that all k-mers have been handled before moving on to the next stage.
     * 
     */
    void waitCalculatingAllKmers();


    /**
     * @brief Computes the average and maximum memory usage of classes.
     * 
     * This function calculates the average and maximum memory usage of the classes 
     * based on their sizes in bytes. The results are provided through reference parameters.
     * 
     * @param class_avg_bytes Output parameter for the average memory usage (in bytes) of the classes.
     * @param class_max_bytes Output parameter for the maximum memory usage (in bytes) of a single class.
     * 
     */
    void getClassMemoryUsage(uint64_t& class_avg_bytes, uint64_t& class_max_bytes);


    /**
     * @brief Extracts the headers from a given FASTA file and saves them in binary files.
     * 
     * This function reads a FASTA file to extract headers, which begin with the '>' character.
     * It then writes the headers into binary files in the specified temp directory.
     * The headers are grouped into several files based on a predefined limit (`output_max_row`).
     * 
     * @param inputFile Path to the input FASTA file.
     * 
     * @return Returns the total number of headers extracted from the FASTA file. If there's an error during processing, the program terminates.
     *
     */
    size_t extractHeader(string& inputFile);


    /**
     * @brief Sets the output prefix for the generated files.
     * 
     * This function assigns an output prefix which can be utilized by other methods in the class
     * when generating the output files. This ensures that the files have a consistent naming pattern
     * based on the provided prefix.
     * 
     * @param output_prefix A string representing the desired prefix for output files.
     *
     */
    void setOutputPrefix(const string& output_prefix);


    /**
     * @brief Sets the maximum output size for rows and columns.
     * 
     * This function specifies the maximum size for the rows and columns in the output files.
     * It allows control over how large each output file can be, especially useful when dealing 
     * with large datasets or when there's a need to split outputs into manageable chunks.
     * 
     * @param row The maximum number of rows for the output files.
     * @param col The maximum number of columns for the output files.
     *
     */
    void setMaxOutputSize(uint64_t row, uint64_t col);

    static enum Debug{
      NO_LOG,
      LOG_SOME,
      LOG_ALL
    } debug_flag;

    // buffer, seq_index, seq_start_index, seq_end_index
    typedef tuple<vector<char>*, uint64_t, uint64_t, uint64_t> classifyJob;

    static bool OUTPUT_FULL_LOG_LIKELIHOOD;

protected:
    unordered_map<string, Class<int>* > classes;
    unsigned int kmer_size, progress, total, nthreads, cls_size;
    path save_dir;
    queue<Class<int>* > classesToProcess;
    mutex classQueueAccess;
    mutex memoryCntAccess;
    mutex process_update_access;
    
    vector<pair<path, uint64_t>> training_genomes;
    uint64_t load_start_index;
    
    queue<classifyJob> classifyJobs;
    vector<thread> threads;
    bool job_done = false;
    mutex classify_job_lock;
    condition_variable jobUpdateStatus;
    mutex output_modify;

    mutex num_seq_kmer_counted_access;
    condition_variable num_seq_kmer_counted_cv;
    bool waiting_for_kmer_counting = false;

    uint64_t num_seq_kmer_processed = 0;
    uint64_t num_seq_kmer_processing = 0;

    vector<vector<double>> outputs;
    vector<vector<double>> outputs_write;
    
    size_t write_buffer_height;

    string temp_dir;
    uint64_t processed_seq_num = 0;
    uint64_t start_seq_index = 0;
    uint64_t output_class_index = 1;
    uint64_t output_seq_index;
    uint64_t output_max_row;
    uint64_t output_max_col;


    mutex output_mtx;
    condition_variable output_cv;

    condition_variable start_write_cv;
    condition_variable write_done_cv;
    bool finished_writing = true;

    int last_written_class_index = -1;
    string output_prefix;


    /**
     * @brief Handles the processing of classificating 1 sequence.
     * 
     * Each thread invokes this function to check if there are any available jobs in the `classifyJobs` queue.
     * Once a job is dequeued, k-mers are counted and a classification decision is made based on pre-loaded
     * training genome data. Results are stored in the `outputs` matrix.
     * 
     * Internal working:
     * 1. The function waits for a job to be available.
     * 2. Once a job is available, it pops the job, counts the k-mers using `Diskutil::countKmer`.
     * 3. It then checks the k-mers against the training genomes to classify the sequence.
     * 4. The classification results (both the predicted class and its score) are stored in the `outputs` matrix.
     * 5. If there's no k-mer data available, a default value is assigned to indicate an unclassified sequence.
     *
     */
    void classifyThreadController();


    /**
     * @brief Oversees the training process of classification models in a multi-threaded manner.
     * 
     * Each thread invokes this function to check if there are any class instances in the `classesToProcess` queue.
     * Once a class is dequeued, its genomes are loaded, processed, and then saved. This involves populating the class
     * with genome data (kmers), training it, and then serializing the results.
     * 
     * Internal working:
     * 1. The function waits for a class instance to be available in the `classesToProcess` queue.
     * 2. Once a class is available, it's dequeued and checked if it's already loaded.
     * 3. If the class hasn't been loaded yet, it's loaded, trained, and then saved.
     * 4. The class is then unloaded from memory and deleted to free up resources.
     *
     */
    void trainThreadController();


    /**
     * @brief Initializes the buffer used to store the classification results.
     * 
     * The buffer's size and structure depend on the NB::OUTPUT_FULL_LOG_LIKELIHOOD flag.
     * If the flag is set, the buffer will be structured to store the full log likelihoods for each class.
     * Otherwise, it will be sized to only store the max information.
     */
    void initializeOutputBuffer();


    /**
     * @brief Concatenates multiple CSV files column-wise and writes the result to an output CSV file.
     * 
     * The method uses memory-mapped files for efficient reading of large CSV files and then 
     * generates the output CSV file(s) by concatenating the columns from the input files. 
     * If the total columns exceed the allowed max columns for the output, multiple CSV files 
     * will be generated.
     * 
     * @param inputFiles A vector of file paths for the input CSV files.
     * @param output_prefix A prefix used to name the output CSV file(s).
     * @param sequence_num Used to generate unique output file names.
     */
    void concatenateCSVByColumns(const std::vector<std::string>& inputFiles, const std::string& output_prefix, const std::size_t& sequence_num);


    /**
     * @brief Appends content (row-wise) from multiple temporary files (binary format) to a single (or multiple) output file(s) (binary format), 
     *        based on a specified maximum file size.
     * 
     * The function uses memory-mapped I/O to efficiently process large files. The data is divided into columns, 
     * and the function ensures that data written to any output file is aligned by these columns. The class headers
     * for the columns are loaded from the given class_header_file.
     * 
     * Additionally, the function removes the processed temporary input files after they've been processed. 
     * 
     * @param inputFiles A list of paths to input files to be concatenated.
     * @param outputFile The prefix of the output file path(s).
     * @param class_header_file Path to a binary file containing the class headers.
     */
    void fullAppend(std::vector<std::string>& inputFiles, const std::string& outputFile, string class_header);


    /**
     * @brief Processes a list of input temporary (binary format) files to compare and update maximum values across all input files.
     * 
     * The function uses memory-mapped I/O for efficient processing of large files. For each row of data in the input file, 
     * the function compares the values to existing maximum values loaded from prior runs. 
     * If the current row values exceed the maximums, the function updates them using the `compareAndWriteMax` function.
     * 
     * The function also cleans up by removing the temporary input files after processing.
     * 
     * @param inputFiles A list of paths to input CSV files to be processed.
     */
    void maxAppend(std::vector<std::string>& inputFiles);
};

#endif /* NB_hpp */