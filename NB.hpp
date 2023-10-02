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
     * @brief Searches a specified directory for a particular class header file and a set of temporary files, 
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

    void getConcatFiles(std::unordered_map<size_t, std::vector<std::string>>& file_groups, const std::string& directory);

    void loadMaxResults(vector<pair<int,string>>& max_files, vector<vector<double>>& max_list);


    void compareAndWriteMax(vector<double> max_vector, vector<double>& compare_vector, const string &filepath);
    void getStringSeqHeaders(vector<string>& seq_headers, const string& file_name);
    void getStringClassHeaders(vector<string>& cls_headers, vector<pair<void*, size_t>>& maps, vector<size_t>& offsets, vector<size_t>& cols);
    void getSeqHeaderFiles(vector<pair<int,string>>& seq_header_files);
    void getMaxFiles(vector<pair<int,string>>& max_files);
    void outputMaxResults(vector<pair<int,string>>& seq_header_files, vector<pair<int,string>>& max_files);
    void removeFilesWithExtensions(const std::vector<std::string>& extensions);
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
     * Populates NB's list of classes by crawling through the root directory and
     * adding all savefiles available - without loading the class data itself.
     */
    void loadTrain();


    void loadClassify();

    /**
     * Directs all classes to save their data to disk.
     * Classes that haven't been loaded will be skipped.
     */
    void save();

    /**
     * Getter for the path of the save directory.
     * @return A boost::path object pointing to the save directory.
     */
    path getSavedir();

    /**
     * Getter for the kmer size to be used within NB.
     * @return An integer - the current kmer size.
     */
    int getKmerSize();

    /**
     * Getter for the total number of threads to be used within NB.
     * @return An integer - the total number of threads.
     */
    int getThreadNumber();

    /** MARKED-CHANGED
     * Adds a class to the classifier.
     * @param cl The class to be added.
     */
    void addClass(Class<int>* cl);

    /**
     * Looks up a class in the classifier by its ID.
     * @param cl_id Unique class ID.
     */
    Class<int>* getClass(string cl_id);

    /**
     * Adds a class to the update queue, instead of processing updates right
     * away.
     * @param cl Class to add to the queue.
     */
    void addClassToUpdateQueue(Class<int>* cl);

    /**
     * Triggers all queued classes to update.
     */
    void processClassUpdates();

    void setWriteBufferSize(size_t size);
    
    void writeToCSV();

    void loadClassesThreadProcess(uint64_t &memoryLimit, uint64_t &used_memory);

    void loadClasses(uint64_t &memoryLimit);

    void unloadClassesThreadProcess();

    void unloadClasses();

    bool hasLoadedAll();

    void queueClassifyJob(vector<char>& buffer, tuple<uint64_t, uint64_t, uint64_t>& sequence_data);

    void startClassifyThreads();

    void joinClassifyThreads();

    void setTempDir(const string& temp_dir);

    void concatenateCSVs(string& output_name);

    void waitCalculatingAllKmers();

    void getClassMemoryUsage(uint64_t& class_avg_bytes, uint64_t& class_max_bytes);

    size_t extractHeader(string& inputFile);

    void setOutputPrefix(const string& output_prefix);

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

    void classifyThreadController();
    void trainThreadController();
    void initializeOutputBuffer();
    void concatenateCSVByColumns(const std::vector<std::string>& inputFiles, const std::string& output_prefix, const std::size_t& sequence_num);
    void fullAppend(std::vector<std::string>& inputFiles, const std::string& outputFile, string class_header);
    void maxAppend(std::vector<std::string>& inputFiles);
};

#endif /* NB_hpp */