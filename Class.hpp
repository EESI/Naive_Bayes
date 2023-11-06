//
//  Class.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#ifndef Class_hpp
#define Class_hpp

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <queue>
#include <utility>
#include <cmath>
#include <cassert>
#include <cstring>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <sys/resource.h>
using namespace boost::filesystem;
using namespace std;

class Genome;

template <class T>
class Class {
private:
public:

    /**
     * Instantiates a Class instance.
     * @param id_       The unique identifier for this class.
     * @param kmer_size The kmer size used to create the savefile.
     * @param savefile  The path to the savefile on disk.
     */
    Class(string id_, int kmer_size, path savefile);

    ~Class();

    /**
     * Getter for the class ID.
     * @return A string containing the class unique ID.
     */
    string getId();

    /**
     * Prints the current instance in machine-readable format.
     * @return The serialized data of this class, typically contained in
     *         a savefile.
     */
    void serialize(std::ofstream& out);

    /**
     * Reads a savefile and reconstructs a class from machine-readable format.
     * @param in A stream pointing to the savefile.
     */
    void deserialize(std::ifstream& in);

    /**
     * Readies the instance for use and loads data from specified savefile.
     * @param source_path The path to the savefile on disk.
     */
    void load(path source_path);

    /**
     * Readies the instance for use and loads data from the set savefile.
     */
    void load();

    /**
     * Releases resources pertaining to kmer counts and any loaded DNA
     * sequences. These must be reloaded if they are to be used after this
     * function call.
     */
    void unload();

    /**
     * Writes the serialized class to a savefile.
     * @param destination_path The path to the savefile on disk.
     */
    void save(path destination_path);

    /**
     * Writes the serialized class to the set savefile.
     */
    void save();

    /**
     * Compute numerators for all genomes awaiting classification.
     * @param genomes Array of genomes for which to compute numerator.
     */
    void computeBatchNumerators(vector<Genome*> genomes);

    /**
     * Adds a genome to the class, training it to fit its pattern.
     * @param genome        The genome to be added to the class.
     * @param confidence_lg Our confidence in the genome's class membership.
     */
    void addGenome(Genome* genome, double confidence_lg = 0.0);

    /**
     * Queues a genome to be added to this class (with confidence = 1).
     * @param genome The genome to be trained on.
     */
    void queueGenome(Genome* genome);

    /**
     * Train this class on all of its queued genomes, waiting to be added.
     */
    void addGenomesInQueue();

    /**
     * Getter for the NGenomes_lg Naive Bayes parameter.
     * @return The log of the number of genomes present in this class.
     */
    double getNGenomes_lg();

    /**
     * Getter for the FreqCount_lg Naive Bayes parameter.
     * @param  feature_ The feature (kmer) for which we're looking up the count.
     * @return          The (smoothed and logged) count.
     */
    double getFreqCount_lg(T feature_);

    /**
     * Getter for the SumFreq_lg Naive Bayes parameter.
     * @return The total count of all kmers, smoothed and logged.
     */
    double getSumFreq_lg();

    /**
     * Gets a path to the savefile of this class.
     * @return Path to savefile.
     */
    path getSavefilePath();

    /**
     * @return True if the class is ready for use (has been initialized /
     *         loaded from disk).
     */
    bool loaded();

  /**
   * @brief Computes the size in bytes of a class instance and its associated data.
   * 
   * This function calculates the memory consumption of the current class instance. 
   * The calculated size is determined based on the size of the `Class<T>` object itself and 
   * the size of its associated data (`freqcnt_lg` and its elements).
   * 
   * @return The total size in bytes of the class instance.
   */
    size_t getClassSizeInBytes();
    

  /**
   * @brief Estimates the memory consumption in bytes for a hypothetical class instance based on the number of elements.
   * 
   * This function provides an estimation of how much memory a class instance would consume if it held 
   * a given number of elements. It calculates the size based on:
   * 1. Memory occupied by the keys (of type T).
   * 2. Memory occupied by the values (of type double_wflag).
   * 3. Memory occupied by the buckets containing the key-value pairs, with a bucketing factor of 1.5.
   * 
   * @param num_elements The number of key-value elements we are estimating memory for.
   * @return The estimated total size in bytes.
   */
    static size_t getEstimatedClassBytes(size_t& num_elements);



  /**
   * @brief Calculates the size in bytes of a single map element consisting of a key-value pair.
   * 
   * This function calculates the total memory consumed by an individual map element. The size is based on:
   * 1. Memory occupied by the key (of type T).
   * 2. Memory occupied by the value (of type double).
   * 
   * @return The combined size in bytes of both key and value.
   */
    static size_t getMapElementSize();

protected:
  bool isLoaded = false;
  string id;
  path savefile;
  vector<Genome*> genomes;
  queue<Genome*> queuedGenomes;
    
  typedef pair<double, bool> double_wflag;

  // Log data cached for classifying
  double_wflag ngenomes_lg, sumfreq_lg;
  unordered_map<T, double_wflag> *freqcnt_lg;

  // Plain data, for updates
  int ngenomes;
  long long int sumfreq;
  unordered_map<T, int> *freqcnt;

  void addNGenomes_lg(double ngenomes_);
  void addSumfreq_lg(double sumfreq_);
  void addFreqCount_lg(T feature_, double count_);
  void addNGenomes(int ngenomes_);
  void addSumfreq(long long int sumfreq_);
  void addFreqCount(T feature_, int count_);
  bool existsInFreqCount(T feature_);

  unordered_map<T, double_wflag>& getFreqcnt_lg();
  unordered_map<T, int>& getFreqcnt();

  static double logAdd(vector<double> exponents);

};

#include "Class.cpp"

#endif /* Class_hpp */
