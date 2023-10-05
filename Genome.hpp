//
//  Genome.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 01/05/2017.
//
//

#ifndef Genome_hpp
#define Genome_hpp

#include <vector>
#include <utility>
#include <unordered_map>
#include <queue>
#include <mutex>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics.hpp>
#include "Diskutil.hpp"
using namespace std;
using namespace boost::accumulators;


template <class T>
class Class;

class Genome {
private:
public:
  typedef pair<double, Class<int>* > score;
  typedef priority_queue<score, vector<score>, less<score> > pqueue;

  /**
   * Instantiates a Genome instance.
   * @param  kmr_path Path to file containing serialized k-mers.
   * @param  sequence_path Path to file containing DNA sequence.
   */
  Genome(path kmr_path, path sequence_path);

  ~Genome();

  /**
   * Loads kmer counts from the file specified in the "kmr_path" instantiation
   * parameter. When triggered by multiple threads (on the same object), will
   * wait for the initial thread to finish loading data and then simply return
   * to all threads waiting.
   */
  void loadKmerCounts();

  /**
   * Reads FASTA sequence from disk using the file specified in the
   * "sequence_path" instantiation parameter.
   * NOT THREAD SAFE, at the moment.
   */
  void loadSequence();

  /**
   * Deletes the sequence and/or the kmer counts from memory.
   * Leaves the numerator queue intact.
   */
  void unload();

  /**
   * Computes optimal sequence alignment between the current instance and the
   * Genome supplied in the parameter. Using the Smith-Waterman algorithm.
   * @param  seq Genome instance to compare.
   * @return     The pair's alignment score (greater value is better).
   */
  long long int computeAlignmentScore(Genome* seq);

  /**
   * Takes in a Class instance, and returns the prior numerator; a score
   * correlated to the probability of the genome instance pertaining to the
   * supplied class.
   * @param  class Instance of Class object to use for the computation.
   * @return       Score correlated to probability of class membership.
   */
  double computeClassificationNumerator(Class<int>* cls);

  /**
   * Gets the character located at position "pos" in the DNA sequence.
   * @param  pos Integer, representing the position we want to query.
   * @return     The base associated to that position.
   */
  string charAt(int pos);

  /**
   * Getter for the FASTA sequence path.
   * @return     The path to the FASTA file.
   */
  string& getSequence();


  /**
   * Returns the path on disk of the kmr count file registered for this read.
   * @return The file's path.
   **/
  path getKmrPath();

  /**
   * @return     The length of the FASTA sequence associated to this genome.
   */
  int size();

  /**
   * @return     A hashmap containing the counts (values) associated to each
   *             kmer (key).
   */
  unordered_map<int, int>& getKmerCounts();


  /**
   * @brief Sets the k-mer counts for the genome and marks them as loaded.
   * 
   * This function updates the k-mer counts data member with the provided 
   * unordered map and sets the `kmersLoaded` flag to true. It's typically
   * called after the k-mer counts for the genome have been computed and 
   * are ready to be stored within the object.
   * 
   * @param _kmer_counts Pointer to an unordered_map containing k-mer counts 
   * where the key represents the k-mer (as an integer) and the value is its count.
   */
  void setKmerCounts(unordered_map<int, int>* _kmer_counts);


  /**
   * @brief Resets the k-mer counts status for the genome.
   * 
   * This function resets the `kmersLoaded` flag to false, indicating that the
   * k-mer counts for the genome are no longer available or loaded. Note that
   * it does not clear the actual k-mer counts data, only the status flag.
   */
  void resetKmerCounts();

  static bool STORE_ALL_NUMERATORS;
protected:
  bool kmersLoaded = false, sequenceLoaded = false;
  pqueue numerator;
  mutex numeratorAccess, loadKmersLock, loadSequenceLock;
  path sequence_path, kmr_path;
  string *sequence;
  unordered_map<int, int> *kmer_counts;
  long long int max(long long int a,
                    long long int b,
                    long long int c,
                    long long int d);

  double maximumNumerator;
  Class<int>* maximumNumeratorClass = NULL;
  
};

#include "Class.hpp"
#include "NB.hpp"

#endif /* Genome_hpp */
