//
//  Diskutil.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 05/05/2017.
//
//

#ifndef Diskutil_hpp
#define Diskutil_hpp

#include <vector>
#include <utility>
#include <iostream>
#include <boost/filesystem.hpp>
#include <unordered_map>
using namespace boost::filesystem;
using namespace std;

class Diskutil {
private:
public:


  /**
   * @brief Retrieves the paths of all items within a specified directory.
   *
   * This function iterates through the specified directory and collects the paths 
   * of all items within it. If the specified path does not exist or is not a directory, 
   * an error message is output to the console, and the program exits with a status of 1.
   *
   * @param parent The path of the directory to be scanned.
   * @return A vector of paths, each pointing to an item within the specified directory.
   */
  static vector<path> getItemsInDir(path parent);



  /**
   * @brief Checks if a given directory path points to a valid existing directory.
   *
   * This function verifies if the specified directory path exists and represents a directory
   * on the file system.
   *
   * @param dir The path of the directory to check.
   * @return True if the directory exists and is a valid directory; false otherwise.
   */
  static bool isFolderValid(path dir);



  /**
   * @brief Checks whether a specified file has a particular file extension.
   *
   * This function first checks if the path provided points to an existing regular file. 
   * If it does, it then checks whether the file has the specified extension.
   *
   * @param dir The path of the file to check.
   * @param ext The file extension to check for (including the leading dot, e.g., ".txt").
   * @return True if the file exists, is a regular file, and has the specified extension; false otherwise.
   */
  static bool hasFileExtension(path dir, string ext);



  /**
   * @brief Counts the occurrences of k-mers within a specified range of a buffer.
   *
   * This function iterates through the given buffer of characters, computing and counting
   * the occurrences of k-mers from the start to end positions. The count of each k-mer is stored
   * in a hash map where the key is an integer representation of the k-mer and the value is the
   * count of occurrences. The function skips newline characters and handles invalid characters by
   * avoiding the processing of k-mers containing them. Both the k-mer and its reverse complement
   * are considered, but only the lexicographically smaller one is counted.
   *
   * @param my_map A reference to an unordered_map to store the count of each k-mer.
   * @param kmer_size The size of the k-mers to be counted.
   * @param buffer A reference to a vector of characters containing the sequence data.
   * @param start The starting position in the buffer from where to begin counting k-mers.
   * @param end The ending position in the buffer until where to count k-mers.
   */
  static void countKmer(unordered_map<int, int>& my_map, int kmer_size, const vector<char>& buffer, const uint64_t& start, const uint64_t& end);



  /**
   * @brief Retrieves the paths of genome training data within a specified source directory.
   *
   * This function traverses through the specified source directory to find and organize genome training data.
   * Each class of data should reside in its own sub-directory under the source directory. The genome files
   * in these sub-directories should have the specified file extension.
   *
   * @param source_folder The directory path where the genome training data is stored.
   * @param extension The file extension of the genome files to look for.
   * @return A vector of tuples, where each tuple contains the class name, the path of the genome file, 
   * and the path of the sequence file (which is derived from the genome file name).
   */
  static vector<tuple<string, path, path> >
    getTrainingGenomePaths(path source_folder, string extension);



  /**
   * @brief Returns the size of the specified file in bytes.
   *
   * If the file does not exist, this function will return 0.
   *
   * @param file The path to the file whose size is to be retrieved.
   * @return The size of the file in bytes. Returns 0 if the file does not exist.
   */
  static size_t getFileSize(path file);

  static const string SEQUENCE_FILE_EXT; // = ".fasta"
  static const string SAVE_FILE_EXT; // = ".dat"

protected:
};

#endif /* Class_hpp */
