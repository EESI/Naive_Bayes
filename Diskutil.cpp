//
//  Diskutil.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 05/05/2017.
//
//

#include "Diskutil.hpp"

const string Diskutil::SEQUENCE_FILE_EXT = ".fasta";
const string Diskutil::SAVE_FILE_EXT = ".dat";

vector<path> Diskutil::getItemsInDir(path parent) {
  vector<path> result;

  if (!exists(parent) || !is_directory(parent)) {
    cout << "Could not open directory path: "<< parent.native() <<"\n";
    exit(1);
  }
  directory_iterator iter(parent);
  for(;iter != directory_iterator(); iter++) {
    result.push_back(iter->path());
  }
  return result;
}

size_t Diskutil::getFileSize(path file) {
  if (!exists(file)) {
    return 0;
  }
  
  return file_size(file);
}

bool Diskutil::isFolderValid(path dir) {
  if (dir.native().compare(".") == 0 || dir.native().compare("..") == 0) {
    return false;
  }
  return exists(dir) && is_directory(dir);
}

bool Diskutil::hasFileExtension(path dir, string ext) {
  if (exists(dir) && is_regular_file(dir)) {
    return dir.has_extension() && (dir.extension().compare(ext) == 0);
  }else{
    return false;
  }
}

void Diskutil::countKmer(unordered_map<int, int>& my_map, int kmer_size, const vector<char>& buffer, const uint64_t& start, const uint64_t& end) {
  
  const int k = 2; // number of bits per base

  if (end - start + 1 < kmer_size) return;

  for (uint64_t i = start; i <= end - kmer_size + 1; i++) {
    
    if (buffer[i] == '\n') continue;

    int result = 0;
    int skippedChars = 0;
    bool validKmer = true;

    for (int j = 0; j < kmer_size && (i + j + skippedChars) < buffer.size(); j++) {
      char base = buffer[i + j + skippedChars];

      if ((j == kmer_size - 1) && base == '\n') { // checking if the last character in k-mer is newline
        validKmer = false;
        break;
      }
      
      while((i + j + skippedChars) < buffer.size() && base == '\n') {
        skippedChars++;
        base = buffer[i + j + skippedChars];
      }

      if ((i + j + skippedChars) >= buffer.size()) {
        validKmer = false;
        break;
      }

      result <<= k;

      switch (base) {
          case 'A': case 'a': result |= 0b00; break;
          case 'C': case 'c': result |= 0b01; break;
          case 'G': case 'g': result |= 0b10; break;
          case 'T': case 't': result |= 0b11; break;
          default: 
            validKmer = false;
            break; // handle invalid characters
      }

      if (!validKmer) break; 
    }

    // If not valid or not enough bases are read, continue to the next position
    if (!validKmer) continue;

    
    
    auto getComplementBits = [](int bits) -> int {
      // Bits for A=00, C=01, G=10, T=11
      switch(bits) {
        case 0b00: return 0b11;  // A -> T
        case 0b01: return 0b10;  // C -> G
        case 0b10: return 0b01;  // G -> C
        case 0b11: return 0b00;  // T -> A
      }
      return 0;
    };

    int tempResult = result;
    int resultRevComp = 0;
    for (int j = 0; j < kmer_size; ++j) {
      resultRevComp <<= k;
      resultRevComp |= getComplementBits(tempResult & 0b11);
      tempResult >>= k;
    }

    // Use the lexicographically smaller kmer
    if (result > resultRevComp) {
      result = resultRevComp;
    }

    

    my_map[result]++;
  }
}


vector<tuple<string, path, path> >
  Diskutil::getTrainingGenomePaths(path source_folder, string extension) {
    vector<tuple<string, path, path> > result;
    vector<path> class_paths = getItemsInDir(source_folder);
    for(vector<path>::iterator iter=class_paths.begin();
        iter != class_paths.end(); iter++) {
      if (!isFolderValid(*iter)) {
        continue;
      }

      string cls(iter->native());
      string cls_s;
      try{
        cls_s = cls.substr(cls.rfind(path::preferred_separator)+1);
      }catch(invalid_argument) {
        continue;
      }

      vector<path> genomes = getItemsInDir(*iter);
      for(vector<path>::iterator gen=genomes.begin();
          gen != genomes.end(); gen++) {
        if (!hasFileExtension(*gen, extension)) {
          continue;
        }

        string seqFile = gen->native();
        seqFile = seqFile.substr(0, seqFile.rfind("."));
        seqFile += SEQUENCE_FILE_EXT;

        result.push_back(make_tuple(cls_s, *gen, path(seqFile)));
      }
    }

    return result;
  }