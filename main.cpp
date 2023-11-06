//
//  main.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//
#include "NB.hpp"
#include "Diskutil.hpp"
#include <cstring>
#include <cstdint>
#include <boost/program_options.hpp>

using namespace std;
namespace p_opt = boost::program_options;

string usageMsg("./NB.run [mode: train/classify/benchmark] [source dir] [options]\n");
const int CONF_THRESHOLD = 100000;
const int DEF_KMER_SIZE = 6;
const string GEN_SEQ_EXTEN(".fna");
const string KMER_EXTEN(".kmr");
const string PROG_VER("NB v. 0.1.5a-dev.");
const string DEF_SAVEDIR("./NB_save");
const string KMER_NUM_FILE("kmer.in");
uint64_t max_row;
uint64_t max_col;

unordered_map<string, vector<double>* > confidence_list;
vector<path> read_filenames;
string result_file("");
string tempDir("");
bool full_result = false;

void trainNB(NB &nb, path srcdir, string extension, unsigned int nbatch,
             uint64_t memoryLimit){
              
  unsigned int count = 1, counter = 0;
  uint64_t usedMemory = 0;
  vector<tuple<string, path, path> > result =
    Diskutil::getTrainingGenomePaths(srcdir, extension);
  string cls_s="-1"; Class<int> *current = NULL;
  for(vector<tuple<string, path, path> >::iterator iter = result.begin();
      iter != result.end(); iter++, counter++){
        unsigned int savefileSize = 0;
        bool loadedNewClass = false;
        
        // FASTA files not needed for training, so just add up the kmr file size
        size_t genomeSize = Diskutil::getFileSize(get<1>(*iter));

        if(cls_s.compare(get<0>(*iter)) != 0){

          cls_s = get<0>(*iter);

          current = nb.getClass(cls_s);

          if(current == NULL){
            path save_file = path(nb.getSavedir().native()
                                  + path::preferred_separator
                                  + cls_s + "-save.dat");
            current = new Class<int>(cls_s,
                                     nb.getKmerSize(),
                                     save_file);

            nb.addClass(current);
          }

          savefileSize = Diskutil::getFileSize(current->getSavefilePath());
          loadedNewClass = true;
        }

        if(memoryLimit != 0 && usedMemory + genomeSize + savefileSize> memoryLimit){
          nb.processClassUpdates();
          usedMemory = 0;
          cls_s = "-1"; // This will force the next iteration to add this savefile's size again
        }

        Genome *genome = new Genome(get<1>(*iter), get<2>(*iter));
        current->queueGenome(genome);
        if(loadedNewClass){
          nb.addClassToUpdateQueue(current);
          usedMemory += savefileSize;
        }
        usedMemory += genomeSize;

        if(nbatch != 0 && counter % nbatch == 0){
          nb.processClassUpdates();
        }
        
      }
      
      nb.processClassUpdates();
}

// Find the optimal memory distribution for the memory allocation for loading classes and sequences.
void allocateMemoryDistribution(NB& nb, string& input_file, uint64_t& total_memory_limit, uint64_t& class_memory_limit, uint64_t& buffer_memory_limit, size_t& total_seq){
  
  const double DYNAMIC_LIST_CAPACITY_CONSTANT = 2;
  const size_t SEQ_PER_THREAD_INPUT_BUFFER = 200;
  const size_t SEQ_PER_THREAD_OUTPUT_BUFFER = 10 * SEQ_PER_THREAD_INPUT_BUFFER;
  const int MAX_SEQ_LENGTH = 2000000;

  size_t class_avg_bytes = 0;
  size_t class_max_bytes = 0;
  nb.getClassMemoryUsage(class_avg_bytes, class_max_bytes);

  size_t filesize = Diskutil::getFileSize(input_file);
  size_t estimate_seq_avg_length = filesize / total_seq;

  if(estimate_seq_avg_length > 1.5 * MAX_SEQ_LENGTH){
    cout << "Error: estimated average sequence length is larger than the max 2MB length." << endl;
    exit(1);
  }

  size_t num_classes_load = 0;
  buffer_memory_limit = (estimate_seq_avg_length * SEQ_PER_THREAD_INPUT_BUFFER * nb.getThreadNumber());
  if(buffer_memory_limit > filesize){
    buffer_memory_limit = filesize;
  }

  size_t remaining_memory = total_memory_limit;
  size_t substract_part = buffer_memory_limit + (2 * sizeof(double) * SEQ_PER_THREAD_OUTPUT_BUFFER * nb.getThreadNumber());
  
  remaining_memory -= substract_part;

  if(remaining_memory > total_memory_limit){
    cout << "Error: not enough memory. Memory needs to be greater than: " << substract_part << " bytes." << endl;
    exit(1);
  }

  if(full_result){
    num_classes_load = remaining_memory / 
            (sizeof(double) * SEQ_PER_THREAD_OUTPUT_BUFFER * nb.getThreadNumber() + class_avg_bytes);
  }else{
    num_classes_load = remaining_memory / (class_avg_bytes);
  }

  class_memory_limit = num_classes_load * class_avg_bytes;
  
  if(class_memory_limit < class_max_bytes){
    cout << "Error: Not enough memory to load the largest class: " << class_max_bytes << " bytes." << endl;
    exit(1);
  }

  size_t write_buffer_size = SEQ_PER_THREAD_OUTPUT_BUFFER * nb.getThreadNumber();
  nb.setWriteBufferSize(write_buffer_size);

  //check if can load entire output into memory
  const int SINGLE_CELL_LENGTH = 10;
  size_t out_file_size = 0;

  if(total_seq > max_row){
    out_file_size = SINGLE_CELL_LENGTH * max_row * sizeof(char);
  }else{
    out_file_size = SINGLE_CELL_LENGTH * total_seq * sizeof(char);
  }
  
  if(full_result){
    if(num_classes_load > max_col){
      out_file_size *= max_col;
    }else{
      out_file_size *= num_classes_load;
    }
  }else{
    out_file_size *= 2;
  }

  if(out_file_size > total_memory_limit){
    cout << "Error: not enough memory to load output file into memory. Please increase memory limit. Estimated output file size: " << out_file_size << endl;
    exit(1);
  }
}

void classifyNB(NB &nb, path srcdir, string &extension, unsigned int &nbatch,
                uint64_t &memory_limit, string output_file, const bool& full_result, const string& temp_dir){

  if (!exists(temp_dir)) {
    if (!create_directory(temp_dir)) {
      std::cout << "Failed To Create Temporary Directory: " << temp_dir << std::endl;
      exit(1);
    }
  }

  if(nb.getThreadNumber() < 1){
    cout<<"Error: at least 2 threads are required for classification.\n";
    exit(1);
  }
  
  NB::OUTPUT_FULL_LOG_LIKELIHOOD = full_result;

  nb.setTempDir(temp_dir);
  nb.setOutputPrefix(output_file);

  string inputFile = Diskutil::getItemsInDir(srcdir)[0].native();
  size_t total_seq = nb.extractHeader(inputFile);

  uint64_t class_memory_limit, buffer_memory_limit = 0;
  allocateMemoryDistribution(nb, inputFile, memory_limit, class_memory_limit, buffer_memory_limit, total_seq);

  std::vector<char> buffer(buffer_memory_limit, 0);
  
  thread writer(&NB::writeToCSV, &nb);

  while(!nb.hasLoadedAll()){
    std::ifstream file(inputFile);

    // Check if the file is opened successfully.
    if (!file.is_open()) {
      std::cout << "Failed to open file." << std::endl;
      exit(1);
    }

    unsigned int buffer_read_start_index = 0;
    uint64_t sequence_number = 1;
    vector< tuple<uint64_t, uint64_t, uint64_t> > sequences_position;

    std::streamsize bytesRead;

    nb.loadClasses(class_memory_limit);
    nb.startClassifyThreads();
    
    bool in_sequence = false;
    size_t sequence_start = 0;

    while (file) {
      file.read(buffer.data() + buffer_read_start_index, buffer_memory_limit - buffer_read_start_index);
      bytesRead = file.gcount();

      if (bytesRead == 0) {
          break;
      }

      size_t buffer_end = bytesRead + buffer_read_start_index;

      for (size_t i = 0; i < buffer_end; ++i) {
          char c = buffer[i];

          if (c == '>') { // Starting a new header
              if (in_sequence) {
                  sequences_position.emplace_back(sequence_start, i - 1, sequence_number++);
                  nb.queueClassifyJob(buffer, sequences_position.back());
                  sequences_position.clear();
              }
              in_sequence = false; // Currently reading a header, not a sequence
              continue;
          }

          if (!in_sequence && c == '\n') {
              // We've reached the end of the header and the start of the sequence
              sequence_start = i + 1;
              in_sequence = true;
          }
      }

      // Wait for all k-mers to be calculated before moving to the next buffer
      nb.waitCalculatingAllKmers();

      // If buffer ends with an incomplete sequence, adjust for next read
      if (in_sequence) {
          size_t leftover = buffer_end - sequence_start;
          for (size_t i = 0; i < leftover; ++i) {
              buffer[i] = buffer[sequence_start + i];
          }
          sequence_start = 0;
          buffer_read_start_index = leftover;
      } else {
          buffer_read_start_index = 0;
      }
    }

    // Process any remaining sequence in the buffer
    if (in_sequence && buffer_read_start_index > 0) {
        sequences_position.emplace_back(0, buffer_read_start_index - 1, sequence_number++);
        nb.queueClassifyJob(buffer, sequences_position.back());
    }

    file.close();
    nb.joinClassifyThreads();
    nb.unloadClasses();
  }

  writer.join();
}

int main(int argc, char* argv[]){
  unsigned int nbatch, nthreads, kmersize;
  uint64_t memLimit;
  string kmer_ext, srcdir, mode, savedir;

  p_opt::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "Print help message")
    ("version,v", "Print version information");

  p_opt::options_description hidden("Hidden options");
  hidden.add_options()
    ("mode", "Sets mode of program, train or classify.")
    ("srcdir", "Path to source folder");

  p_opt::options_description visible("Allowed options");
  visible.add_options()
    ("savedir,s", p_opt::value<string>(&savedir)->default_value(DEF_SAVEDIR),
                  "Path to save folder")
    ("kmersize,k", p_opt::value<unsigned int>(&kmersize)->default_value(DEF_KMER_SIZE),
                   "Kmer size used in count files")
    ("memlimit,m", p_opt::value<uint64_t>(&memLimit)->default_value(0),
                   "Cap memory use to a predefined value (MBs).")
    ("nthreads,t", p_opt::value<unsigned int>(&nthreads)->default_value(1),
                 "Number of threads to spawn, 1 by default")
    ("ext,e", p_opt::value<string>(&kmer_ext)->default_value(KMER_EXTEN),
            "Extension of kmer count files, \".kmr\" by default")
    ("nbatch,n", p_opt::value<unsigned int>(&nbatch)->default_value(0),
               "Number of genomes to load at one time in memory for training, \
all at once by default")
    ("full_result,f", p_opt::bool_switch(&full_result)->default_value(false),
               "Print log likelihoods for class for every classified read.")
    ("output_prefix,o", p_opt::value<string>(&result_file)->default_value("log_likelihood"),
            "Output path log.")
    ("temp-dir,d", p_opt::value<string>(&tempDir)->default_value("/tmp"),
            "Temporary (working) directory path")
    ("row,r", p_opt::value<uint64_t>(&max_row)->default_value(150000),
            "Maximum number of rows (classify mode))")
    ("col,c", p_opt::value<uint64_t>(&max_col)->default_value(150000),
            "Maximum number of columns (classify mode)")
  ;
  
  p_opt::positional_options_description pos_args;
  pos_args.add("mode", 1);
  pos_args.add("srcdir", 1);

  p_opt::options_description cmdline_options;
  cmdline_options.add(generic).add(visible).add(hidden);

  p_opt::variables_map opt_map;
  p_opt::store(
    p_opt::command_line_parser(argc, argv).options(cmdline_options)
                                          .positional(pos_args).run(),
    opt_map);
  p_opt::notify(opt_map);

  if(opt_map.count("version")){
    cout<<PROG_VER<<"\n";
    return 1;
  }

  if(opt_map.count("help") || opt_map.count("mode") == 0
     || opt_map.count("srcdir") == 0){
    cout<<usageMsg<<"\n"<<generic<<"\n"<<visible<<"\n";
    return 1;
  }

  srcdir = opt_map["srcdir"].as<string>();
  mode = opt_map["mode"].as<string>();

  create_directories(savedir);
  memLimit *= 1000 * 1000;
  NB nb(kmersize, path(savedir), nthreads);

  Genome::STORE_ALL_NUMERATORS = false;

  if(mode.compare("train") == 0){
    cout<<"Train mode.\n";
    nb.loadTrain();
    trainNB(nb, path(srcdir), kmer_ext, nbatch, memLimit);

    std::ofstream outputFile(savedir + "/" + KMER_NUM_FILE);
    if (outputFile.is_open()) {
        outputFile << nb.getKmerSize();
        outputFile.close();
    } else {
        std::cout << "Failed to open the file for writing." << std::endl;
    }

    cout<<"Training complete.\n";

  }else if(mode.compare("classify") == 0){

    cout<<"Classify mode.\n";
  
    std::ifstream inputFile(savedir + "/" + KMER_NUM_FILE);
    if (inputFile.is_open()) {
        string firstElement;
        if (inputFile >> firstElement) {
            if (firstElement != to_string(nb.getKmerSize())){
                std::cout << "The kmer number used for training (" << firstElement << ") is different from the one used for classification (" << nb.getKmerSize() << ")." << std::endl;
                exit(1);
            }
        } else {
            std::cout << "Failed to extract the kmer number used for training." << std::endl;
            exit(1);
        }
        inputFile.close();
    } else {
        std::cout << "Failed to extract the kmer number used for training." << std::endl;
        exit(1);
    }

    nb.setMaxOutputSize(max_row, max_col);

    nb.loadClassify();
    classifyNB(nb, path(srcdir), kmer_ext, nbatch, memLimit, result_file, full_result, tempDir);

  }else{
    cout<<usageMsg<<"\n"<<generic<<"\n"<<visible<<"\n";
    return 1;
  }

  return 0;
}

