//
//  NB.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#include "NB.hpp"

bool NB::OUTPUT_FULL_LOG_LIKELIHOOD = true;
NB::Debug NB::debug_flag = NB::Debug::NO_LOG;

NB::NB(int _kmer_size, path _save_dir, int _nthreads)  {
  nthreads = _nthreads;
  debug_flag = Debug::NO_LOG;
  kmer_size = _kmer_size;
  save_dir = _save_dir;

  load_start_index = 0;
}

NB::~NB() {
  for (unordered_map<string, Class<int>* >::iterator iter = classes.begin();
    iter != classes.end(); iter++)
      delete iter->second;
}

void NB::setMaxOutputSize(uint64_t row, uint64_t col) {
  output_max_row = row;
  output_max_col = col;
}

void NB::loadTrain() {
  vector<path> res = Diskutil::getItemsInDir(save_dir);
  for (vector<path>::iterator iter = res.begin(); iter != res.end(); iter++) {

    if (!Diskutil::hasFileExtension(*iter, Diskutil::SAVE_FILE_EXT)) {
      continue;
    }

    string filename = iter->filename().native();
    string cls_s = filename.substr(0,filename.rfind('-'));
    addClass(new Class<int>(cls_s, kmer_size, *iter));
  }
}

void NB::loadClassify() {
  nthreads -= 1;
  
  vector<path> res = Diskutil::getItemsInDir(save_dir);
  for (vector<path>::iterator iter = res.begin(); iter != res.end(); iter++) {

    if (!Diskutil::hasFileExtension(*iter, Diskutil::SAVE_FILE_EXT)) {
      continue;
    }
    training_genomes.push_back(make_pair(*iter,0));
  }

  if (training_genomes.size() == 0) {
    cout << "Error: No training data found." << endl;
    exit(1);
  }

  std::sort(training_genomes.begin(), training_genomes.end(),
          [](const std::pair<path, uint64_t>& a,
              const std::pair<path, uint64_t>& b) {
              return a.first < b.first;
          });

  // calculate the estimated size of each training data
  for (auto& pair : training_genomes) {
    auto path = pair.first;
    uint64_t file_size = Diskutil::getFileSize(path);

    size_t num_kmer = file_size / (Class<int>::getMapElementSize());
    pair.second = Class<int>::getEstimatedClassBytes(num_kmer);
  }

}

size_t NB::extractHeader(string& inputFile) {
  const string input_fasta_file = inputFile;
  const string output_header_prefix = temp_dir + "/" + temp_prefix + "_0_";
  const string output_header_extension = "_f.hd";
  const string output_max_extension = ".max";
  size_t total_headers = 0;

  int inFileDescriptor = open(input_fasta_file.c_str(), O_RDONLY);
  if (inFileDescriptor == -1) {
      std::cerr << "Error opening file: " << input_fasta_file << endl;
      return 0;
  }

  // Get the file size
  size_t fileSize = lseek(inFileDescriptor, 0, SEEK_END);
  lseek(inFileDescriptor, 0, SEEK_SET);

  // Map the input file into memory
  char* fileMemory = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, inFileDescriptor, 0));
  if (fileMemory == MAP_FAILED) {
      std::cerr << "Error mapping file into memory." << endl;
      close(inFileDescriptor);
      return 0;
  }

  // Extract headers
  std::vector<string> headers;
  char* lineStart = fileMemory;
  for (size_t offset = 0; offset < fileSize; ++offset) {
      if (fileMemory[offset] == '\n') {
          // Found a newline character
          string line(lineStart, fileMemory + offset - lineStart);
          if (!line.empty() && line[0] == '>') {
              // Found a header line
              headers.push_back(line.substr(1)); // Remove the '>' character
              total_headers++;
          }
          lineStart = fileMemory + offset + 1; // Move the lineStart pointer
      }
  }

  // Write headers to the output binary files
  size_t currentHeadersWritten = 0;
  size_t fileIndex = 1;
  std::ofstream headerFile;
  std::ofstream maxFile;
  size_t remainingHeaders = headers.size();
  
  for (const string& header : headers) {
      if (currentHeadersWritten == 0) {
          // Need to open a new output file
          if (headerFile.is_open()) {
              headerFile.close();
          }

          string currentMaxOutputFile = temp_dir + "/" + temp_prefix + "_" + std::to_string(fileIndex) + output_max_extension;
          maxFile.open(currentMaxOutputFile, std::ios::binary);

          string currentHeaderOutputFile = output_header_prefix + std::to_string(fileIndex) + output_header_extension;
          headerFile.open(currentHeaderOutputFile, std::ios::binary);
          if (!headerFile.is_open() || !maxFile.is_open()) {
              std::cerr << "Error opening output file: " << currentHeaderOutputFile << " and " << currentMaxOutputFile << endl;
              munmap(fileMemory, fileSize);
              close(inFileDescriptor);
              return 0;
          }
          fileIndex += output_max_row;
      }

      size_t headerSize = header.size();
      headerFile.write(reinterpret_cast<const char*>(&headerSize), sizeof(headerSize));
      headerFile.write(header.c_str(), headerSize);

      ++currentHeadersWritten;
      --remainingHeaders;

      // Check if the maximum limit for headers per file is reached
      if (currentHeadersWritten >= output_max_row || remainingHeaders == 0) {
          std::vector<double> matrix_data(currentHeadersWritten * 2, 2.0);
          maxFile.write(reinterpret_cast<const char*>(matrix_data.data()), matrix_data.size() * sizeof(double));
          maxFile.close();

          headerFile.close();
          currentHeadersWritten = 0;
      }
  }

  // Clean up
  munmap(fileMemory, fileSize);
  close(inFileDescriptor);

  return total_headers;
}

void NB::setTempDir(const string& _temp_dir) {
  temp_dir = _temp_dir;
}

void NB::save() {
  for (unordered_map<string, Class<int>* >::iterator iter = classes.begin();
      iter != classes.end(); iter++) {
        if (iter->second->loaded()) {
          iter->second->save();
        }
      }
}

void NB::addClass(Class<int>* cl) {
    classes[cl->getId()] = cl;
}

Class<int>* NB::getClass(string cl_id) {
  if (classes.find(cl_id) == classes.end()) {
    return NULL;
  }
  return classes[cl_id];
}

int NB::getKmerSize() {
  return kmer_size;
}

int NB::getThreadNumber() {
  return nthreads;
}

path NB::getSavedir() {
  return save_dir;
}

void NB::addClassToUpdateQueue(Class<int>* cl) {
  classesToProcess.push(cl);
}

void NB::setWriteBufferSize(size_t size) {
  write_buffer_height = size;
}

void NB::initializeOutputBuffer() {
  if (NB::OUTPUT_FULL_LOG_LIKELIHOOD) {
    outputs.resize(write_buffer_height, vector<double>(cls_size + 2, 2));
  }else{
    outputs.resize(write_buffer_height, vector<double>(2, 2));
  }
  
}

void NB::classifyThreadController() {
  while(true) {

    classifyJob job;

    {
      unique_lock<std::mutex> job_lock(classify_job_lock);
      if (classifyJobs.empty()) {
        jobUpdateStatus.wait(job_lock, [this] { return !classifyJobs.empty() || job_done;});
      }

      if (job_done && classifyJobs.empty()) {
        return;
      }

      job = classifyJobs.front();
      classifyJobs.pop();
    }

    unordered_map<int, int> *kmer = new unordered_map<int, int>();
    Diskutil::countKmer(*kmer, kmer_size, *(get<0>(job)), get<2>(job), get<3>(job));

    {
      lock_guard<mutex> lock(num_seq_kmer_counted_access);
      
      num_seq_kmer_processed += 1;
      if (waiting_for_kmer_counting && num_seq_kmer_processed == num_seq_kmer_processed) {
        num_seq_kmer_counted_cv.notify_one();
      }
    }

    if (kmer->size() == 0) {
      NB::outputs[(get<1>(job) - start_seq_index) % write_buffer_height][0] = -1;
    } else {

      Genome genome(".", ".");
      pair<int, double> predicted_class = make_pair(0,0);
      genome.setKmerCounts(kmer);

      for (uint64_t i = 0; i < cls_size; i++) {
        string filename = training_genomes[load_start_index-cls_size+i].first.filename().native();
        string cls_s = filename.substr(0,filename.rfind('-'));
        Class<int>* cl = classes[cls_s];

        double result = genome.computeClassificationNumerator(cl);
        if (NB::OUTPUT_FULL_LOG_LIKELIHOOD) {
          NB::outputs[(get<1>(job) - start_seq_index) % write_buffer_height][i+2] = result;
        }

        if (result > predicted_class.second || predicted_class.second == 0) {
          predicted_class.first = load_start_index-cls_size+i;
          predicted_class.second = result;
        }
      }

      NB::outputs[(get<1>(job) - start_seq_index) % write_buffer_height][0] = predicted_class.first;
      NB::outputs[(get<1>(job) - start_seq_index) % write_buffer_height][1] = predicted_class.second;

      genome.resetKmerCounts();
    }

    delete kmer;
    

    lock_guard<mutex> output_update_lock(output_mtx);

    processed_seq_num++;

    if (processed_seq_num % write_buffer_height == 0) {
      output_cv.notify_one();
    }
  }
}

void NB::startClassifyThreads() {
  initializeOutputBuffer();
  job_done = false;
  for (uint64_t i=0; i<nthreads; i++) {
    threads.push_back(thread(&NB::classifyThreadController, this));
  }
}

void NB::queueClassifyJob(vector<char>& buffer, tuple<uint64_t, uint64_t, uint64_t>& sequence_data) {

  if (start_seq_index == 0) {
    start_seq_index = get<2>(sequence_data);
  }

  num_seq_kmer_processing += 1;

  // write to CSV if exceeds buffer size
  if ((get<2>(sequence_data) - start_seq_index) % write_buffer_height == 0 && get<2>(sequence_data) != start_seq_index) {

    unique_lock<mutex> output_lock(output_mtx);
    
    if (processed_seq_num % write_buffer_height != 0) {
      output_cv.wait(output_lock, [this]{ return processed_seq_num % write_buffer_height == 0;});
    }

    if (!finished_writing) {
      unique_lock<mutex> writing_lock(output_modify);
      write_done_cv.wait(writing_lock, [this]{ return finished_writing;});
    }
    
    swap(outputs, outputs_write);

    lock_guard<mutex> start_write_lock(start_write_mtx);
    finished_writing = false;
    start_write_cv.notify_one();

    initializeOutputBuffer();

  }
  
  classifyJob job(&buffer, get<2>(sequence_data), get<0>(sequence_data), get<1>(sequence_data));

  lock_guard<mutex> job_lock(classify_job_lock);

  classifyJobs.push(job);
  jobUpdateStatus.notify_one();
}

void NB::getClassMemoryUsage(uint64_t& class_avg_bytes, uint64_t& class_max_bytes) {
  for (auto it = training_genomes.begin(); it != training_genomes.end(); it++) {
    uint64_t class_bytes = it->second;
    class_avg_bytes += class_bytes;
    if (class_bytes > class_max_bytes) {
      class_max_bytes = class_bytes;
    }
  }
  class_avg_bytes /= training_genomes.size();
}

void NB::getStringSeqHeaders(vector<string>& seq_headers, const string& file_name) {
  // Open the binary file for reading
  std::ifstream file(file_name, std::ios::binary);
  if (!file.is_open()) {
      std::cerr << "Error opening file: " << file_name << endl;
      return;
  }

  while (!file.eof()) {
      // Read the size of the next string
      size_t size;
      if (!file.read(reinterpret_cast<char*>(&size), sizeof(size))) {
          if (file.eof()) {
              break;  // Reached the end of the file
          } else {
              std::cerr << "Error reading size from file." << endl;
              return;
          }
      }

      // Read the string data
      std::vector<char> buffer(size);
      if (!file.read(buffer.data(), size)) {
          std::cerr << "Error reading string data from file." << endl;
          return;
      }

      // Append the extracted string to the vector
      seq_headers.emplace_back(buffer.begin(), buffer.end());
  }

  // Close the file
  file.close();
}

void NB::getStringClassHeaders(vector<string>& cls_headers, vector<pair<void*, size_t>>& maps, vector<size_t>& offsets, vector<size_t>& cols) {
  for (size_t i = 0; i < maps.size(); i++) {
    size_t total_headers = *reinterpret_cast<size_t*>(maps[i].first) - 2 ;
    cols[i] = total_headers;
    char* data_ptr = reinterpret_cast<char*>(maps[i].first) + sizeof(size_t); // Move past the total_headers value

    for (size_t i = 0; i < total_headers; ++i) {
        size_t header_length = *reinterpret_cast<size_t*>(data_ptr);

        data_ptr += sizeof(size_t); // Move past the header_length value

        string header(data_ptr, header_length);
        data_ptr += header_length; // Move past the actual header data

        cls_headers.push_back(header);
    }

    offsets[i] = data_ptr - reinterpret_cast<char*>(maps[i].first);
  }
}

void NB::concatenateCSVByColumns(const std::vector<string>& inputFiles, const string& output_prefix, const std::size_t& sequence_num) {

  vector<string> seq_headers;

  getStringSeqHeaders(seq_headers, inputFiles[0]);

  // load all the data into mmap
  vector<pair<void*, size_t>> fileMaps(inputFiles.size() - 1);
  for (size_t i = 0; i < fileMaps.size(); i++) {
    string filename = inputFiles[i+1];

    // Open the file
    int fileDescriptor = open(filename.c_str(), O_RDONLY);
    if (fileDescriptor == -1) {
        std::cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Get the file size
    struct stat fileInfo;
    if (fstat(fileDescriptor, &fileInfo) == -1) {
        std::cerr << "Error getting file size: " << filename << endl;
        close(fileDescriptor);
        return;
    }

    // Map the file into memory
    void* fileMap = mmap(nullptr, fileInfo.st_size, PROT_READ, MAP_PRIVATE, fileDescriptor, 0);
    if (fileMap == MAP_FAILED) {
        std::cerr << "Error mapping file to memory: " << filename << endl;
        close(fileDescriptor);
        return;
    }
    fileMaps[i] = pair<void*, size_t>(fileMap, fileInfo.st_size);

    // Close the file descriptor
    close(fileDescriptor);
  }

  // extract the class headers
  vector<string> cls_headers;
  vector<size_t> maps_start(inputFiles.size()-1,0);
  vector<size_t> file_col_num(inputFiles.size()-1,0);
  getStringClassHeaders(cls_headers, fileMaps, maps_start, file_col_num);

  int num_csv = cls_headers.size() / output_max_col;
  
  bool is_frag = cls_headers.size() % output_max_col != 0;
  if (is_frag) num_csv++;

  size_t start_file_index = 0;
  size_t skip_cols = 0;

  // pair<start_index, num_cols>
  vector<pair<size_t, size_t>> write_start_pos(fileMaps.size(), pair<size_t, size_t>(0,0));
  for (size_t j = 0; j < file_col_num.size(); j++) {
    write_start_pos[j] = make_pair(maps_start[j] + (2 * sizeof(double)), file_col_num[j]);
  }

  for (double i = 0; i < num_csv; i++) {
    size_t class_num = last_written_class_index + 2 + i*output_max_col;
    string output_name = output_prefix + "_" + std::to_string(class_num) + "_" + to_string(sequence_num) + ".csv";
    size_t new_start_file_index = start_file_index;

    size_t num_files = 0;
    size_t cols_cnt = 0;
    if (skip_cols != 0) {
      write_start_pos[start_file_index].first += (skip_cols * sizeof(double));
      size_t offset = (write_start_pos[start_file_index].first - maps_start[start_file_index]) / sizeof(double) - 2;
      write_start_pos[start_file_index].second = file_col_num[start_file_index] - offset;
      skip_cols = 0;
    }

    for (size_t j = start_file_index; j < fileMaps.size(); j++) {
      if (write_start_pos[j].second + cols_cnt > output_max_col) {
        num_files++;
        
        size_t col_to_write = output_max_col - cols_cnt;
        write_start_pos[j].second = col_to_write;
        skip_cols = col_to_write;

        cols_cnt = output_max_col;

        break;
      } 
      else{

        cols_cnt += write_start_pos[j].second;
        new_start_file_index++;

        num_files++;
      }
    }

    bool is_binary = false;
    // the last csv file whose size is smaller than output_max_col
    if (i +  1 >= num_csv && is_frag && load_start_index < training_genomes.size()) {
      output_name = temp_dir + "/" + temp_prefix + "_concat_" + to_string(sequence_num) + ".tmp";

      std::ifstream file(output_name);
      if (file.good()) {
          file.close();
          if (std::remove(output_name.c_str()) != 0) {
              cout << "Error: Failed to delete the file: " << output_name << endl;
          }
      }

      is_binary = true;
    }

    std::ofstream outputFile(output_name);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << endl;
        return;
    }

    if (is_binary) {
      size_t index = fileMaps.size() - 1;
      void* file_to_write = fileMaps[index].first;

      size_t cls_cnt = 0;
      for (size_t i = 0; i < file_col_num.size(); i++) {
        cls_cnt += file_col_num[i];
      }
      size_t cls_num = (cls_cnt % output_max_col) + 2;

      outputFile.write(reinterpret_cast<const char*>(&cls_num), sizeof(cls_num));

      for (size_t j = file_col_num[index] + 2 - cls_num; j < file_col_num[index]; j++) {
        size_t header_size = cls_headers[j].size();
        outputFile.write(reinterpret_cast<const char*>(&header_size), sizeof(header_size));
        outputFile.write(cls_headers[j].c_str(), header_size);
      }

      outputFile.write(reinterpret_cast<char*>(file_to_write + maps_start[index]), 2 * sizeof(double));
  
      for (size_t row = 0; row < seq_headers.size(); row++) {
        size_t skipping_bytes = (2 + file_col_num[index]) * sizeof(double) * row;
        skipping_bytes += write_start_pos[index].first;
        outputFile.write(reinterpret_cast<char*>(file_to_write + skipping_bytes), (2 *  sizeof(double)));
        outputFile.write(reinterpret_cast<char*>(file_to_write + skipping_bytes), (write_start_pos[index].second) * sizeof(double));
      }

      outputFile.close();

      for (size_t i = 0; i < fileMaps.size(); i++) {
        if (munmap(fileMaps[i].first, fileMaps[i].second) == -1) {
            std::cerr << "Error unmapping memory for file at index " << i << endl;
        }
      }

      return;
    }

    outputFile << "sequence_header,";
    for (size_t j = i * output_max_col; j < cls_headers.size() && j < (i+1) * output_max_col; j++) {
      outputFile << cls_headers[j] << ",";
    }
    outputFile << "\n";
    for (size_t row = 0; row < seq_headers.size(); row++) {
      outputFile << seq_headers[row] << ",";
      size_t total_cols_written = 0;

      size_t iter = start_file_index + num_files;
      for (size_t current_map = start_file_index; current_map < iter; current_map++) {

        size_t skipping_bytes = (2 + file_col_num[current_map]) * sizeof(double) * row;
        double* data_ptr = reinterpret_cast<double*>(fileMaps[current_map].first + write_start_pos[current_map].first + skipping_bytes);
        size_t cols_to_write = write_start_pos[current_map].second;

        for (size_t col = 0; col < cols_to_write; col++) {

          if (*data_ptr > 0) {
            if (current_map == start_file_index) {
              outputFile << "invalid sequence" << ",\n";
            }
            
            data_ptr+=cols_to_write;
            total_cols_written = 0;
            break;
          }

          outputFile << *data_ptr << ",";
          data_ptr ++;
          total_cols_written++;
        }

        if (total_cols_written == cols_cnt) {
          outputFile << "\n";
          total_cols_written = 0;
        }
      }
    }

    outputFile.close();
    start_file_index = new_start_file_index;
  }


  for (size_t i = 0; i < fileMaps.size(); i++) {
    if (munmap(fileMaps[i].first, fileMaps[i].second) == -1) {
        std::cerr << "Error unmapping memory for file at index " << i << endl;
    }
  }

  for (size_t i = 1; i < inputFiles.size(); i++) {
    string concat_filename = temp_prefix + "_concat_" + to_string(sequence_num) + ".tmp";
    if (inputFiles[i].find(concat_filename) != string::npos) {

      if (remove(inputFiles[i].c_str()) != 0) {
        std::cerr << "Error: Failed to delete temporary file: " << inputFiles[i] << endl;
        return;
      }
      
      break;
    }
  }
  
}

void NB::waitCalculatingAllKmers() {
  if (num_seq_kmer_processing != 0) {
    unique_lock<mutex> lock(num_seq_kmer_counted_access);
    waiting_for_kmer_counting = true;
    if (num_seq_kmer_processing != num_seq_kmer_processed) {
      num_seq_kmer_counted_cv.wait(lock, [this] { return num_seq_kmer_processing == num_seq_kmer_processed; });
    }
    
    waiting_for_kmer_counting = false;
  }
}

string NB::getAppendFiles(int filter_number, std::vector<string>& file_names, const string& directory) {
    string class_header_file = "";
    DIR* dir = opendir(directory.c_str());
    if (dir) {
        struct dirent* ent;
        while ((ent = readdir(dir)) != nullptr) {
            string file_name = ent->d_name;

            if (class_header_file == "") {
                string pattern = "^" + temp_prefix + "_" + std::to_string(output_class_index) + R"(\.clshd)";
                std::regex regexPattern(pattern);

                if (std::regex_match(file_name, regexPattern)) {
                    class_header_file = temp_dir + "/" + file_name;
                }

            } 
            
            string target_name = temp_prefix + "_" + std::to_string(filter_number) + "_";
            if (file_name != "." && file_name != ".." && file_name.find(".tmp") != string::npos && file_name.find(target_name) == 0) {
                file_name = temp_dir + "/" + file_name;
                file_names.push_back(file_name);
            }
        }
        closedir(dir);
          
      std::sort(file_names.begin(), file_names.end(), [this](const string& a, const string& b) {
          string aStr = a.substr(0, a.rfind(".tmp"));
          string bStr = b.substr(0, b.rfind(".tmp"));

          int aNum = std::stoi(aStr.substr(aStr.find_last_of('_') + 1));
          int bNum = std::stoi(bStr.substr(bStr.find_last_of('_') + 1));

          return aNum < bNum;
      });

    } else {
        std::cerr << "Error opening directory." << endl;
    }
  
  return class_header_file;
}

void NB::loadMaxResults(vector<pair<int,string>>& max_files, vector<vector<double>>& max_list) {
      for (const auto& element : max_files) {
        string file_path = element.second;
        // Open the file for reading
        int fd = open(file_path.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "Error opening file: " << file_path << endl;
            continue;
        }

        // Get the file size
        off_t file_size = lseek(fd, 0, SEEK_END);
        lseek(fd, 0, SEEK_SET);

        // Map the file into memory
        void* file_memory = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (file_memory == MAP_FAILED) {
            std::cerr << "Error mapping file: " << file_path << endl;
            close(fd);
            continue;
        }

        size_t max_size = file_size / sizeof(double);
        double* max_ptr = reinterpret_cast<double*>(file_memory);
        
        max_list.emplace_back(max_size); // Add an empty vector<double>
        std::vector<double>& max = max_list.back(); // Get a reference to the added vector

        // Populate the new vector
        for (size_t i = 0; i < max_size; ++i) {
          max[i] = max_ptr[i];
        }

        munmap(file_memory, file_size); // Unmap the file

        // Close the file (the mapped memory remains valid)
        close(fd);
    }
}


void NB::compareAndWriteMax(vector<double> max_vector, vector<double>& compare_vector, const string& filepath) {
    
    // Open the output file for writing
    std::ofstream outFile(filepath, std::ios::binary);
    if (!outFile) {
        std::cerr << "Failed to open output file: " << filepath << endl;
        return;
    }

    // Iterate through odd-indexed values, compare, and modify compare_vector
    for (size_t i = 1; i < max_vector.size(); i += 2) {
        double index = max_vector[i-1];
        double value = max_vector[i];

        if (value > compare_vector[i] && value <= 0) {
            // If max_vector's odd-indexed value is greater, modify compare_vector
            compare_vector[i] = value;
            compare_vector[i - 1] = index;
        }
    }

    // Write the modified compare_vector to the output file
    outFile.write(reinterpret_cast<const char*>(compare_vector.data()), sizeof(double) * max_vector.size());

    // Close the output file
    outFile.close();
}

void NB::fullAppend(std::vector<string>& inputFiles, const string& outputFile, string class_header_file) {

  size_t col = 2;
  uint64_t bytesWritten = 0;
  int outputFileIndex = 1;

  size_t header_size;
  vector<char> buffer;
  
  std::ifstream classHeaderFile(class_header_file, std::ios::binary | std::ios::ate);

  // Get the file size
  header_size = classHeaderFile.tellg();
  classHeaderFile.seekg(0, std::ios::beg);

  // Read the entire file content into a vector
  buffer.resize(header_size);
  if (!classHeaderFile.read(buffer.data(), header_size)) {
      std::cerr << "Error reading file." << endl;
      return;
  }

  // Process the file content from the buffer
  const char* ptr = buffer.data();
  const char* endPtr = ptr + header_size;

  while (ptr < endPtr) {
      size_t size;
      std::memcpy(&size, ptr, sizeof(size));
      ptr += sizeof(size);
      if (ptr + size <= endPtr) {
          col++;
          ptr += size;
      } else {
          // Incomplete string data, handle the error as needed
          break;
      }
  }

  classHeaderFile.close();

  vector<pair<int,string>> max_files;
  getMaxFiles(max_files);
  vector<vector<double>> max_list;
  loadMaxResults(max_files, max_list);
  size_t max_result_index = 0;
  vector<double> new_max_result(output_max_row * 2, 0);
  size_t new_result_size = 0;

  uint64_t MAX_FILE_SIZE = (output_max_row * col) * sizeof(double);
  string currentOutputFile = outputFile + "_" + std::to_string(outputFileIndex) + "_f.tmp";

  // Append the input files to the output file(s) with a byte limit
  std::ofstream outFile(currentOutputFile, std::ios::binary | std::ios::app);
  if (outFile.is_open()) {
      // Write the number of strings as uint64_t
      outFile.write(reinterpret_cast<const char*>(&col), sizeof(col));
      outFile.write(buffer.data(), header_size);

      // Use mmap to load each input file and append its contents
      for (const string& inputFile : inputFiles) {
          int fd = open(inputFile.c_str(), O_RDONLY);
          if (fd != -1) {
              struct stat stat_buf;
              fstat(fd, &stat_buf);
              size_t fileSize = static_cast<size_t>(stat_buf.st_size);

              void* fileMemory = mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
              size_t num_row = fileSize / (col * sizeof(double));
              
              double* doubleArray = static_cast<double*>(fileMemory);

              for (size_t i = 0; i < num_row; i++) {
                double  first_column = doubleArray[col * i];
                double second_column = doubleArray[col * i + 1];

                new_max_result[2 * new_result_size] = first_column;
                new_max_result[2 * new_result_size + 1] = second_column;
                new_result_size++;

                if (new_result_size == output_max_row) {
                  compareAndWriteMax(max_list[max_result_index], new_max_result, max_files[max_result_index].second);
                  new_result_size = 0;
                  max_result_index++;
                }
              }

              size_t totalBytesWritten = 0;
              size_t file_bytes_left = fileSize;
              while (file_bytes_left != 0) {
                  // Determine how many bytes can be written to the current output file
                  size_t bytesToWrite = std::min(file_bytes_left, MAX_FILE_SIZE - bytesWritten);

                  if (bytesToWrite % col != 0 && bytesToWrite < fileSize) {
                    // Adjust bytesToWrite to the nearest multiple of col
                    bytesToWrite = (bytesToWrite / col) * col;
                  }

                  // Write the data to the current output file
                  outFile.write(reinterpret_cast<const char*>(fileMemory) + totalBytesWritten, bytesToWrite);

                  // Update the total bytes written and check if a new output file is needed
                  bytesWritten += bytesToWrite;
                  file_bytes_left -= bytesToWrite;
                  totalBytesWritten += bytesToWrite;

                  if (bytesWritten >= MAX_FILE_SIZE && file_bytes_left != 0) {
                      outFile.close();

                      bytesWritten = 0;
                      outputFileIndex += output_max_row;
                      currentOutputFile = outputFile + "_" + std::to_string(outputFileIndex) + "_f.tmp";

                      outFile.open(currentOutputFile, std::ios::binary | std::ios::app);

                      if (!outFile.is_open()) {
                          std::cerr << "Failed to open a new output file." << endl;
                          return;
                      }

                      outFile.write(reinterpret_cast<const char*>(&col), sizeof(col));
                      outFile.write(buffer.data(), header_size);
                  }
              }
              
              munmap(fileMemory, fileSize);

              close(fd);
              if (remove(inputFile.c_str()) != 0) {
                  std::cerr << "Error: Failed to delete temporary file: " << inputFile << endl;
              }

          } else {
              std::cerr << "Failed to open input file: " << inputFile << endl;
              return;
          }
      }

      if (max_result_index < max_files.size()) {
        compareAndWriteMax(max_list[max_result_index], new_max_result, max_files[max_result_index].second);
      }

      outFile.close();
  }

  if (remove(class_header_file.c_str()) != 0) {
    std::cerr << "Error: Failed to delete temporary file: " << class_header_file << endl;
  }
}

void NB::maxAppend(std::vector<string>& inputFiles) {
    size_t col = 2;

    vector<pair<int,string>> max_files;
    getMaxFiles(max_files);
    vector<vector<double>> max_list;
    loadMaxResults(max_files, max_list);
    size_t max_result_index = 0;
    vector<double> new_max_result(output_max_row * 2, 0);
    size_t new_result_size = 0;

    uint64_t MAX_FILE_SIZE = (output_max_row * col) * sizeof(double);

    // Use mmap to load each input file and append its contents
    for (const string& inputFile : inputFiles) {
        int fd = open(inputFile.c_str(), O_RDONLY);
        if (fd != -1) {
            struct stat stat_buf;
            fstat(fd, &stat_buf);
            size_t fileSize = static_cast<size_t>(stat_buf.st_size);

            void* fileMemory = mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
            size_t num_row = fileSize / (col * sizeof(double));
            
            double* doubleArray = static_cast<double*>(fileMemory);

            for (size_t i = 0; i < num_row; i++) {
              double  first_column = doubleArray[col * i];
              double second_column = doubleArray[col * i + 1];

              new_max_result[2 * new_result_size] = first_column;
              new_max_result[2 * new_result_size + 1] = second_column;

              new_result_size++;

              if (new_result_size == output_max_row) {
                compareAndWriteMax(max_list[max_result_index], new_max_result, max_files[max_result_index].second);
                new_result_size = 0;
                max_result_index++;
              }
            }
            
            munmap(fileMemory, fileSize);

            close(fd);
            if (remove(inputFile.c_str()) != 0) {
                std::cerr << "Error: Failed to delete temporary file: " << inputFile << endl;
            }

        } else {
            std::cerr << "Failed to open input file: " << inputFile << endl;
            return;
        }
    }

    if (max_result_index < max_files.size()) {
      compareAndWriteMax(max_list[max_result_index], new_max_result, max_files[max_result_index].second);
    }
}

void NB::getConcatFiles(std::unordered_map<size_t, std::vector<string>>& file_groups, const string& directory) {
    DIR* dir = opendir(directory.c_str());
    if (dir) {
        struct dirent* ent;
        while ((ent = readdir(dir)) != nullptr) {
            string file_name = ent->d_name;

            std::smatch match;

            int group_number = -1;
            string pattern = "^" + temp_prefix + R"(_concat_(\d+)\.tmp)";
            if (std::regex_search(file_name, match, std::regex(pattern))) {
                group_number = std::stoi(match[1]);
            } else {
              pattern = "^" + temp_prefix + R"(_\d+_(\d+)_f\.(tmp|hd))";
              if (std::regex_search(file_name, match, std::regex(pattern))) {
                  group_number = std::stoi(match[1]);
              }
            }
            
            if (group_number >= 0) {
              file_groups[group_number].push_back(file_name);
            }
        }
        closedir(dir);

        // Sort files within each group
        for (auto& kv : file_groups) {
          std::sort(kv.second.begin(), kv.second.end(), [](const string& file1, const string& file2) {
          float num1, num2 = 0;
          std::smatch match1, match2;

          string pattern = R"(_concat_(\d+)\.tmp)";
          if (std::regex_search(file1, match1, std::regex(pattern))) {
              num1 = 0.5;
          } else {
            pattern = R"(_(\d+)_\d+_f\.(tmp|hd))";
            if (std::regex_search(file1, match1, std::regex(pattern))) {
                num1 = std::stoi(match1[1]);
            }
          }

          pattern = R"(_concat_(\d+)\.tmp)";
          if (std::regex_search(file2, match2, std::regex(pattern))) {
              num2 = 0.5;
          } else {
            pattern = R"(_(\d+)_\d+_f\.(tmp|hd))";
            if (std::regex_search(file2, match2, std::regex(pattern))) {
                num2 = std::stoi(match2[1]);
            }
          }

          return num1 < num2;
        });
      }

      for (auto& kv : file_groups) {
        for (auto& file : kv.second) {
          file = directory + "/" + file;
        }
      }

    } else {
        std::cerr << "Error opening directory." << endl;
    }
}

void NB::concatenateCSVs(string& output_name) {

  unordered_map<size_t, vector<string>> file_groups;
  getConcatFiles(file_groups, temp_dir);

  for (auto it = file_groups.begin(); it != file_groups.end(); it++) {
    string output_prefix = output_name;
    concatenateCSVByColumns(it->second, output_prefix, it->first);
  }

  // Loop through each group of files
  for (const auto& group : file_groups) {
      const std::vector<string>& filenames = group.second;

      for (int i = 1; i < filenames.size(); ++i) {
          const string& filename = filenames[i];
          
          string concat_filename = temp_prefix + "_concat_" + to_string(group.first) + ".tmp";
          if (filename.find(concat_filename) != string::npos) {
              continue;
          }

          if (std::remove(filename.c_str()) != 0) {
              std::cerr << "Error deleting file: " << filename << endl;
          }
      }

  }
  
}

void NB::joinClassifyThreads() {

  job_done = true;

  {
    lock_guard<std::mutex> job_lock(classify_job_lock);
    jobUpdateStatus.notify_all();
  }

  for (uint64_t i=0; i<threads.size(); i++) {
    threads[i].join();
  }
  threads.clear();

}

void NB::processClassUpdates() {
  progress = 1;
  total = classesToProcess.size();

  if (NB::debug_flag == NB::Debug::LOG_SOME) {
    cout<<"Started training on batch...\n";
    cout.flush();
  }

  int c_nthreads = nthreads;
  if (nthreads > classes.size()) {
    c_nthreads = classes.size();
  }

  vector<thread*> workers(c_nthreads);
  for (int i=0; i < c_nthreads; i++) {
    workers[i] = new thread(&NB::trainThreadController, this);
  }
  for (int i=0; i < c_nthreads; i++) {
    workers[i]->join();
    delete workers[i];
  }
}

void NB::trainThreadController() {
  bool exitCondition = false;
  while(!exitCondition) {
    classQueueAccess.lock();
    if (classesToProcess.empty()) {
      exitCondition = true;
      classQueueAccess.unlock();
      continue;
    }

    Class<int>* cl = classesToProcess.front();
    classesToProcess.pop();

    if (NB::debug_flag == NB::Debug::LOG_SOME
       || NB::debug_flag == NB::Debug::LOG_ALL) {
      cout<<"("<<progress++<<"/"<<total<<") Training genomes for class ";
      cout<<cl->getId()<<"\n";
      cout.flush();
    }
    classQueueAccess.unlock();

    if (cl->loaded()) {
      cout<<"Attempted double-load for Class " << cl->getId()<<"\n";
    } else {
      cl->load();
      cl->addGenomesInQueue();
      cl->save();
      cl->unload();
    }
    
    std::lock_guard<std::mutex> lock(process_update_access);
    classes.erase(cl->getId());

    delete cl;

  }
}

void NB::writeToCSV() {
  
  while(true) {

    {
      unique_lock<mutex> start_write_lock(start_write_mtx);
      start_write_cv.wait(start_write_lock, [this]{return outputs_write.size() > 0 || (job_done && hasLoadedAll() && outputs.size() == 0);});
    }

    if (job_done && hasLoadedAll() && outputs_write.size() == 0) {
      return;
    }

    if (output_seq_index == 1 && NB::OUTPUT_FULL_LOG_LIKELIHOOD) {
      string header_name = temp_dir + "/" + temp_prefix + "_" + to_string(output_class_index) + ".clshd";

      std::ofstream header_file(header_name);
      if (!header_file) {
          cout << "Error Creating Output File." << "  " << header_name << endl;
          return;
      }

      for (size_t i = 0; i < cls_size; ++i) {
        string filename = training_genomes[load_start_index-cls_size+i].first.filename().native();
        string cls_s = filename.substr(0,filename.rfind('-'));
        
        size_t size = cls_s.size();
        header_file.write(reinterpret_cast<const char*>(&size), sizeof(size));

        header_file.write(cls_s.c_str(), size);
      }

      header_file.close();
    }

    string filename = temp_dir + "/" + temp_prefix + "_" + to_string(output_class_index) + "_" + to_string(output_seq_index) + ".tmp";
    std::ofstream outputFile(filename);
    if (!outputFile) {
        cout << "Error Creating Output File." << "  " << filename << endl;
        return;
    }

    uint64_t cols = outputs_write[0].size();
    for (const auto& row : outputs_write) {
        if (row[0] > 0 && row[1] > 0) {
          break;
        }

        outputFile.write(reinterpret_cast<const char*>(row.data()), cols * sizeof(double));
    }

    outputFile.close();
    output_seq_index ++;

    lock_guard<mutex> writing_lock(output_modify);
    finished_writing = true;
    outputs_write.clear();
    write_done_cv.notify_one();
  }
}

bool NB::hasLoadedAll() {
  return load_start_index >= training_genomes.size();
}

void NB::unloadClassesThreadProcess() {
  while(true) {
    unique_lock<mutex> lock(memoryCntAccess);
    if (classesToProcess.empty()) {
      lock.unlock();
      return;
    }else{
      Class <int> *cl = classesToProcess.front();
      classesToProcess.pop();
      lock.unlock();

      string cls_s = cl->getId();
      cl->unload();
      
      std::lock_guard<std::mutex> lock(process_update_access);
      classes.erase(cls_s);

      delete cl;

      //cout << "removed: " << cls_s << endl;
    }
    
  }
}

void NB::unloadClasses() {

  for (auto cls : classes) {
    classesToProcess.push(cls.second);
  }

  for (size_t i = 0; i < nthreads; ++i)
  {
    threads.push_back(thread(&NB::unloadClassesThreadProcess, this));
  }

  for (auto& thread : threads)
  {
    thread.join();
  }

  threads.clear();

  wrapUpBatch();

  output_class_index += cls_size;

  if (load_start_index >= training_genomes.size()) {
    vector<pair<int,string>> seqHeaderFiles;
    getSeqHeaderFiles(seqHeaderFiles);

    vector<pair<int,string>> maxFiles;
    getMaxFiles(maxFiles);

    outputMaxResults(seqHeaderFiles, maxFiles);

    vector<string> extensions = {".hd", ".max"};
    removeFilesWithExtensions(extensions);

    lock_guard<mutex> start_write_lock(start_write_mtx);
    outputs.clear();
    outputs_write.clear();
    start_write_cv.notify_one();
  }

  cls_size = 0;

  cout << "Current Batch Finished unloading" << endl;
}

void NB::wrapUpBatch() {

  if (!finished_writing) {
    unique_lock<mutex> writing_lock(output_modify);
    write_done_cv.wait(writing_lock, [this]{ return finished_writing;});
  }

  bool finished = outputs[0][0] == 2 && outputs[0][1] == 2;

  if (!finished) {
    swap(outputs, outputs_write);

    {
      lock_guard<mutex> start_write_lock(start_write_mtx);
      finished_writing = false;
      start_write_cv.notify_one();
    }

    if (!finished_writing) {
      unique_lock<mutex> writing_lock(output_modify);
      write_done_cv.wait(writing_lock, [this]{ return finished_writing;});
    }
  }

  std::vector<string> files;

  string class_header_file = getAppendFiles(output_class_index, files, temp_dir);
  string append_output_name = temp_dir + "/" + temp_prefix + "_" + std::to_string(output_class_index);
 
  if (NB::OUTPUT_FULL_LOG_LIKELIHOOD) {
    fullAppend(files, append_output_name, class_header_file);

    if (load_start_index - last_written_class_index > output_max_col || load_start_index >= training_genomes.size()) {
      concatenateCSVs(output_prefix);
      last_written_class_index = (load_start_index / output_max_col) * output_max_col - 1;
    }
  }else{
    maxAppend(files);
  }

  
  start_seq_index = 0;
  processed_seq_num = 0;
  num_seq_kmer_processed = 0;
  num_seq_kmer_processing = 0;

}

void NB::loadClassesThreadProcess(uint64_t &max_memory, uint64_t &used_memory) {

  while(true) {
    unique_lock<mutex> lock(memoryCntAccess);

    if (used_memory + training_genomes[load_start_index].second > max_memory
        || load_start_index >= training_genomes.size()) {

      lock.unlock();
      return;

    }else{
      uint64_t load_index = load_start_index;
      used_memory += training_genomes[load_start_index].second;    

      load_start_index += 1;
      path file = training_genomes[load_index].first;
      string filename = file.filename().native();
      string cls_s = filename.substr(0,filename.rfind('-'));

      Class<int> *cl = new Class<int>(cls_s, kmer_size, file);

      addClass(cl);
      
      lock.unlock();

      cl->load();

      std::lock_guard<std::mutex> lock(process_update_access);
      //cout << "loaded: " << cl->getId() << endl;
    }
    
  }

}

void NB::loadClasses(uint64_t &max_memory) {
    
  uint64_t memory_used = 0;

  for (int i = 0; i < nthreads; i++) {
    threads.push_back(thread(&NB::loadClassesThreadProcess, this, ref(max_memory), ref(memory_used)));
  }

  for (int i = 0; i < nthreads; i++) {
    threads[i].join();
  }
  threads.clear();

  if (classes.size() == 0) {
    cout << "Error: Failed Loading Classes." << endl;
    exit(1);
  }

  output_seq_index = 1;
  cls_size = classes.size();

  cout << "Current Batch Finished loading. Size: " << cls_size << endl;
}

void NB::getSeqHeaderFiles(vector<pair<int,string>>& seq_header_files) {
  DIR* dir = opendir(temp_dir.c_str());
  if (dir) {
      struct dirent* ent;
      while ((ent = readdir(dir)) != nullptr) {
          string file_name = ent->d_name;

          std::smatch match;

          int seq_num = 0;

          string pattern_str = R"(\b)" + temp_prefix + "_" + R"((\d+)_(\d+)_f\.hd\b)";
          if (std::regex_search(file_name, match, std::regex(pattern_str))) {
            seq_num = std::stoi(match[2]);
          }
          
          if (seq_num > 0) {
            seq_header_files.push_back(make_pair(seq_num, file_name));
          }
      }
      closedir(dir);

    std::sort(seq_header_files.begin(), seq_header_files.end(), [](const std::pair<int, string>& a, const std::pair<int, string>& b) {
      return a.first < b.first;
    });

    for (auto& kv : seq_header_files) {
      kv.second = temp_dir + "/" + kv.second;
    }

  } else {
      std::cerr << "Error opening directory." << endl;
  }
}

void NB::getMaxFiles(vector<pair<int,string>>& max_files) {
  DIR* dir = opendir(temp_dir.c_str());
  if (dir) {
    struct dirent* ent;
    while ((ent = readdir(dir)) != nullptr) {
        string file_name = ent->d_name;

        string pattern_str = R"(\b)" + temp_prefix + "_" + R"((\d+)\.max\b)";
        std::regex pattern(pattern_str);
        std::smatch matches;

        int seq_num = 0;
        if (std::regex_match(file_name, matches, pattern)) {
            seq_num = std::stoi(matches[1]); // Convert the captured group to an integer
        }
        
        if (seq_num > 0) {
          max_files.push_back(make_pair(seq_num, file_name));
        }
    }
    closedir(dir);

    std::sort(max_files.begin(), max_files.end(), [](const std::pair<int, string>& a, const std::pair<int, string>& b) {
      return a.first < b.first;
    });

    for (auto& kv : max_files) {
      kv.second = temp_dir + "/" + kv.second;
    }

  } else {
      std::cerr << "Error opening directory." << endl;
  }
}

void NB::outputMaxResults(vector<pair<int,string>>& seq_header_files, vector<pair<int,string>>& max_files) {
  for (size_t i = 0; i < seq_header_files.size(); i++) {

    vector<string> seq_headers;
    getStringSeqHeaders(seq_headers, seq_header_files[i].second);

    string filename = max_files[i].second;

    int fd = open(filename.c_str(), O_RDWR); // Open the file

    if (fd == -1) {
        std::cerr << "Failed to open file." << endl;
        return;
    }

    // Get the file size
    off_t file_size = lseek(fd, 0, SEEK_END);

    // Map the file to memory
    void* mapped_data = mmap(nullptr, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    string output_name = output_prefix + "_max_" + std::to_string(seq_header_files[i].first) + ".csv";
    std::ofstream output_file(output_name);

    if (!output_file.is_open()) {
      cout << "Unable to open the file for writing." << endl;
    }

    double* data = (double*)mapped_data;
    for (size_t j = 0; j < seq_headers.size(); j++) {
      if (*data < 0) {
        output_file << seq_headers[j] << "," << "sequence contains no valid kmers" << ","<< endl;
        data++;
        data++;
        continue;
      }

      string filename = training_genomes[*data].first.filename().native();
      string cls_s = filename.substr(0,filename.rfind('-'));
      data++;

      output_file << seq_headers[j] << "," << cls_s << "," << *data << endl;
      data++;
    }

    munmap(mapped_data, file_size);

    output_file.close();
  }
}

void NB::setOutputPrefix(const string& _output_prefix) {
  output_prefix = _output_prefix;

  const char pathSeparator = '/';

  size_t pos = output_prefix.find_last_of(pathSeparator);
  
  if (pos != string::npos) {
      temp_prefix = output_prefix.substr(pos + 1);
  } else {
    temp_prefix = output_prefix;
  }

}

void NB::removeFilesWithExtensions(const std::vector<string>& extensions) {
    DIR* dir = opendir(temp_dir.c_str());
    if (dir == nullptr) {
        std::cerr << "Error opening directory." << endl;
        return;
    }

    dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_REG) { // Check if it's a regular file
            string filename = entry->d_name;

            if (filename.find(temp_prefix + "_") != 0) {
              continue;
            }

            for (const string& ext : extensions) {
                if (filename.size() >= ext.size() &&
                    filename.compare(filename.size() - ext.size(), ext.size(), ext) == 0) {
                    string filepath = temp_dir + "/" + filename;
                    if (unlink(filepath.c_str()) != 0) {
                        std::cerr << "Error removing: " << filepath << endl;
                    }
                }
            }
        }
    }

    closedir(dir);
}