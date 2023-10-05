//
//  Class.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#include "Genome.hpp"
#include "NB.hpp"

template <class T>
Class<T>::Class(string id_, int kmer_size, path _savefile)  {
  ngenomes_lg = make_pair(0, false);
  sumfreq_lg = make_pair(0, false);

  ngenomes=0;

  sumfreq=1;
  sumfreq<<=2*kmer_size;

  id = id_;

  savefile = _savefile;
}

template <class T>
Class<T>::~Class() {
  for(vector<Genome*>::iterator iter = genomes.begin(); iter != genomes.end();
      iter++){
        delete *iter;
      }
}

template <class T>
void Class<T>::queueGenome(Genome* genome){
  queuedGenomes.push(genome);
}

template <class T>
void Class<T>::addGenomesInQueue(){
  while(!queuedGenomes.empty()){
    Genome* current = queuedGenomes.front();
    addGenome(current);
    current->unload();
    queuedGenomes.pop();
  }
}

template <class T>
double Class<T>::logAdd(vector<double> exponents){
  double result = 0, max;
  int maxId = -1;
  for(vector<double>::iterator num=exponents.begin();
    num != exponents.end(); num++){
    if(maxId == -1 || max < *num){
      maxId = 0;
      max = *num;
    }
  }

  for(vector<double>::iterator num=exponents.begin();
    num != exponents.end(); num++){
    result += exp(*num - max);
  }
  result = max + log(result);
  
  return result;
}

template <class T>
void Class<T>::addNGenomes(int ngenomes_){
  if(ngenomes_lg.second){
    ngenomes = 0;
    ngenomes_lg.second = false;
  }
  ngenomes += ngenomes_;
}

template <class T>
void Class<T>::addSumfreq(long long int sumfreq_){
  if(sumfreq_lg.second){
    sumfreq = 0;
    sumfreq_lg.second = false;
  }
  sumfreq += sumfreq_;
}

template <class T>
size_t Class<T>::getEstimatedClassBytes(size_t& num_elements){
  size_t key_memory = num_elements * sizeof(T);
  size_t value_memory = num_elements * sizeof(double_wflag);
  size_t bucket_memory = (num_elements * 1.5) * sizeof(std::list<std::pair<T, double_wflag>>);
  size_t approx_memory = key_memory + value_memory + bucket_memory;

  return approx_memory;
}

template <class T>
size_t Class<T>::getMapElementSize(){
  return sizeof(T) + sizeof(double);
}

template <class T>
void Class<T>::addFreqCount(T feature_, int count_){
  if(getFreqcnt().find(feature_) == getFreqcnt().end()){
    getFreqcnt()[feature_] = count_;
  }else{
    getFreqcnt()[feature_] += count_;
  }
  if(getFreqcnt_lg().find(feature_) == getFreqcnt_lg().end()){
      getFreqcnt_lg()[feature_] = make_pair(0 ,false);
  }else{
      getFreqcnt_lg()[feature_].second = false;
  }
}

template <class T>
void Class<T>::addSumfreq_lg(double sumfreq_){
  vector<double> items;
  items.push_back(getSumFreq_lg());
  items.push_back(sumfreq_);
  sumfreq_lg.first = logAdd(items);
}

template <class T>
void Class<T>::addNGenomes_lg(double ngenomes_){
  vector<double> items;
  items.push_back(getNGenomes_lg());
  items.push_back(ngenomes_);
  ngenomes_lg.first = logAdd(items);
}

template <class T>
void Class<T>::addFreqCount_lg(T feature_, double count_){
  vector<double> items;
  items.push_back(getFreqCount_lg(feature_));
  items.push_back(count_);
  getFreqcnt_lg()[feature_].first = logAdd(items);
}

template <class T>
path Class<T>::getSavefilePath(){
  return savefile;
}

template <class T>
double Class<T>::getNGenomes_lg(){
  if(!ngenomes_lg.second){
    vector<double> terms;
    terms.push_back(ngenomes_lg.first);
    terms.push_back(log(ngenomes));
    ngenomes = 0;
    ngenomes_lg.first = logAdd(terms);
    ngenomes_lg.second = true;
  }
  return ngenomes_lg.first;
}

template <class T>
double Class<T>::getFreqCount_lg(T feature_){
  if(!existsInFreqCount(feature_)){
    return 0;
  }else{
    if(!getFreqcnt_lg()[feature_].second){
      vector<double> terms;
      terms.push_back(getFreqcnt_lg()[feature_].first);
      terms.push_back(log(getFreqcnt()[feature_]));
      getFreqcnt()[feature_] = 0;
      getFreqcnt_lg()[feature_].first = logAdd(terms);
      getFreqcnt_lg()[feature_].second = true;
    }
    return getFreqcnt_lg()[feature_].first;
  }
}

template <class T>
double Class<T>::getSumFreq_lg(){
  if(!sumfreq_lg.second){
    vector<double> terms;
    terms.push_back(sumfreq_lg.first);
    terms.push_back(log(sumfreq));
    sumfreq = 0;
    sumfreq_lg.first = logAdd(terms);
    sumfreq_lg.second = true;
  }
  return sumfreq_lg.first;
}

template <class T>
bool Class<T>::existsInFreqCount(T feature_){
  return (getFreqcnt().find(feature_) != getFreqcnt().end())
    || (getFreqcnt_lg().find(feature_) != getFreqcnt_lg().end());
}

template <class T>
string Class<T>::getId(){
  return id;
}

template <class T>
uint64_t Class<T>::serializeDouble(double value){
  uint64_t result;
  assert(sizeof(result) == sizeof(value));
  memcpy(&result, &value, sizeof(result));
  return result;
}

template <class T>
double Class<T>::deserializeDouble(uint64_t value){
  double result;
  assert(sizeof(result) == sizeof(value));
  memcpy(&result, &value, sizeof(result));
  return result;
}

template <class T>
void Class<T>::addGenome(Genome* genome, double confidence_lg){

  if(confidence_lg == 0.0)
    addNGenomes(1);
  else
    addNGenomes_lg(confidence_lg + log(1));

  int sum=0, new_features=0;
  unordered_map<int, int> freqs = genome->getKmerCounts();
  for(unordered_map<int, int>::iterator iter = freqs.begin();
    iter != freqs.end(); iter++){
      sum += iter->second;
      if(existsInFreqCount(iter->first)){
        new_features++;
      }
      if(confidence_lg == 0.0)
        addFreqCount(iter->first, iter->second);
      else
        addFreqCount_lg(iter->first, confidence_lg + log(iter->second));
    }

  if(confidence_lg == 0.0)
    addSumfreq(sum);
  else
    addSumfreq_lg(confidence_lg + log(sum));
    
    
  genomes.push_back(genome);
}

template <class T>
void Class<T>::computeBatchNumerators(vector<Genome*> genomes){
  for(vector<Genome*>::iterator iter = genomes.begin(); iter != genomes.end();
    iter++){
      (*iter)->computeClassificationNumerator(this);
  }
}

template <class T>
void Class<T>::load(path source_path){
  if(isLoaded){
    return;
  }else{
    isLoaded = true;
  }

  freqcnt_lg = new unordered_map<T, double_wflag>;
  freqcnt = new unordered_map<T, int>;
  
  std::ifstream in(source_path.native());
  if(in.good()){
    deserialize(in);
  }

  in.close();
}

template <class T>
void Class<T>::save(path destination_path){
  std::ofstream out(destination_path.native());
  // out<<serialize();
  serialize(out);
  out.close();
}

template <class T>
void Class<T>::save(){
  save(savefile);
}

template <class T>
void Class<T>::load(){
  load(savefile);
}

template <class T>
void Class<T>::unload(){

  freqcnt->clear();
  std::unordered_map<T, int>().swap(*freqcnt);

  freqcnt_lg->clear();
  std::unordered_map<T, std::pair<double, bool>>().swap(*freqcnt_lg);

  delete freqcnt;
  delete freqcnt_lg;

  isLoaded = false;
}

template <class T>
bool Class<T>::loaded(){
  return isLoaded;
}

template <class T>
unordered_map<T, pair<double, bool> >& Class<T>::getFreqcnt_lg(){
  return *freqcnt_lg;
}

template <class T>
unordered_map<T, int>& Class<T>::getFreqcnt(){
  return *freqcnt;
}

template <class T>
size_t Class<T>::getClassSizeInBytes(){
  return (*freqcnt_lg).size() * (sizeof(T) + sizeof(double_wflag)) + sizeof(Class<T>);
}

/*
 *  WARNING: Serialization may be platform-dependent, as C++11 does not
 *  guarantee all 'double' types to follow the IEEE-754 64-bit standard.
 */
template <class T>
void Class<T>::serialize(std::ofstream& out) {

    // Serialize ngenomes_lg
    double ngenomes_lg = getNGenomes_lg();
    out.write(reinterpret_cast<const char*>(&ngenomes_lg), sizeof(ngenomes_lg));

    // Serialize sumfreq_lg
    double sumfreq_lg = getSumFreq_lg();
    out.write(reinterpret_cast<const char*>(&sumfreq_lg), sizeof(sumfreq_lg));

    // Serialize freqcnt_lg.size()
    int freqcntSize = getFreqcnt_lg().size();
    out.write(reinterpret_cast<const char*>(&freqcntSize), sizeof(freqcntSize));

    // Serialize freqcnt_lg
    for (const auto& iter : getFreqcnt_lg()) {
        T key = iter.first;
        double value = getFreqCount_lg(iter.first);
        
        // Serialize key
        out.write(reinterpret_cast<const char*>(&key), sizeof(key));
        
        // Serialize value
        out.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }
}

template <class T>
void Class<T>::deserialize(std::ifstream& in) {
    long long int tmp;

    // Read and convert ngenomes_lg
    in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    ngenomes_lg = make_pair(reinterpret_cast<double&>(tmp), true);

    // Read and convert sumfreq_lg
    in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    sumfreq_lg = make_pair(reinterpret_cast<double&>(tmp), true);

    int n, kmer;

    // Read the number of elements in the map
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    
    getFreqcnt_lg().reserve(n);

    for (int i = 0; i < n; i++) {
      // Read the kmer value
      in.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));

      // Read and convert the tmp value
      in.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));

      getFreqcnt_lg()[kmer] = make_pair(reinterpret_cast<double&>(tmp), true);
    }

}
