//counting kmers in R


// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

#include <iostream>
//#include <utility>
#include <cmath>
#include <vector>
#include <map>
// #include <stdlib.h>

using namespace Rcpp;
using namespace std;


//encode binary
inline const short int encode_bp (const char x, unsigned int& currN_val){
  switch(x){
  case 'N':case 'n': //illegal value
    return (currN_val++ % 4);
  case 'A':case 'a':
    return 0;
    break;
  case 'C':case 'c':
    return 1;
    break;
  case 'G':case 'g':
    return 2;
    break;
  default:
    return 3;
  break;
  }
}

inline unsigned int  encode_oligo(const std::string oligo, const int k, const int km1, unsigned int& currN_val){
  unsigned int  result=0;
  for (int i=0; i<k; i++) {
    result+= pow(4, km1-i) * encode_bp(oligo[i],currN_val);
  }
  return result;
}




//decode binary
inline const char decode_bp (const int x){
  switch(x){
  case 0:
    return 'A';
    break;
  case 1:
    return 'C';
    break;
  case 2:
    return 'G';
    break;
  default:
    return 'T';
  break;
  }
}
string decode_oligo(const unsigned int arrSlotNum, const int k, const int km1,char * kmerName){
  // char * kmerName= (char *)calloc(k * sizeof(char),0);    // char * kmerName= new char[k]; //wierd additional char when k>7
  int last2BitMask= 3;
  for (int i=0; i<k; i++) {
    kmerName[i]= decode_bp(arrSlotNum >> (2*(km1-i)) & last2BitMask);    // std::cout <<arrSlotNum<< "  >>" << (2*(km1-i))<<" "<< (arrSlotNum >> (2*(km1-i)))<< " <> " <<(arrSlotNum >> (2*(km1-i)) & last2BitMask)<<" "<< kmerName[i]<< " "<<kmerName<< "\n";
  }
  return string(kmerName);;
}







// N is assigned with ACGT recycled
// collapse only work when all lig the same length
// [[Rcpp::export]]
SEXP kmerCntBit(std::vector<std::string> strings, int k=2, bool diffLen=false, bool collapse=false, bool asDf= true, bool all_possible_k=true, int pseudo=0)
{
  if (pseudo) all_possible_k=true; //return all kmers if using pseudo

  int kmerPerLig= strings[0].length() -k + 1;  //how many kmers per lig
  size_t ligNo = strings.size(); //total lig No.
  unsigned long kmerMask = 0;
      for (int i = 0; i < k; i++)
      {
        kmerMask <<= 2;
        kmerMask |= 3;
      }
  int km1=k-1; //for quick use
  unsigned int currN_val=0;


  std::vector<int> kcnt(pow(4,k), pseudo); //vect for result, ini with pseudo

  std::vector< std::vector<int> > posCnt; //store pos infomation when collapse=F
  if (!collapse) {
    if (!all_possible_k)
      {posCnt.resize(pow(4,k));}
    else
      {posCnt.resize(pow(4,k), std::vector<int>(kmerPerLig,pseudo));} };

  unsigned long curr_k_code=0;
  for (int i = 0; i < ligNo; i++)
  {
    if (diffLen) kmerPerLig= strings[i].length() -km1;
    for (int j=0; j<kmerPerLig; j++)
    {
      if (j == 0)
        {
          curr_k_code= encode_oligo(strings[i].substr(j, k), k, km1, currN_val);
        } else
        {
          int curr_char_code = encode_bp(strings[i][j+km1], currN_val);
          curr_k_code = (curr_k_code << 2 | curr_char_code) & kmerMask;
        }

      if(!collapse)
        {
          if( (!all_possible_k) && posCnt[curr_k_code].empty()) {posCnt[curr_k_code].resize(kmerPerLig,pseudo);} //ini if not ini ed
          posCnt[curr_k_code][j]++;
        }else
        {
          kcnt[curr_k_code]++;  //#std::string sub_s= strings[i].substr(j, k); kcnt[encode_oligo(sub_s,k)]++;
        }
    }
  }



  //------------------
  size_t vectLen=pow(4, k);
  char * kmerName= (char *)calloc(k * sizeof(char),0); //for use in decoding kmerName
  if(all_possible_k)
    {
      vector< string >allNames(pow(4,k));
      for (size_t i=0; i<vectLen; i++) { allNames[i]= decode_oligo(i, k, km1,kmerName); }

      if(!collapse){
        List result= wrap(posCnt);
        result.attr("names")= allNames;
          // wrapping to df is still a bit time consuming
          IntegerVector row_name(kmerPerLig);
          for (int i=0; i< kmerPerLig; i++) {row_name[i]= i+1;} //std::to_string(i+1)
          result.attr("class") = "data.frame";
          result.attr("row.names") = row_name;
        return result;
      }else{
        if (asDf) {return DataFrame::create(Named("kmer")= allNames, Named("counts")= wrap(kcnt), _["stringsAsFactors"] = false);}
        else {IntegerVector result= wrap(kcnt); result.attr("names")= allNames; return result; }
      }
    }
        //else
        if(!collapse)
          {
            std::vector< std::vector<int> >nonEmpty;
            std::vector< std::string >nonEmpty_name;
            for (size_t i=0; i<vectLen; i++)
            {
              if (!posCnt[i].empty()) {nonEmpty.push_back(posCnt[i]); nonEmpty_name.push_back(decode_oligo(i, k, km1,kmerName)); }
            }
            List result= wrap(nonEmpty);
            result.attr("names")= nonEmpty_name;

          // wrapping to df is still a bit time consuming
            IntegerVector row_name(kmerPerLig);
            for (int i=0; i< kmerPerLig; i++) {row_name[i]= i+1;} //std::to_string(i+1)
            result.attr("class") = "data.frame";
            result.attr("row.names") = row_name;

            return result;
          }
        else  //if collapse, then consider asDF
          {
            std::vector< int >nonEmpty;
            std::vector< std::string >nonEmpty_name;
            for (size_t i=0; i<vectLen; i++)
              {
              if (kcnt[i]!=0) {nonEmpty.push_back(kcnt[i]); nonEmpty_name.push_back(decode_oligo(i, k, km1,kmerName));}
              }
            if (asDf) {return DataFrame::create(Named("kmer")= nonEmpty_name, Named("counts")= nonEmpty, _["stringsAsFactors"] = false);}
            else {IntegerVector result= wrap(nonEmpty); result.attr("names")= nonEmpty_name; return result; }
          }

}


/*** R

*/
