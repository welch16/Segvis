//' Creates a list with the matched reads for each regions
//'
//' @param regionStart - start(regions) from GRanges
//' @param regionEnd - end(regions) from GRanges
//' @param readStart - start(reads) from GRanges
//' @param readEnd - end(reads) from GRanges
//' @param strand - as.character(strand(reads)) form GRanges
//' @param fraglen - Numeric value, fragment length
//' @return position in reads vector

#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
List match_reads(IntegerVector regionStart,
                 IntegerVector regionEnd,
                 IntegerVector readStart,
                 IntegerVector readEnd,
                 CharacterVector strand,
                 int fraglen){

  int n = regionStart.size(),m = readStart.size();
  int extReadStart, extReadEnd,nextReadStart,nextReadEnd;
  int j=0;
  Rcpp::List match(n);
  bool cond;
  for(int i =0; i<n;i++){
    Rcpp::IntegerVector reads(0);
    cond = true;
    j =0;
    if(strand[j] == "+"){
      extReadStart = readStart[j];
      extReadEnd = readStart[j] + fraglen -1;
    }else{
      extReadStart = readEnd[j] - fraglen +1;
      extReadEnd = readEnd[j];
    }    
    while(j < m && cond){
      if(strand[j] == "+"){
        nextReadStart = readStart[j+1];
        nextReadEnd = readStart[j+1] + fraglen -1;
      }else{
        nextReadStart = readEnd[j+1] - fraglen +1;
        nextReadEnd = readEnd[j+1];
      }
      if( regionStart[i] <= extReadEnd && regionEnd[i] >= extReadStart){
        reads.push_back(j+1);        
        if( regionStart[i] > nextReadEnd || regionEnd[i] < nextReadStart)cond = false;
      }
      extReadStart = nextReadStart;
      extReadEnd = nextReadEnd;
      j++;
    }
    match(i) =  reads;    
  }
  return match;
}
