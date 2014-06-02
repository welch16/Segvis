//' Find the positions in reads vectors where the regions start
//'
//' @param regionStart - start(regions) from GRanges
//' @param regionEnd - end(regions) from GRanges
//' @param readStart - start(reads) from GRanges
//' @param readEnd - end(reads) from GRanges
//' @param strand - as.character(strand(reads)) form GRanges, we are assuming that they are all the same
//' @param fraglen - Numeric value, fragment length
//' @return position in reads vector

/* Both the reads and regions GRanges object must be ordered
   respect the start positions of ranges */

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector find_anchors(IntegerVector regionStart,
                 IntegerVector regionEnd,
                 IntegerVector readStart,
                 IntegerVector readEnd,
                 CharacterVector strand,
                 int fraglen){
  
  /* The exceptions are in the outside R calling */

  int n = regionStart.size(), m = readStart.size(); 
  int k =0,i=0;
  int extReadStart, extReadEnd;
  Rcpp::IntegerVector startIdx = rep(-1,n);
  bool cond;

  for(k=0;k<n;k++){
    cond = true;
    i=0;
    while((i<m) & cond){ 
      if(strand[i] == "+"){
        extReadStart = readStart[i];
        extReadEnd =  readStart[i] + fraglen -1;
      }else{
        extReadStart = readEnd[i] - fraglen + 1;
        extReadEnd = readEnd[i];
      }
      if(extReadEnd > regionStart[k] && extReadEnd > regionEnd[k]){
        cond = false;
	startIdx[k] = i;
      }
      i++;
    }
  }   
  return startIdx;
}
