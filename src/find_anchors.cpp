#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector find_anchors(IntegerVector regionStart,
                 IntegerVector regionEnd,
                 IntegerVector readStart,
                 IntegerVector readEnd,
                 CharacterVector strand,
                 int fraglen){
  /* The peaks and the reads are sorted */

  /* Need to add exception to check that all the peak and read vectors
     have the same length respectively */
  int n = regionStart.size(), m = readStart.size();
  int k =0,i=0;
  int extReadStart, extReadEnd;
  //  Rcpp::List match(n);
  Rcpp::IntegerVector startIdx(n);
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
      if(extReadStart > regionStart[k] && extReadEnd < regionEnd[k]){
        cond = false;
        startIdx[k] = i;
      }
      i++;
    }
  }   
  return startIdx;
}
