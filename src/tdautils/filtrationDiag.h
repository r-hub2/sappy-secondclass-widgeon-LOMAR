#ifndef __FILTRATIONDIAG_H__
#define __FILTRATIONDIAG_H__

#include <vector>
#include <string>
#include <algorithm>

#include <tdautils/filtrationUtils.h>

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>


// for grid
#include <tdautils/gridUtils.h>

#include <iostream>



// FiltrationDiag
/** \brief Interface for R code, construct the persistence diagram from the
*         filtration.
*
* @param[out] Rcpp::List     A list
* @param[in]  filtration     The input filtration
* @param[in]  maxdimension   Max dimension of the homological features to be
*                            computed.
* @param[in]  library        Either "GUDHI", "Dionysus", or "PHAT"
* @param[in]  location       Are location of birth point, death point, and
*                            representative cycles returned?
* @param[in]  printProgress  Is progress printed?
*/
// TODO: see whether IntegerVector in template is deducible
template< typename VertexVector, typename VectorList, typename RealVector >
inline void filtrationDiagSorted(
    VectorList        & cmplx,
    RealVector        & values,
    const int           maxdimension,
    const std::string & library,
    const bool          location,
    const bool          printProgress,
    const unsigned      idxShift,
    std::vector< std::vector< std::vector< double > > > & persDgm,
    std::vector< std::vector< std::vector< unsigned > > > & persLoc,
    std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {

  if (library[0] == 'D') {
    FiltrationDiagDionysus< Persistence >(
        filtrationTdaToDionysus< VertexVector, Fltr >(
            cmplx, values, idxShift),
        maxdimension, location, printProgress, persDgm, persLoc, persCycle);
  }
 
}



// FiltrationDiag
/** \brief Interface for R code, construct the persistence diagram from the
*         filtration.
*
* @param[out] Rcpp::List     A list
* @param[in]  filtration     The input filtration
* @param[in]  maxdimension   Max dimension of the homological features to be
*                            computed.
* @param[in]  library        Either "GUDHI", "Dionysus", or "PHAT"
* @param[in]  location       Are location of birth point, death point, and
*                            representative cycles returned?
* @param[in]  printProgress  Is progress printed?
*/
// TODO: see whether IntegerVector in template is deducible
template< typename VertexVector, typename VectorList, typename RealVector >
inline void filtrationDiag(
    VectorList        & cmplx,
    RealVector        & values,
    const int           maxdimension,
    const std::string & library,
    const bool          location,
    const bool          printProgress,
    const unsigned      idxShift,
    std::vector< std::vector< std::vector< double > > > & persDgm,
    std::vector< std::vector< std::vector< unsigned > > > & persLoc,
    std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {

  if (std::is_sorted(values.begin(), values.end())) {
    filtrationDiagSorted< VertexVector >(
        cmplx, values, maxdimension, library, location, printProgress,
        idxShift, persDgm, persLoc, persCycle);
  }
  else {
    std::vector< std::vector< unsigned > > cmplxTemp = 
        RcppCmplxToStl< std::vector< unsigned >, VertexVector >(cmplx, 0);
    std::vector< double > valuesTemp(values.begin(), values.end());
    filtrationSort(cmplxTemp, valuesTemp);
    filtrationDiagSorted< std::vector< unsigned > >(
        cmplxTemp, valuesTemp, maxdimension, library, location, printProgress,
        idxShift, persDgm, persLoc, persCycle);
  }
}



# endif // __FILTRATIONDIAG_H__
