/* -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

   rcpp_sfs.cpp -- R wrappers for SFS class
   
   This file is part of SFS.
   
   SFS is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   Foobar is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
   
   (C) 2016 Utz-Uwe Haus, Cray EMEA Research Lab.
*/


#include "RcppArmadillo.h"

#include "SFSMatrix.h"

#include <iostream>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

template <typename T>
static std::ostream&
operator<< (std::ostream& os, const std::vector<T>& v) 
{
  bool first=true;
  os << "[";
  for (typename std::vector<T>::const_iterator ii = v.begin();
       ii != v.end();
       ++ii,first=false)
    os << (first? "" : " ") << *ii;
  os << "]";
  return os;
}

static
SFSMatrix::SpMat *
sfs__tripleform_to_SpMat(const std::vector<int>&    rowidx,
                         const std::vector<int>&    colidx,
                         const std::vector<double>& val) {
    
    // prepare batch insertion format for armadillo
    arma::umat locations(2,rowidx.size());
    arma::vec  values(rowidx.size());
    for(size_t k=0; k<rowidx.size(); k++) {
        locations(0,k) = rowidx[k];
        locations(1,k) = colidx[k];
        values(k) = val[k];
    }

    return new SFSMatrix::SpMat(locations,values,true);
}

static
SFSMatrix::SpMat *
sfs__matrix_to_SpMat(SEXP E, const double zero_eps) {
    Rcpp::NumericMatrix matrix(E);

    // since we don't know the number of nonzeros we first convert everything into triple format:
    std::vector<int> rowidx,colidx;
    std::vector<double> value;
    for(size_t row=0; row<matrix.nrow(); row++) {
        for(size_t col=0; col<matrix.ncol(); col++) {
            double v = matrix(row,col);
            if(fabs(v)>zero_eps) {
                rowidx.push_back(row);
                colidx.push_back(col);
                value.push_back(v);
            }
        }
    }
    return sfs__tripleform_to_SpMat(rowidx,colidx,value);
}

// from seriation/src/lt.h: col/row 1-based
#ifndef LT_POS
#define LT_POS(n, col, row)					\
    (col)==(row) ? 0                                            \
          : (col)<(row)                                         \
                  ? n*((col)-1) - (col)*((col)-1)/2 + (row)-(col) -1    \
                  : n*((row)-1) - (row)*((row)-1)/2 + (col)-(row) -1
#endif

static
SFSMatrix::SpMat *
sfs__dist_to_SpMat(SEXP R_dist, const double zero_eps) {
    // dist objects as used in the seriation package: lower-triangular,
    // without diagonal elements. Typically rather dense; access via
    // LT_POS(size,col_idx,row_idx) (1-based)
    size_t length = LENGTH(R_dist);
    int n = 1 + (int)sqrt(2*length);
    if(length != n*(n-1)/2) {
        throw std::runtime_error("dist object has invalid length");
    }
    const double *dist = REAL(R_dist);

    std::vector<int> rowidx,colidx;
    std::vector<double> value;
    for(size_t col=0; col<n; col++) {
        for(size_t row=0; row<col; row++) {
            double v = dist[LT_POS(n,col+1,row+1)];
            if(fabs(v)>zero_eps) {
                rowidx.push_back(row);
                colidx.push_back(col);
                value.push_back(v);
                // ... and symmetric copy
                rowidx.push_back(col);
                colidx.push_back(row);
                value.push_back(v);
            }
        }
    }

    return sfs__tripleform_to_SpMat(rowidx,colidx,value);
}

static
SFSMatrix::SpMat *
sfs__dataframe_to_SpMat(SEXP E, const double zero_eps) {
    Rcpp::DataFrame D = Rcpp::as<Rcpp::DataFrame>(E);
    // if it has 3 columns assume they're (row, col, val) format
    // ... and fill an SpMat.
    size_t num_cols = D.size();
    if(num_cols!=3) {
        throw std::runtime_error("sfs algorithm can only deal with 3-column dataframes");
    }

    size_t num_rows = D.nrows();
    Rcpp::IntegerVector rowidx = D[0];
    Rcpp::IntegerVector colidx = D[1];
    Rcpp::NumericVector value = D[2];

    // FIXME: kill entries below ZERO_EPS threshold
    
    arma::umat locations(2,num_rows);
    arma::vec  values(num_rows);
    for(size_t k=0; k<num_rows; k++) {
        locations(0,k) = rowidx[k];
        locations(1,k) = colidx[k];
        values(k) = value[k];
    }

    return new SFSMatrix::SpMat(locations,values,true);
}


// exported functionality: Run SFS on various input data formats.
// We support 'matrix', the dense numeric format
//            'dataframe', which is assumed to have 3 columns (row-index, col-index, value)
//                         (column names are ignored)
//            'Matrix', sparse matrices from the Matrix package
//            'dist', lower-triangular distance matrices as used in the seriation package

// [[Rcpp::export]]
arma::Row<int>
sfs(SEXP dissimilarity) {
    SFSMatrix::SpMat *M;
    // when to consider an entry zero. Should be a parameter to sfs call.
    const double zero_eps=1.0e-200; // this should catch most noise for now
    
    if(Rf_isMatrix(dissimilarity)) {
        M = sfs__matrix_to_SpMat(dissimilarity,zero_eps);
    } else if(Rf_inherits(dissimilarity, "data.frame")) {
        M = sfs__dataframe_to_SpMat(dissimilarity,zero_eps);
    } else {
        // some object
        Rcpp::RObject o(dissimilarity);
        if(o.inherits("dist")) {
            M = sfs__dist_to_SpMat(dissimilarity,zero_eps);
        } else {
            std::string c = o.hasAttribute("class")
                    ? "(class " + Rcpp::as<std::string>(o.attr("class")) +")"
                    : "(no class)";
            std::string errmsg = "sfs(): unknown object " + c + " in dissimilarity argument";
            throw std::runtime_error(errmsg);
        }
    }

    if(M->n_cols!=M->n_rows) {
        throw std::runtime_error("Input matrix must be square");
    }

    SFSMatrix S(*M);

    arma::Row<int> permutation(S.solve());
    
    delete M;

    return permutation;
}
