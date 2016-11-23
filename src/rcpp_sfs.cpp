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
sfs__matrix_to_SpMat(SEXP E) {
    Rcpp::NumericMatrix matrix(E);

    // since we don't know the number of nonzeros we first convert everything into triple format:
    std::vector<int> rowidx,colidx;
    std::vector<double> value;
    for(size_t row=0; row<matrix.nrow(); row++) {
        for(size_t col=0; col<matrix.ncol(); col++) {
            double v = matrix(row,col);
            if(v!=0) {
                rowidx.push_back(row);
                colidx.push_back(col);
                value.push_back(v);
            }
        }
    }
    return sfs__tripleform_to_SpMat(rowidx,colidx,value);
}

// from seriation/src/lt.h:
#ifndef LT_POS
#define LT_POS(n, i, j)					\
  (i)==(j) ? 0 : (i)<(j) ? n*((i)-1) - (i)*((i)-1)/2 + (j)-(i) -1	\
        : n*((j)-1) - (j)*((j)-1)/2 + (i)-(j) -1
#endif

static
SFSMatrix::SpMat *
sfs__dist_to_SpMat(SEXP R_dist) {
    // dist objects as used in the seriation package: lower-triangular,
    // without diagonal elements. Typically rather dense; access via
    // LT_POS(size,row_idx,col_idx) (1-based)
    size_t length = LENGTH(R_dist);
    int n = 1 + (int)sqrt(2*length);
    if(length != n*(n-1)/2) {
        throw std::runtime_error("dist object has invalid length");
    }
    const double *dist = REAL(R_dist);

    std::vector<int> rowidx,colidx;
    std::vector<double> value;
    for(size_t row=0; row<n; row++) {
        for(size_t col=0; col<row; col++) {
            double v = dist[LT_POS(length,row+1,col+1)];
            if(v!=0) {
                rowidx.push_back(row);
                colidx.push_back(col);
                value.push_back(v);
                // ... and symmetric copy
                rowidx.push_back(col);
                rowidx.push_back(row);
                value.push_back(v);
            }
        }
    }
    return sfs__tripleform_to_SpMat(rowidx,colidx,value);
}

static
SFSMatrix::SpMat *
sfs__dataframe_to_SpMat(SEXP E) {
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
sfs(SEXP E) {
    SFSMatrix::SpMat *M;
    
    if(Rf_isMatrix(E)) {
        M = sfs__matrix_to_SpMat(E);
    } else if(Rf_inherits(E, "data.frame")) {
        M = sfs__dataframe_to_SpMat(E);
    } else if(Rf_inherits(E, "dist")) {
        M = sfs__dist_to_SpMat(E);
    } else {
        std::cerr << "Unsupported argument type to rcpp_sfs()" << std::endl;
    }

    if(M->n_cols!=M->n_rows) {
        throw std::runtime_error("Input matrix must be square");
    }

    SFSMatrix S(*M);

    arma::Row<int> permutation(S.solve());
    
    delete M;

    return permutation;
}
