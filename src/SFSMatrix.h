/* -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

   SFSMatrix.h
   
   This file is part of SFS.
   
   SFS is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SFS is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
   
   (C) 2014-2016 Matteo Seminaroti, CWI Amsterdam.
   R integration (C) 2016 Utz-Uwe Haus, Cray EMEA Research Lab.
*/

#ifndef SFSMATRIX_H
#define SFSMATRIX_H

#include <vector>
#include <ostream>
#include <istream>
#include <string>
#include <fstream>
#include <limits>
#include <assert.h>
#include <random>
#include <list>
#include <armadillo>
#include <time.h>

class SFSMatrix
{
  public:
    //general types
    typedef double data_type;
    typedef std::size_t size_t;
    typedef arma::uword uword;
    typedef std::string string;
    typedef std::vector<std::vector<data_type> > adjacency_matrix;
    typedef std::vector<int> IntVector;
    typedef std::vector<data_type> vector;
    
    //Block
    typedef std::list<uword> Block; //used to define the classes of the similarity partition
    typedef Block::iterator iterator; // iterator in the block (used to pass through a block)
    typedef std::vector<iterator> map_iterator; // vector of iterator ot the position of each vertex in the blocks in the queue
    
    //WeakLinearOrder
    typedef std::list<Block> WeakLinearOrder; // represents the similarity partittion (N.B. each block is a vector of integer)
    typedef WeakLinearOrder::iterator block_iterator; //iterator among the blocks (used to pass through a weak liner order)
    
    //Queue
    typedef std::pair<Block, uword> Queue_Block; //each pair is a class of the queue of unvisited vertices
    typedef std::list<Queue_Block> Queue; //the bool value is true if the block is involved in the partition refinement (used for the queue of unvisited vertices)
    typedef Queue::iterator queue_iterator; //iterator among the blocks of the queue of unvisited vertices Q
    typedef std::vector<queue_iterator> map_queue_iterator; //map from the vertices to the corresponding block of the queue of unvisited vertices
    typedef std::list<queue_iterator> Queue_visited_list; //list of blocks of the queue Q whose elements are in the neighborhood of the current pivot
    
    //Neighborhood (for sorting the neighborhood)
    typedef std::pair<uword, data_type> neighbor; //neighbhor of a given pivot p (vertex + similarity)
    typedef std::list<neighbor> neighborhood; //list of neighbors

    struct sortByHighest {
        bool operator()(const neighbor& lhs, const neighbor& rhs) {
            return lhs.second > rhs.second;
        }
    };
    
    //sparse matrix Armadillo
    typedef arma::SpMat<data_type> SpMat; //sparse matrix
    typedef arma::mat Mat;
    typedef arma::vec vec;
    typedef arma::uvec uvec;
    
    //Vertices properties
    struct VertexProperty
    {
        uword nominal;
        bool visited;
        bool neighbor; //used to check if a vertex is adjacent to a vertex in the connected component
    };
    typedef std::vector<VertexProperty> permutation; //our map of vertices
    
    
    //class functions
  public:
    // SFSMatrix(SpMat& A)
    // : _A()
    // , _n()
    // , _m()
    // , _tau_inv()
    // , _binary()
    // , _file_name()
    // , _coco()
    // , _epsilon()
    // {
    //   assert(A.n_rows==A.n_cols);
    //   std::cerr << A << std::endl;
    //     _A = A;
    //     std::cerr << "copied input" <<std::endl;
    //     _n = _A.n_rows;
    //     _m = _A.n_nonzero;
    //     _epsilon = 0;
    //     _binary = binary();
    //     _tau_inv.resize(_n);
    //     std::cerr << "resized" <<std::endl;
    //     VertexProperty vp;
        
    //     for (int i = 0; i < _n; ++i)
    //     {
    //         vp.nominal = i;
    //         vp.visited = false;
    //         vp.neighbor = false;
    //         _tau_inv[i] = vp;
    //     }
    //     std::cerr << "vp filled" <<std::endl;

    // }

  public:
    SFSMatrix::SpMat read(const std::string& file) {
        SFSMatrix::IntVector row, column;
        SFSMatrix::vector value;
        std::ifstream myfile;
        myfile.open(file);
        SFSMatrix::string line;
        if (myfile.is_open())
        {
            std::getline(myfile,line);
            int size_line = std::count(line.begin(), line.end(), ' ');
            if (size_line <= 3) //sparse format
            {
                while (!line.empty())
                {
                    std::stringstream ss(line);
                    int r = -1;
                    int c = -1;
                    data_type v = -1;
                    ss >> r;
                    row.push_back(r-1);
                    ss >> c;
                    column.push_back(c-1);
                    ss >> v;
                    value.push_back(v-1);
                    if (row.back() != column.back()) //add symmetric edges
                    {
                        row.push_back(c-1);
                        column.push_back(r-1);
                        value.push_back(v-1);
                    }
                    std::getline(myfile,line);
                }
            }
            else //matrix format
            {
                int n = size_line;
                for(int i =0; i < n; ++i)
                {
                    std::stringstream ss(line);
                    int tmp = -1;
                    for (int j = 0; j < n; ++j)
                    {
                        ss >> tmp;
                        if (tmp > 0)
                        {
                            row.push_back(i);
                            column.push_back(j);
                            value.push_back(tmp);
                        }
                    }
                    std::getline(myfile,line);
                }
            }
        }
        else {
            throw std::runtime_error("The file could not be opened");
        }
    
        //eliminate identical values
        int i = 0;
        int j = 0;
    
        while (i < row.size() -1)
        {
            j = i + 1;
            while (j < row.size())
            {
                if (row[i] == row[j] && column[i] == column[j])
                {
                    row.erase(row.begin()+j);
                    column.erase(column.begin()+j);
                }
                else
                {
                    j++;
                }
            }
            i++;
        }
    
        //create matrix
        arma::umat locations(2,row.size());
        SFSMatrix::vec values(row.size());
        for (int k=0; k < row.size(); ++k)
        {
            locations(0,k) = row[k];
            locations(1,k) = column[k];
            values(k) = 1;//value[k];
        }
    
        return SFSMatrix::SpMat(locations, values, true);
    }
  public:
    SFSMatrix(const SpMat& A) {
        _A = A;

        assert(_A.n_rows==_A.n_cols);

        _n = _A.n_rows;
        _m = _A.n_nonzero;
        _epsilon = 0;
        _binary = binary();
        _tau_inv.resize(_n);
        VertexProperty vp;
        
        for (int i = 0; i < _n; ++i)
        {
            vp.nominal = i;
            vp.visited = false;
            vp.neighbor = false;
            _tau_inv[i] = vp;
        }
    }

    //class variables
protected:
    SpMat _A;
    int _n;
    int _m;
    permutation _tau_inv;
    bool _binary;
    string _file_name;
    int _coco;
    data_type _epsilon;
    
public:
    //COMMON METHODS
    void preprocessing(); //only for SFS and LBFS
    void shift();
    typename SFSMatrix::IntVector solve();
    // void solve();
    void subgraph(Block& cc_it, SFSMatrix& G);
    void reorder_graph(WeakLinearOrder& CC, bool reversed);
    WeakLinearOrder SFS();
    WeakLinearOrder Neighborhood (int p); //compute the ordered neighborhood of p
    void partition_refinement(Block& S, Queue& Q, map_queue_iterator& map_queue, map_iterator& map_queue_position);
    void concatenate(WeakLinearOrder& Phi);
    IntVector WeakLinearOrder_to_IntVector(WeakLinearOrder& Phi, bool nominal);
    bool binary(); //check if the given matrix is binary or not

    void SerialRank();
    
    //SFS routines
    void Robinson(IntVector& pi);
    int multisweep(IntVector& pi_opt);
    bool isReversed(IntVector& sigma_inv, IntVector& sigma_loop); //cehck if two permutation are reversed
    data_type isEpsilon_Robinson(); //check if the given matrix is epsilon Robinson (only connected comp)
    data_type check_isEpsilon_Robinson();
    
    //Atkins Routines
    void Atkins(IntVector& pi, WeakLinearOrder& CC);
    void spectral(IntVector& pi); //sort fiedler vector
    void repeated (vec& Fiedler, IntVector& pi);
    
    
    //debugging tools
    void print_adjacency_list_tau(SpMat& A);
    void print_adjacency_list_nominal();
    void print_permutation(const permutation& sigma_inv);
    void print_permutation_reversed(const permutation& sigma_inv);
    void print_queue(Queue& Q);
    void print_neighborhood(WeakLinearOrder& N);
    void check_tau();
    void print_vertices_properties();
    void print_graph_properties();
    void print_refine_tau(WeakLinearOrder& phi);
    void print_refine_nominal(WeakLinearOrder& phi);
    bool is_symmetric(SpMat& A);
    bool is_permutation(IntVector& pi);
    bool is_Robinson(SpMat& A);
    int get_n(){return _n;};
    
    
    //printing files
    void print_permutation(IntVector& pi, string& file);
    void print_log(clock_t& t);
    void print_ordered_matrix (IntVector& pi, string& file);
};

#endif // MATRIX_H
