// Two dimensional array class
#ifndef ARRAY2D_H
#define ARRAY2D_H
#include <vector>
#include <cassert>

template <class T>
class Array2D
{
    private:
        int nrow_;
        int ncol_;
        int size_; // number of entries

    public:
        std::vector<T> data_;
        Array2D(int nrow, int ncol, const T& initVal);
        Array2D(const Array2D<T>& i);
        int numRows() const { return nrow_; }
        int numCols() const { return ncol_; }
        int size() const { return size_; }

        // Accessing elements / setting elements
        T& operator()(int i, int j)
        {
            assert( (i>=0) && (i<nrow_));
            assert( (j>=0) && (j<ncol_));
            return data_[i*ncol_ + j];
        }
        T& get(int i, int j)
        { 
            assert( (i>=0) && (i<nrow_));
            assert( (j>=0) && (j<ncol_));
            return data_[i*ncol_ + j];
        }

        // Accessing elements
        T operator()(int i, int j) const {return data_[i*ncol_ +j];}
        T get(int i, int j) const { return data_[i*ncol_ + j];}
};

template <class T>
Array2D<T>::Array2D(const Array2D<T>& i) :
    nrow_(i.nrow_), ncol_(i.ncol_), size_(i.size_),
    data_(i.data_) {};

template <class T>
Array2D<T>::Array2D(int nrow, int ncol, const T& initVal = T())
{
    nrow_ = nrow;
    ncol_ = ncol;
    size_ = nrow*ncol;
    data_ = std::vector<T>(size_, initVal);
    return;
}

#endif

