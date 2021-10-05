/**
 * \file include/shg/matrix.h
 *  Matrix class and functions.
 * \date Created on 9 November 2013.
 */

#ifndef SHG_MATRIX_H
#define SHG_MATRIX_H

#include <cstdlib>
#include <ios>
#include <limits>
#include <utility>
#include <vector.h>

namespace SHG {

/** \addtogroup vector_and_matrix */
/** \{ */

/**
 * %Matrix class.
 */
    template <class T>
    class Matrix {
    public:
        typedef T value_type;
        typedef std::size_t size_type;

        /** \name Construct/destroy. */
        /** \{ */

        /** Constructs an empty matrix. */
        inline Matrix();

        /**
         * Constructs an <em>m &times; n</em> matrix with
         * value-initialized elements.
         */
        Matrix(std::size_t m, std::size_t n);

        /**
         * Constructs an <em>m &times; n</em> matrix with all elements
         * equal to \a a.
         */
        inline Matrix(std::size_t m, std::size_t n, const T& a);

        /**
         * Constructs from C memory block.
         */
        inline Matrix(std::size_t m, std::size_t n, const T* a);

        /**
         * Constructs from C two-dimensional array.
         */
        Matrix(std::size_t m, std::size_t n, const T* const* a);

        /**
         * Constructs from Vector.
         */
        Matrix(std::size_t m, std::size_t n, const Vector<T>& v);

        /**
         * Move constructor from Vector.
         */
        Matrix(std::size_t m, std::size_t n, Vector<T>&& v);

        /**
         * Constructs a matrix from an initializer list. Matrix elements
         * are initialized by rows. If the length of the list is less
         * than \a m &times; \a n, given elements are used cyclically. If
         * the list is empty, nothing happens.
         */
        Matrix(std::size_t m, std::size_t n,
               std::initializer_list<T> il);

        /**
         * Copy constructor.
         */
        Matrix(const Matrix& a);

        /**
         * Move constructor.
         */
        inline Matrix(Matrix&& a) noexcept;

        /**
         * Destructor.
         */
        virtual ~Matrix();

        /** \} */
        /** \name Assignment. */
        /** \{ */

        /**
         * Copy assignment.
         */
        Matrix& operator=(const Matrix& a);

        /**
         * Move assignment.
         */
        inline Matrix& operator=(Matrix&& a) noexcept(
        noexcept(std::is_nothrow_move_constructible<T>::value &&
                 std::is_nothrow_move_assignable<T>::value));

        /**
         * Assigns \a a to all elements
         */
        inline Matrix& operator=(const T& a);

        /**
         * Assigns an initializer list. Matrix elements are assigned by
         * rows. If the length of the list is less than \a m &times; \a
         * n, given elements are used cyclically. If the list is empty,
         * nothing happens.
         */
        Matrix& operator=(std::initializer_list<T> il);

        /** \} */
        /** \name Element access. */
        /** \{ */

        /**
         * Returns a pointer to the \a i-th row.
         */
        inline T* operator[](std::size_t i);

        /**
         * Returns a constant pointer to the \a i-th row.
         */
        inline const T* operator[](std::size_t i) const;

        /**
         * Returns element in row \a i and column \a j.
         */
        inline T& operator()(std::size_t i, std::size_t j);

        /**
         * Returns constant element in row \a i and column \a j.
         */
        inline const T& operator()(std::size_t i, std::size_t j) const;

        /**
         * Returns a reference to the element <em>(i, j)</em> with range
         * checking.
         *
         * \exception std::out_of_range if !(i < nrows() && j < ncols())
         */
        T& at(std::size_t i, std::size_t j);

        /**
         * Returns a constant reference to the element <em>(i, j)</em>
         * with range checking.
         *
         * \exception std::out_of_range if !(i < nrows() && j < ncols())
         */
        const T& at(std::size_t i, std::size_t j) const;

        /** \} */
        /** \name Member functions. */
        /** \{ */

        /**
         * Returns the number of rows in the matrix.
         */
        inline std::size_t nrows() const;

        /**
         * Returns the number of columns in the matrix.
         */
        inline std::size_t ncols() const;

        /**
         * Changes the dimensions of the matrix to \a m by \a n. The
         * elements are value-initialized.
         */
        void resize(std::size_t m, std::size_t n);

        /**
         * Changes the dimensions of the matrix to \a m by \a n and then
         * assigns to each element the value \a a.
         */
        void assign(std::size_t m, std::size_t n, const T& a);

        /**
         * Returns a pointer to the memory block where the data of matrix
         * are stored.
         */
        inline T* c_vec();

        /**
         * Returns a constant pointer to the memory block where the data
         * of matrix are stored.
         */
        inline const T* c_vec() const;

        /**
         * Returns this matrix as a C-style matrix.
         */
        inline T* const* c_mat();

        /**
         * Returns this matrix as a constant C-style matrix.
         */
        inline T const* const* c_mat() const;

        /**
         * Returns this matrix as a Vector.
         */
        inline Vector<T>& vector();

        /**
         * Returns this matrix as a constant Vector.
         */
        inline const Vector<T>& vector() const;

        /**
         * Exchanges values of \a *this and \a a.
         */
        void swap(Matrix& a) noexcept(
        noexcept(std::is_nothrow_move_constructible<T>::value &&
                 std::is_nothrow_move_assignable<T>::value));

        /** \} */

    protected:
        Vector<T> v_;   /**< data */
        Vector<T*> p_;  /**< pointers to rows */
        std::size_t m_; /**< the number of rows */
        std::size_t n_; /**< the number of columns */
    };

/** \name Matrix typedefs for some fundamental types. */
/** \{ */

    typedef Matrix<bool> Matbool;
    typedef Matrix<char> Matchar;
    typedef Matrix<short> Matshort;
    typedef Matrix<int> Matint;
    typedef Matrix<long> Matlong;
    typedef Matrix<long long> Matlonglong;
    typedef Matrix<float> Matfloat;
    typedef Matrix<double> Matdouble;
    typedef Matrix<long double> Matlongdouble;

/** \} */ /* matrix typedefs for some fundamental types */

/** \name Matrix non-member functions. */
/** \{ */

/**
 * Compares two matrices.
 *
 * \returns true if the two matrices have the same dimensions and a(i,
 * j) == b(i, j) for all i, j, false otherwise.
 */
    template <class T>
    bool equal(const Matrix<T>& a, const Matrix<T>& b);

/**
 * \copydoc equal(const Matrix<T>& a, const Matrix<T>& b)
 */
    template <class T>
    inline bool operator==(const Matrix<T>& a, const Matrix<T>& b);

/**
 * Returns the sum of all the elements of the matrix. Operator \a +=
 * must be defined for the type \a T. If the matrix is empty, the
 * behaviour is undefined.
 */
    template <class T>
    inline T sum(const Matrix<T>& a);

/**
 * Returns the minimum value contained in \a *this. The value for an
 * empty matrix is undefined.
 */
    template <class T>
    inline T min(const Matrix<T>& a);

/**
 * Returns the maximum value contained in \a *this. The value for an
 * empty matrix is undefined.
 */
    template <class T>
    inline T max(const Matrix<T>& a);

/**
 * Returns the minimum and maximum values contained in \a *this. The
 * value for an empty matrix is undefined.
 */
    template <class T>
    inline std::pair<T, T> minmax(const Matrix<T>& a);

/**
 * Returns index of the smallest element. Uses operator< to compare
 * values. If the matrix is empty, 0 is returned.
 */
    template <class T>
    inline std::pair<std::size_t, std::size_t> minloc(const Matrix<T>& a);

/**
 * Returns index of the greatest element. Uses operator< to compare
 * values. If the matrix is empty, 0 is returned.
 */
    template <class T>
    inline std::pair<std::size_t, std::size_t> maxloc(const Matrix<T>& a);

/**
 * Returns indexes of the smallest and the greatest element. Uses
 * operator< to compare values. If the matrix is empty, 0 is returned.
 */
    template <class T>
    inline std::pair<std::pair<std::size_t, std::size_t>,
            std::pair<std::size_t, std::size_t>>
    minmaxloc(const Matrix<T>& a);

/**
 * Changes the dimensions of the matrix to 0 by 0.
 */
    template <class T>
    inline void clear(Matrix<T>& a);

/**
 * Exchanges values of \a *this and \a a.
 */
    template <class T>
    inline void swap(Matrix<T>& a, Matrix<T>& b) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value));

/**
 * Outputs matix to a text stream. First, the dimensions of the matrix
 * are output with width zero followed by a newline character. Then
 * all elements are printed with each row on a line with the width set
 * just before a call to this function. stream.fail() should be
 * checked after return.
 */
    template <class T>
    std::ostream& operator<<(std::ostream& stream, const Matrix<T>& a);

/**
 * Inputs matrix from a text stream. The input should be as output
 * from std::ostream& operator<<(std::ostream&, const Matrix<T>&).
 * stream.fail() should be checked after return. If the operation
 * failed, previous content of \a a is untouched.
 */
    template <class T>
    std::istream& operator>>(std::istream& stream, Matrix<T>& a);

/**
 * Outputs matrix to a text stream in the form of initializer list.
 * stream.fail() should be checked after return.
 */
    template <class T>
    void print(const Matrix<T>& a, std::ostream& stream);

/**
 * Writes this matrix to a binary stream. Check f.fail() to check if
 * the operation was successful.
 */
    template <class T>
    void write(const Matrix<T>& a, std::ostream& f);

/**
 * Reads this matrix from a binary stream. Check f.fail() to check if
 * the operation was successful. If the operation fails, the matrix
 * remains unchanged.
 */
    template <class T>
    void read(Matrix<T>& a, std::istream& f);

/**
 * Returns maximum norm distance between two matrices. The norm is
 * defined as \f$ \max \{|a_{ij} - b_{ij}| \colon 0 \leq i < m, 0 \leq
 * j < n \} \f$, where \f$m\f$ and \f$n\f$ are the common dimensions
 * of the matrices. If the matrices are empty, 0 is returned.
 *
 * \exception std::invalid_argument if !(a.nrows() == b.nrows() &&
 * a.ncols() == b.ncols())
 */
    template <class T>
    T maximum_norm_distance(const Matrix<T>& a, const Matrix<T>& b);

/**
 * Returns diagonal matrix with \a c on the main diagonal.
 */
    template <class T>
    Matrix<T> diagonal_matrix(std::size_t n, const T& c = 1);

/**
 * Returns matrix transposed to \a a.
 */
    template <class T>
    Matrix<T> transpose(const Matrix<T>& a);

/**
 * Transposes the square matrix \a a in situ.
 * \returns reference to \a a
 * \exception std::invalid_argument if (a.nrows() != a.ncols())
 */
    template <class T>
    Matrix<T>& transpose_in_situ(Matrix<T>& a);

/**
 * Returns the multiplication of \a a and \a b.
 * \exception std::invalid_argument if a.ncols() != b.nrows()
 */
    template <class T>
    Matrix<T> multiply(const Matrix<T>& a, const Matrix<T>& b);

/**
 * Performs operation \f$ a \leftarrow a \times b \f$. \f$b\f$ must be
 * a square matrix with the number of rows equal to the number of
 * columns in \f$a\f$. The matrices must be different.
 *
 * \exception std::invalid_argument if &a == &b || a.ncols() !=
 * b.nrows() || b.nrows() != b.ncols()
 */
    template <class T>
    void right_multiply_and_assign(Matrix<T>& a, const Matrix<T>& b);

/**
 * Returns the matrix \f$a^Ta\f$.
 */
    template <class T>
    Matrix<T> left_multiply_by_transposition(const Matrix<T>& a);

/**
 * Inverts matrix by Cholesky decomposition. Calculates in situ the
 * inverse matrix for a symmetric and positive definite matrix by
 * Cholesky decomposition. The function does not check if the matrix
 * is symmetric and uses only the upper-right triangle of the matrix.
 *
 * \tparam T floating-point type
 * \param [inout] a matrix to invert
 * \param [in] eps if a number is less than or equal to \a eps, it is
 * treated as 0
 * \exception std::invalid_argument if a.nrows() != a.ncols() || eps <
 * 0 \exception std::range_error if the matrix is not positive
 * definite
 *
 * \note In case of exception elements of the matrix are undefined.
 *
 * \implementation See \cite marciniak-gregulec-kaczmarek-1991, p.
 * 90-92.
 */
    template <class T>
    void cholesky(Matrix<T>& a,
                  T eps = 4 * std::numeric_limits<T>::epsilon());

/**
 * Returns an <em>n &times; n</em> Hilbert matrix. The Hilbert matrix
 * is defined by \f$a_{ij} = 1 / (i + j + 1) \; i, j = 0, \ldots, n -
 * 1\f$.
 */
    template <class T>
    Matrix<T> hilbert_matrix(std::size_t n);

/**
 * Returns the vector \f$av\f$.
 * \exception std::invalid_argument if a.ncols() != v.size()
 */
    template <class T>
    Vector<T> multiply(const Matrix<T>& a, const Vector<T>& v);

/**
 * Returns the vector \f$a^Tv\f$.
 * \exception std::invalid_argument if a.nrows() != v.size()
 */
    template <class T>
    Vector<T> multiply_transposed(const Matrix<T>& a, const Vector<T>& v);

/** \} */ /* end of matrix non-member functions */

/** \} */  // end of group vector_and_matrix

}  // namespace SHG

#include <matrix-inl.h>

#endif