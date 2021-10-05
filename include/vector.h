/**
 * \file include/shg/vector.h
 * Vector class and functions.
 * \date Created on 9 November 2013.
 */

#ifndef SHG_VECTOR_H
#define SHG_VECTOR_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <initializer_list>
#include <istream>
#include <iterator>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

namespace SHG {

/**
 * \defgroup vector_and_matrix Vector and Matrix
 *
 * Classes Vector and Matrix and related functions.
 *
 * \{
 */

/**
 * %Vector class.
 */
    template <class T>
    class Vector {
    public:
        typedef T value_type;
        typedef std::size_t size_type;
        typedef T* iterator;
        typedef const T* const_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator>
                const_reverse_iterator;

        /** \name Construct/destroy. */
        /** \{ */

        /** Constructs an empty vector. */
        inline Vector();

        /** Constructs a vector with \a n value-initialized elements. */
        inline explicit Vector(std::size_t n);

        /** Constructs a vector with \a n elements equal to \a a. */
        inline Vector(std::size_t n, const T& a);

        /**
         * Constructs a vector initialized to the first \a n elements of
         * \a *a.
         */
        inline Vector(std::size_t n, const T* a);

        /**
         * Constructs a vector with the contents of the initializer list
         * \a il.
         */
        inline Vector(std::initializer_list<T> il);

        /** Copy constructor. */
        inline Vector(const Vector& v);

        /** Move constructor. */
        inline Vector(Vector&& v) noexcept;

        /** Constructs a vector from std::vector \a v. */
        explicit Vector(const std::vector<T>& v);

        /** Constructs a vector from std::basic_string \a s. */
        explicit Vector(const std::basic_string<T>& s);

        /** Constructs a vector from std::valarray \a v. */
        explicit Vector(const std::valarray<T>& v);

        /** Destructor. */
        inline virtual ~Vector();

        /** \} */
        /** \name Assignment. */
        /** \{ */

        /** Copy assignment. */
        Vector& operator=(const Vector& v);

        /** Move assignment. */
        inline Vector& operator=(Vector&& v) noexcept(
        noexcept(std::is_nothrow_move_constructible<T>::value &&
                 std::is_nothrow_move_assignable<T>::value));

        /** Assigns to all elements of this vector value \a a. */
        inline Vector& operator=(const T& a);

        /** Assigns elements of initializer list \a il to this vector. */
        inline Vector& operator=(std::initializer_list<T> il);

        /** Assigns to this vector elements of std::vector \a v. */
        Vector& operator=(const std::vector<T>& v);

        /** Assigns to this vector elements of std::basic_string \a s. */
        Vector& operator=(const std::basic_string<T>& s);

        /** Assigns to this vector elements of std::valarray \a v. */
        Vector& operator=(const std::valarray<T>& v);

        /** \} */
        /** \name Element access. */
        /** \{ */

        /** Returns a reference to the \a i-th element. */
        inline T& operator[](std::size_t i);

        /** Returns a constant reference to the \a i-th element. */
        inline const T& operator[](std::size_t i) const;

        /** \copydoc operator[](std::size_t) */
        inline T& operator()(std::size_t i);

        /** \copydoc operator[](std::size_t) const */
        inline const T& operator()(std::size_t i) const;

        /**
         * Returns a reference to the \a i-th element.
         * \exception std::out_of_range if !(i < size())
         */
        T& at(std::size_t i);

        /**
         * Returns a constant reference to the \a i-th element.
         * \exception std::out_of_range if !(i < size())
         */
        const T& at(std::size_t i) const;

        /** \} */
        /** \name Member functions. */
        /** \{ */

        /** Returns the number of elements in the vector. */
        inline std::size_t size() const;

        /**
         * Changes the length of the vector to \a n. All elements are
         * value-initialized.
         */
        void resize(std::size_t n);

        /**
         * Resizes the vector to \a n elements and assigns \a a to each
         * element.
         */
        void assign(std::size_t n, const T& a);

        /**
         * Returns a pointer to the memory block where the data of this
         * vector are stored.
         */
        inline T* c_vec();

        /**
         * Returns a constant pointer to the memory block where the data
         * of this vector are stored.
         */
        inline const T* c_vec() const;

        /**
         * \a *this obtains the value of \a v, \a v obtains the value of
         * \a *this.
         */
        inline void swap(Vector& v) noexcept(
        noexcept(std::is_nothrow_move_constructible<T>::value &&
                 std::is_nothrow_move_assignable<T>::value));

        /** Conversion operator to std::vector. */
        inline explicit operator std::vector<T>() const;

        /** Conversion operator to std::basic_string. */
        inline explicit operator std::basic_string<T>() const;

        /** Conversion operator to std::valarray. */
        inline explicit operator std::valarray<T>() const;

        /**
         * Writes this vector to a binary stream. Check f.fail() to check
         * if the operation was successful.
         */
        void write(std::ostream& f) const;

        /**
         * Reads this vector from a binary stream. Check f.fail() to
         * check if the operation was successful. If the operation fails,
         * the vector remains unchanged.
         */
        void read(std::istream& f);

        /** \} */
        /** \name Iterators. */
        /** \{ */

        inline iterator begin();
        inline iterator end();
        inline const_iterator begin() const;
        inline const_iterator end() const;
        inline const_iterator cbegin() const;
        inline const_iterator cend() const;
        inline reverse_iterator rbegin();
        inline reverse_iterator rend();
        inline const_reverse_iterator rbegin() const;
        inline const_reverse_iterator rend() const;
        inline const_reverse_iterator crbegin() const;
        inline const_reverse_iterator crend() const;

        /** \} */

    protected:
        inline void allocate(std::size_t n);
        inline void deallocate();

        T* v_;          /**< data */
        std::size_t n_; /**< number of elements */
    };

/** \name Vector typedefs for some fundamental types. */
/** \{ */

    typedef Vector<bool> Vecbool;
    typedef Vector<char> Vecchar;
    typedef Vector<short> Vecshort;
    typedef Vector<int> Vecint;
    typedef Vector<long> Veclong;
    typedef Vector<long long> Veclonglong;
    typedef Vector<float> Vecfloat;
    typedef Vector<double> Vecdouble;
    typedef Vector<long double> Veclongdouble;
    typedef Vector<std::complex<float>> Veccomplexfloat;
    typedef Vector<std::complex<double>> Veccomplex;
    typedef Vector<std::complex<long double>> Veccomplexlongdouble;

/** \} */ /* vector typedefs for some fundamental types */

/** \name Vector non-member functions. */
/** \{ */

/**
 * Compares two vectors.
 *
 * \returns true if the two vectors have the same size and a[i] ==
 * b[i] for all i in [0, a.size()), false otherwise.
 */
    template <class T>
    inline bool equal(const Vector<T>& a, const Vector<T>& b);

/**
 * \copydoc equal(const Vector<T>& a, const Vector<T>& b)
 */
    template <class T>
    inline bool operator==(const Vector<T>& a, const Vector<T>& b);

/**
 * Returns the sum of all the elements of the vector. Operator \a +=
 * must be defined for the type \a T. If the vector has length 0, the
 * behaviour is undefined.
 */
    template <class T>
    inline T sum(const Vector<T>& v);

/**
 * Returns the minimum value contained in \a *this. The value for a
 * vector of length 0 is undefined.
 */
    template <class T>
    inline T min(const Vector<T>& v);

/**
 * Returns the maximum value contained in \a *this. The value for a
 * vector of length 0 is undefined.
 */
    template <class T>
    inline T max(const Vector<T>& v);

/**
 * Returns the minimum and maximum values contained in \a *this. Uses
 * operator< to compare values. If the vector is of length 0, the
 * result is undefined.
 */
    template <class T>
    inline std::pair<T, T> minmax(const Vector<T>& v);

/**
 * Returns index of the smallest element. Uses operator< to compare
 * values. If the vector is of length 0, 0 is returned. If several
 * elements are the smallest, the smallest index is returned.
 */
    template <class T>
    inline std::size_t minloc(const Vector<T>& v);

/**
 * Returns index of the greatest element. Uses operator< to compare
 * values. If the vector is of length 0, 0 is returned. If several
 * elements are the greatest, the smallest index is returned.
 */
    template <class T>
    inline std::size_t maxloc(const Vector<T>& v);

/**
 * Returns indexes of the smallest element and the greatest element.
 * Uses operator< to compare values. If the vector is of length 0, 0
 * is returned. The first minimum and maximum elements found are
 * returned.
 */
    template <class T>
    inline std::pair<std::size_t, std::size_t> minmaxloc(
            const Vector<T>& v);

/**
 * Changes the length of the vector to 0.
 */
    template <class T>
    inline void clear(Vector<T>& v);

/**
 * Swaps values of two vectors.
 */
    template <class T>
    inline void swap(Vector<T>& a, Vector<T>& b) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value));

/**
 * Sorts the elements into ascending order using operator<.
 */
    template <class T>
    inline void sort(Vector<T>& v);

/**
 * Sorts the elements in the range [first, last) into ascending
 * order using operator<.
 */
    template <class T>
    inline void sort(Vector<T>& v, std::size_t first, std::size_t last);

/**
 * Sorts the elements into descending order using operator<.
 */
    template <class T>
    inline void reverse_sort(Vector<T>& v);

/**
 * Sorts the elements in the range [first, last) into descending
 * order using operator<.
 */
    template <class T>
    inline void reverse_sort(Vector<T>& v, std::size_t first,
                             std::size_t last);

/**
 * Outputs vector to a text stream. First, the size of the vector is
 * output with width zero followed by a newline character. Then all
 * elements are printed one on a line with the width set just before a
 * call to this function. stream.fail() should be checked after
 * return. For example, the code
 *
 * \code
 * Vector<int> v {9, 99, 999};
 * cout << setw(3) << v;
 * \endcode
 *
 * produces
 *
 * <pre>
 * 3
 *   9
 *  99
 * 999
 * </pre>
 */
    template <class T>
    std::ostream& operator<<(std::ostream& stream, const Vector<T>& v);

/**
 * Inputs vector from a text stream. The input should be as output
 * from std::ostream& operator<<(std::ostream&, const Vector<T>&).
 * stream.fail() should be checked after return. If the operation
 * failed, previous content of \a v is untouched.
 */
    template <class T>
    std::istream& operator>>(std::istream& stream, Vector<T>& v);

/**
 * Outputs vector to a text stream in the form of initializer list.
 * fail() should be checked after return.
 */
    template <class T>
    void print(const Vector<T>& v, std::ostream& stream);

/**
 * \copydoc Vector::write()
 */
    template <class T>
    inline void write(const Vector<T>& v, std::ostream& f);

/**
 * \copydoc Vector::read()
 */
    template <class T>
    inline void read(Vector<T>& v, std::istream& f);

/**
 * Returns the maximum norm distance between two vectors. The norm is
 * defined as \f$ \max \{|a_i - b_i| \colon 0 \leq i < n \} \f$, where
 * \f$n\f$ is a common size of the vectors. If vectors are of length
 * zero, zero is returned.
 *
 * \exception std::invalid_argument if !(a.size() == b.size())
 */
    template <class T>
    T maximum_norm_distance(const Vector<T>& a, const Vector<T>& b);

/**
 * Constructs a vector containing \a n terms of arithmetic sequence
 * with initial term \a a and difference \a d. Operators \a + and \a *
 * must be defined for the type \a T.
 */
    template <class T>
    Vector<T> arithmetic_progression(std::size_t n, const T& a,
                                     const T& d);

/** \} */ /* end of vector non-member functions */

/** \} */ /* end of group vector_and_matrix */

}  // namespace SHG

#include <vector-inl.h>

#endif