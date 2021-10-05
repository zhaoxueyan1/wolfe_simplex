/**
 * \file include/shg/matrix-inl.h
 * Implementation of inline functions and templates from matrix.h.
 */

#ifndef SHG_MATRIX_INL_H
#define SHG_MATRIX_INL_H

namespace SHG {

    template <class T>
    Matrix<T>::Matrix() : v_(), p_(), m_(0), n_(0) {}

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n)
            : v_(m * n), p_(n > 0 ? m : 0), m_(m), n_(n) {
        if (m_ > 0 && n_ > 0) {
            T* p = p_[0] = &(v_[0]);
            std::generate(p_.begin() + 1, p_.end(),
                          [&p, this]() { return p += this->n_; });
        } else {
            m_ = n_ = 0;
        }
    }

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n, const T& a)
            : Matrix(m, n) {
        v_ = a;
    }

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n, const T* a)
            : Matrix(m, n) {
        std::copy_n(a, m_ * n_, v_.begin());
    }

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n, const T* const* a)
            : Matrix(m, n) {
        for (std::size_t i = 0; i < m; i++)
            for (std::size_t j = 0; j < n; j++)
                p_[i][j] = a[i][j];
    }

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n, const Vector<T>& v)
            : Matrix(m, n) {
        if (v.size() != m * n)
            throw std::invalid_argument(__func__);
        v_ = v;
    }

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n, Vector<T>&& v)
            : Matrix() {
        if (m > 0 && n > 0) {
            if (v.size() != m * n)
                throw std::invalid_argument(__func__);
            v_ = std::move(v);
            p_.resize(m);
            m_ = m;
            n_ = n;
            T* p = p_[0] = &(v_[0]);
            std::generate(p_.begin() + 1, p_.end(),
                          [&p, this]() { return p += this->n_; });
        } else {
            m_ = n_ = 0;
        }
    }

    template <class T>
    Matrix<T>::Matrix(std::size_t m, std::size_t n,
                      std::initializer_list<T> il)
            : Matrix(m, n) {
        typename std::initializer_list<T>::const_iterator it =
                il.begin();
        if (it != il.end()) {
            for (std::size_t i = 0; i < v_.size(); i++) {
                v_[i] = *it;
                if (++it == il.end())
                    it = il.begin();
            }
        }
    }

    template <class T>
    Matrix<T>::Matrix(const Matrix& a)
            : v_(a.v_), p_(a.p_.size()), m_(a.m_), n_(a.n_) {
        if (m_ > 0 && n_ > 0) {
            T* p = p_[0] = &(v_[0]);
            std::generate(p_.begin() + 1, p_.end(),
                          [&p, this]() { return p += this->n_; });
        }
    }

    template <class T>
    Matrix<T>::Matrix(Matrix&& a) noexcept : Matrix() {
        swap(a);
    }

    template <class T>
    Matrix<T>::~Matrix() {}

    template <class T>
    Matrix<T>& Matrix<T>::operator=(const Matrix& a) {
        if (this != &a) {
            resize(a.m_, a.n_);
            v_ = a.v_;
        }
        return *this;
    }

    template <class T>
    Matrix<T>& Matrix<T>::operator=(Matrix&& a) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value)) {
        swap(a);
        return *this;
    }

    template <class T>
    Matrix<T>& Matrix<T>::operator=(const T& a) {
        v_ = a;
        return *this;
    }

    template <class T>
    Matrix<T>& Matrix<T>::operator=(const std::initializer_list<T> il) {
        typename std::initializer_list<T>::const_iterator it =
                il.begin();
        if (it != il.end()) {
            for (std::size_t i = 0; i < v_.size(); i++) {
                v_[i] = *it;
                if (++it == il.end())
                    it = il.begin();
            }
        }
        return *this;
    }

    template <class T>
    T* Matrix<T>::operator[](std::size_t i) {
        return p_[i];
    }

    template <class T>
    const T* Matrix<T>::operator[](std::size_t i) const {
        return p_[i];
    }

    template <class T>
    T& Matrix<T>::operator()(std::size_t i, std::size_t j) {
        return *(p_[i] + j);
    }

    template <class T>
    const T& Matrix<T>::operator()(std::size_t i, std::size_t j) const {
        return *(p_[i] + j);
    }

    template <class T>
    T& Matrix<T>::at(std::size_t i, std::size_t j) {
        if (i < nrows() && j < ncols())
            return *(p_[i] + j);
        throw std::out_of_range(__func__);
    }

    template <class T>
    const T& Matrix<T>::at(std::size_t i, std::size_t j) const {
        if (i < nrows() && j < ncols())
            return *(p_[i] + j);
        throw std::out_of_range(__func__);
    }

    template <class T>
    std::size_t Matrix<T>::nrows() const {
        return m_;
    }

    template <class T>
    std::size_t Matrix<T>::ncols() const {
        return n_;
    }

    template <class T>
    void Matrix<T>::resize(std::size_t m, std::size_t n) {
        if (m != m_ || n != n_) {
            v_.resize(0);
            p_.resize(0);
            m_ = 0;
            n_ = 0;
            Matrix<T> a(m, n);
            swap(a);
        }
    }

    template <class T>
    void Matrix<T>::assign(std::size_t m, std::size_t n, const T& a) {
        resize(m, n);
        v_ = a;
    }

    template <class T>
    T* Matrix<T>::c_vec() {
        return v_.c_vec();
    }

    template <class T>
    const T* Matrix<T>::c_vec() const {
        return v_.c_vec();
    }

    template <class T>
    T* const* Matrix<T>::c_mat() {
        return p_.c_vec();
    }

    template <class T>
    T const* const* Matrix<T>::c_mat() const {
        return p_.c_vec();
    }

    template <class T>
    Vector<T>& Matrix<T>::vector() {
        return v_;
    }

    template <class T>
    const Vector<T>& Matrix<T>::vector() const {
        return v_;
    }

    template <class T>
    void Matrix<T>::swap(Matrix& a) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value)) {
        v_.swap(a.v_);
        p_.swap(a.p_);
        std::swap(m_, a.m_);
        std::swap(n_, a.n_);
    }

    template <class T>
    bool equal(const Matrix<T>& a, const Matrix<T>& b) {
        return (a.nrows() == b.nrows()) && (a.ncols() == b.ncols()) &&
               equal(a.vector(), b.vector());
    }

    template <class T>
    bool operator==(const Matrix<T>& a, const Matrix<T>& b) {
        return equal(a, b);
    }

    template <class T>
    T sum(const Matrix<T>& a) {
        return sum(a.vector());
    }

    template <class T>
    T min(const Matrix<T>& a) {
        return min(a.vector());
    }

    template <class T>
    T max(const Matrix<T>& a) {
        return max(a.vector());
    }

    template <class T>
    std::pair<T, T> minmax(const Matrix<T>& a) {
        return minmax(a.vector());
    }

    template <class T>
    std::pair<std::size_t, std::size_t> minloc(const Matrix<T>& a) {
        const auto d = std::div(
                static_cast<std::ptrdiff_t>(minloc(a.vector())), a.ncols());
        return std::make_pair(d.quot, d.rem);
    }

    template <class T>
    std::pair<std::size_t, std::size_t> maxloc(const Matrix<T>& a) {
        const auto d = std::div(
                static_cast<std::ptrdiff_t>(maxloc(a.vector())), a.ncols());
        return std::make_pair(d.quot, d.rem);
    }

    template <class T>
    std::pair<std::pair<std::size_t, std::size_t>,
    std::pair<std::size_t, std::size_t>>
    minmaxloc(const Matrix<T>& a) {
        const std::pair<std::size_t, std::size_t> p =
                minmaxloc(a.vector());

        const auto d1 =
                std::div(static_cast<std::ptrdiff_t>(p.first), a.ncols());
        const auto d2 =
                std::div(static_cast<std::ptrdiff_t>(p.second), a.ncols());
        return std::make_pair(std::make_pair(d1.quot, d1.rem),
                              std::make_pair(d2.quot, d2.rem));
    }

    template <class T>
    void clear(Matrix<T>& a) {
        a.resize(0, 0);
    }

    template <class T>
    void swap(Matrix<T>& a, Matrix<T>& b) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value)) {
        a.swap(b);
    }

    template <class T>
    std::ostream& operator<<(std::ostream& stream, const Matrix<T>& a) {
        const std::streamsize w = stream.width(0);
        stream << a.nrows() << ' ' << a.ncols() << '\n';
        for (std::size_t i = 0; i < a.nrows(); i++) {
            if (a.ncols() > 0) {
                stream.width(w);
                stream << a[i][0];
                for (std::size_t j = 1; j < a.ncols(); j++) {
                    stream << ' ';
                    stream.width(w);
                    stream << a[i][j];
                }
            }
            stream << '\n';
        }
        return stream;
    }

    template <class T>
    std::istream& operator>>(std::istream& stream, Matrix<T>& a) {
        std::size_t m, n;
        stream >> m >> n;
        if (stream.fail())
            return stream;
        Matrix<T> b(m, n);
        for (std::size_t i = 0; i < b.nrows(); i++)
            for (std::size_t j = 0; j < b.ncols(); j++) {
                stream >> b(i, j);
                if (stream.fail())
                    return stream;
            }
        a.swap(b);
        return stream;
    }

    template <class T>
    void print(const Matrix<T>& a, std::ostream& stream) {
        print(a.vector(), stream);
    }

    template <class T>
    void write(const Matrix<T>& a, std::ostream& f) {
        const std::size_t m = a.nrows();
        const std::size_t n = a.ncols();
        f.write(reinterpret_cast<const char*>(&m), sizeof m);
        f.write(reinterpret_cast<const char*>(&n), sizeof n);
        write(a.vector(), f);
    }

    template <class T>
    void read(Matrix<T>& a, std::istream& f) {
        std::size_t m, n;
        f.read(reinterpret_cast<char*>(&m), sizeof m);
        f.read(reinterpret_cast<char*>(&n), sizeof n);
        if (f.fail())
            return;
        Vector<T> v;
        read(v, f);
        if (f.fail())
            return;
        Matrix<T> b(m, n, std::move(v));
        a.swap(b);
    }

    template <class T>
    T maximum_norm_distance(const Matrix<T>& a, const Matrix<T>& b) {
        if (a.nrows() == b.nrows() && a.ncols() == b.ncols())
            return maximum_norm_distance(a.vector(), b.vector());
        throw std::invalid_argument(__func__);
    }

    template <class T>
    Matrix<T> diagonal_matrix(std::size_t n, const T& c) {
        Matrix<T> a(n, n, static_cast<T>(0));
        for (std::size_t i = 0; i < n; i++)
            a[i][i] = c;
        return a;
    }

    template <class T>
    Matrix<T> transpose(const Matrix<T>& a) {
        Matrix<T> b(a.ncols(), a.nrows());
        for (std::size_t i = 0; i < a.nrows(); i++)
            for (std::size_t j = 0; j < a.ncols(); j++)
                b(j, i) = a(i, j);
        return b;
    }

    template <class T>
    Matrix<T>& transpose_in_situ(Matrix<T>& a) {
        if (a.nrows() != a.ncols())
            throw std::invalid_argument(__func__);
        for (std::size_t i = 0; i < a.nrows(); i++)
            for (std::size_t j = i + 1; j < a.ncols(); j++)
                std::swap(a(i, j), a(j, i));
        return a;
    }

    template <class T>
    Matrix<T> multiply(const Matrix<T>& a, const Matrix<T>& b) {
        if (a.ncols() != b.nrows())
            throw std::invalid_argument(__func__);
        Matrix<T> c(a.nrows(), b.ncols());
        for (std::size_t i = 0; i < a.nrows(); i++) {
            const T* const p = a[i];
            T* const q = c[i];
            for (std::size_t j = 0; j < b.ncols(); j++) {
                T s = 0;
                for (std::size_t k = 0; k < a.ncols(); k++)
                    s += p[k] * b[k][j];
                q[j] = s;
            }
        }
        return c;
    }

    template <class T>
    void right_multiply_and_assign(Matrix<T>& a, const Matrix<T>& b) {
        if (&a == &b || a.ncols() != b.nrows() || b.nrows() != b.ncols())
            throw std::invalid_argument(__func__);
        Vector<T> z(a.ncols());
        std::size_t i, j, k;
        T s, *p;
        for (i = 0; i < a.nrows(); i++) {
            p = a[i];
            for (j = 0; j < a.ncols(); j++)
                z[j] = p[j];
            for (j = 0; j < a.ncols(); j++) {
                s = 0;
                for (k = 0; k < a.ncols(); k++)
                    s += z[k] * b[k][j];
                p[j] = s;
            }
        }
    }

    template <class T>
    Matrix<T> left_multiply_by_transposition(const Matrix<T>& a) {
        Matrix<T> b(a.ncols(), a.ncols());
        std::size_t i, j, k;
        T s;
        for (i = 0; i < a.ncols(); i++) {
            for (j = 0; j < i; j++) {
                s = 0;
                for (k = 0; k < a.nrows(); k++)
                    s += a(k, i) * a(k, j);
                b(i, j) = b(j, i) = s;
            }
            s = 0;
            for (k = 0; k < a.nrows(); k++)
                s += [](T x) { return x * x; }(a(k, i));
            b(i, i) = s;
        }
        return b;
    }

    template <class T>
    void cholesky(Matrix<T>& a, T eps) {
        static_assert(std::is_floating_point<T>::value,
                      "cholesky requires floating-point matrix");
        if (a.nrows() != a.ncols() || eps < 0)
            throw std::invalid_argument(__func__);
        const std::size_t n = a.nrows();
        std::size_t i, j, k;
        T x, z;

        // Find the matrix L such that LL^T is equal to the given matrix
        // and the upper-right triangle of L contains only zeros. Put it
        // in the lower-left triangle of a, inverting the diagonal
        // elements of L.

        for (i = 0; i < n; i++) {
            z = a[i][i];
            for (k = 0; k < i; k++)
                z -= [](T x) { return x * x; }(a[i][k]);
            if (z <= eps)
                throw std::range_error(__func__);
            a[i][i] = z = 1.0 / std::sqrt(z);
            if (!std::isfinite(z))
                throw std::range_error(__func__);
            for (j = i + 1; j < n; j++) {
                x = a[i][j];
                for (k = 0; k < i; k++)
                    x -= a[j][k] * a[i][k];
                a[j][i] = x * z;
            }
        }

        // Calculate L^{-1} and put it in the lower-left triangle of a.

        for (i = 1; i < n; i++)
            for (j = 0; j < i; j++) {
                x = 0.0;
                for (k = j; k < i; k++)
                    x -= a[i][k] * a[k][j];
                a[i][j] = x * a[i][i];
            }

        // Calculate (L^{-1})^T L^{-1} and put it in the whole array a.

        for (i = 0; i < n; i++)
            for (j = i; j < n; j++) {
                x = 0.0;
                for (k = j; k < n; k++)
                    x += a[k][i] * a[k][j];
                a[i][j] = a[j][i] = x;
            }
    }

    template <class T>
    Matrix<T> hilbert_matrix(std::size_t n) {
        static_assert(std::is_floating_point<T>::value,
                      "hilbert_matrix requires floating-point type");
        Matrix<T> h(n, n);
        const std::size_t n2 = 2 * n;
        for (std::size_t k = 1, k1 = 0; k < n2; k++, k1++) {
            const T z = static_cast<T>(1) / k;
            const std::size_t min = k < n ? k : n;
            for (std::size_t i = k - min; i < min; i++)
                h(i, k1 - i) = z;
        }
        return h;
    }

    template <class T>
    Vector<T> multiply(const Matrix<T>& a, const Vector<T>& v) {
        if (a.ncols() != v.size())
            throw std::invalid_argument(__func__);
        Vector<T> w(a.nrows());
        for (std::size_t i = 0; i < w.size(); i++) {
            const T* const p = a[i];
            T s = 0;
            for (std::size_t j = 0; j < v.size(); j++)
                s += p[j] * v[j];
            w[i] = s;
        }
        return w;
    }

    template <class T>
    Vector<T> multiply_transposed(const Matrix<T>& a,
                                  const Vector<T>& v) {
        if (a.nrows() != v.size())
            throw std::invalid_argument(__func__);
        Vector<T> w(a.ncols());
        for (std::size_t i = 0; i < w.size(); i++) {
            T s = 0;
            for (std::size_t j = 0; j < v.size(); j++)
                s += a[j][i] * v[j];
            w[i] = s;
        }
        return w;
    }

}  // namespace SHG

#endif