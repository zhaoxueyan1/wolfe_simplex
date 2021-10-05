/**
 * \file include/shg/vector-inl.h
 * Implementation of inline functions and templates from vector.h.
 */

#ifndef SHG_VECTOR_INL_H
#define SHG_VECTOR_INL_H

namespace SHG {

    template <class T>
    Vector<T>::Vector() : v_(nullptr), n_(0) {}

    template <class T>
    Vector<T>::Vector(std::size_t n)
            : v_(n > 0 ? new T[n] : nullptr), n_(n) {}

    template <class T>
    Vector<T>::Vector(std::size_t n, const T& a) : Vector(n) {
        std::fill(begin(), end(), a);
    }

    template <class T>
    Vector<T>::Vector(std::size_t n, const T* a) : Vector(n) {
        std::copy_n(a, n, begin());
    }

    template <class T>
    Vector<T>::Vector(std::initializer_list<T> il) : Vector(il.size()) {
        std::copy(il.begin(), il.end(), begin());
    }

    template <class T>
    Vector<T>::Vector(const Vector& v) : Vector(v.size()) {
        std::copy(v.begin(), v.end(), begin());
    }

    template <class T>
    Vector<T>::Vector(Vector&& v) noexcept : v_(v.v_), n_(v.n_) {
        v.v_ = nullptr;
        v.n_ = 0;
    }

    template <class T>
    Vector<T>::Vector(const std::vector<T>& v) : Vector(v.size()) {
        std::copy(v.begin(), v.end(), begin());
    }

    template <class T>
    Vector<T>::Vector(const std::basic_string<T>& s) : Vector(s.size()) {
        s.copy(c_vec(), s.size());
    }

    template <class T>
    Vector<T>::Vector(const std::valarray<T>& v) : Vector(v.size()) {
        std::copy(std::begin(v), std::end(v), begin());
    }

    template <class T>
    Vector<T>::~Vector() {
        deallocate();
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(const Vector& v) {
        if (this != &v) {
            resize(v.size());
            std::copy(v.begin(), v.end(), begin());
        }
        return *this;
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(Vector&& v) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value)) {
        swap(v);
        return *this;
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(const T& a) {
        std::fill(begin(), end(), a);
        return *this;
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(std::initializer_list<T> il) {
        resize(il.size());
        std::copy(il.begin(), il.end(), begin());
        return *this;
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(const std::vector<T>& v) {
        resize(v.size());
        std::copy(v.begin(), v.end(), begin());
        return *this;
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(const std::basic_string<T>& s) {
        resize(s.size());
        s.copy(c_vec(), s.size());
        return *this;
    }

    template <class T>
    Vector<T>& Vector<T>::operator=(const std::valarray<T>& v) {
        resize(v.size());
        std::copy(std::begin(v), std::end(v), begin());
        return *this;
    }

    template <class T>
    T& Vector<T>::operator[](std::size_t i) {
        return v_[i];
    }

    template <class T>
    const T& Vector<T>::operator[](std::size_t i) const {
        return v_[i];
    }

    template <class T>
    T& Vector<T>::operator()(std::size_t i) {
        return v_[i];
    }

    template <class T>
    const T& Vector<T>::operator()(std::size_t i) const {
        return v_[i];
    }

    template <class T>
    T& Vector<T>::at(std::size_t i) {
        if (i < size())
            return v_[i];
        throw std::out_of_range(__func__);
    }

    template <class T>
    const T& Vector<T>::at(std::size_t i) const {
        if (i < size())
            return v_[i];
        throw std::out_of_range(__func__);
    }

    template <class T>
    std::size_t Vector<T>::size() const {
        return n_;
    }

    template <class T>
    void Vector<T>::resize(std::size_t n) {
        if (n != n_) {
            deallocate();
            allocate(n); /* On exception, v_ == nullptr, n_ == 0. */
        }
    }

    template <class T>
    void Vector<T>::assign(std::size_t n, const T& a) {
        resize(n);
        std::fill(begin(), end(), a);
    }

    template <class T>
    T* Vector<T>::c_vec() {
        return v_;
    }

    template <class T>
    const T* Vector<T>::c_vec() const {
        return v_;
    }

    template <class T>
    void Vector<T>::swap(Vector& v) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value)) {
        std::swap(v_, v.v_);
        std::swap(n_, v.n_);
    }

    template <class T>
    Vector<T>::operator std::vector<T>() const {
        return std::vector<T>(begin(), end());
    }

    template <class T>
    Vector<T>::operator std::basic_string<T>() const {
        return std::basic_string<T>(c_vec(), size());
    }

    template <class T>
    Vector<T>::operator std::valarray<T>() const {
        return std::valarray<T>(c_vec(), size());
    }

    template <class T>
    void Vector<T>::write(std::ostream& f) const {
        const std::size_t n = size();
        f.write(reinterpret_cast<const char*>(&n), sizeof n);
        if (n > 0)
            f.write(reinterpret_cast<const char*>(c_vec()),
                    static_cast<std::streamsize>(n) * sizeof(T));
    }

    template <class T>
    void Vector<T>::read(std::istream& f) {
        std::size_t n;
        f.read(reinterpret_cast<char*>(&n), sizeof n);
        if (f.fail())
            return;
        Vector<T> w(n);
        if (n > 0) {
            const std::streamsize s =
                    static_cast<std::streamsize>(n) * sizeof(T);
            f.read(reinterpret_cast<char*>(w.c_vec()), s);
            if (f.fail())
                return;
        }
        swap(w);
    }

    template <class T>
    typename Vector<T>::iterator Vector<T>::begin() {
        return v_;
    }

    template <class T>
    typename Vector<T>::iterator Vector<T>::end() {
        return v_ + n_;
    }

    template <class T>
    typename Vector<T>::const_iterator Vector<T>::begin() const {
        return v_;
    }

    template <class T>
    typename Vector<T>::const_iterator Vector<T>::end() const {
        return v_ + n_;
    }

    template <class T>
    typename Vector<T>::const_iterator Vector<T>::cbegin() const {
        return v_;
    }

    template <class T>
    typename Vector<T>::const_iterator Vector<T>::cend() const {
        return v_ + n_;
    }

    template <class T>
    typename Vector<T>::reverse_iterator Vector<T>::rbegin() {
        return reverse_iterator(end());
    }

    template <class T>
    typename Vector<T>::reverse_iterator Vector<T>::rend() {
        return reverse_iterator(begin());
    }

    template <class T>
    typename Vector<T>::const_reverse_iterator Vector<T>::rbegin() const {
        return const_reverse_iterator(end());
    }

    template <class T>
    typename Vector<T>::const_reverse_iterator Vector<T>::rend() const {
        return const_reverse_iterator(begin());
    }

    template <class T>
    typename Vector<T>::const_reverse_iterator Vector<T>::crbegin()
    const {
        return const_reverse_iterator(end());
    }

    template <class T>
    typename Vector<T>::const_reverse_iterator Vector<T>::crend() const {
        return const_reverse_iterator(begin());
    }

    template <class T>
    void Vector<T>::allocate(std::size_t n) {
        v_ = n > 0 ? new T[n] : nullptr;
        n_ = n;
    }

    template <class T>
    void Vector<T>::deallocate() {
        delete[] v_;
        v_ = nullptr;
        n_ = 0;
    }

    template <class T>
    bool equal(const Vector<T>& a, const Vector<T>& b) {
        return std::equal(a.begin(), a.end(), b.begin(), b.end());
    }

    template <class T>
    bool operator==(const Vector<T>& a, const Vector<T>& b) {
        return equal(a, b);
    }

    template <class T>
    T sum(const Vector<T>& v) {
        return std::accumulate(v.begin(), v.end(), 0);
    }

    template <class T>
    T min(const Vector<T>& v) {
        return *std::min_element(v.begin(), v.end());
    }

    template <class T>
    T max(const Vector<T>& v) {
        return *std::max_element(v.begin(), v.end());
    }

    template <class T>
    std::pair<T, T> minmax(const Vector<T>& v) {
        const std::pair<typename Vector<T>::const_iterator,
                typename Vector<T>::const_iterator>
                i = std::minmax_element(v.begin(), v.end());
        return std::make_pair(*i.first, *i.second);
    }

    template <class T>
    std::size_t minloc(const Vector<T>& v) {
        return std::distance(v.begin(),
                             std::min_element(v.begin(), v.end()));
    }

    template <class T>
    std::size_t maxloc(const Vector<T>& v) {
        return std::distance(v.begin(),
                             std::max_element(v.begin(), v.end()));
    }

    template <class T>
    std::pair<std::size_t, std::size_t> minmaxloc(const Vector<T>& v) {
        const std::pair<typename Vector<T>::const_iterator,
                typename Vector<T>::const_iterator>
                i = std::minmax_element(v.begin(), v.end());
        return std::make_pair(std::distance(v.begin(), i.first),
                              std::distance(v.begin(), i.second));
    }

    template <class T>
    void clear(Vector<T>& v) {
        v.resize(0);
    }

    template <class T>
    void swap(Vector<T>& a, Vector<T>& b) noexcept(
    noexcept(std::is_nothrow_move_constructible<T>::value &&
             std::is_nothrow_move_assignable<T>::value)) {
        a.swap(b);
    }

    template <class T>
    void sort(Vector<T>& v) {
        std::sort(v.begin(), v.end());
    }

    template <class T>
    void sort(Vector<T>& v, std::size_t first, std::size_t last) {
        std::sort(v.begin() + first, v.begin() + last);
    }

    template <class T>
    void reverse_sort(Vector<T>& v) {
        std::sort(v.begin(), v.end(),
                  [](const T& a, const T& b) { return a > b; });
    }

    template <class T>
    void reverse_sort(Vector<T>& v, std::size_t first, std::size_t last) {
        std::sort(v.begin() + first, v.begin() + last,
                  [](const T& a, const T& b) { return a > b; });
    }

    template <class T>
    std::ostream& operator<<(std::ostream& stream, const Vector<T>& v) {
        const std::streamsize w = stream.width(0);
        stream << v.size() << '\n';
        for (const auto& e : v) {
            stream.width(w);
            stream << e << '\n';
        }
        return stream;
    }

    template <class T>
    std::istream& operator>>(std::istream& stream, Vector<T>& v) {
        std::size_t n;
        if ((stream >> n)) {
            Vector<T> w(n);
            for (auto& e : w)
                if (!(stream >> e))
                    return stream;
            w.swap(v);
        }
        return stream;
    }

    template <class T>
    void print(const Vector<T>& v, std::ostream& stream) {
        stream.width(0);
        stream << "{";
        if (v.size() > 0) {
            stream << v[0];
            for (std::size_t i = 1; i < v.size(); i++)
                stream << ", " << v[i];
        }
        stream << "}";
    }

    template <class T>
    void write(const Vector<T>& v, std::ostream& f) {
        v.write(f);
    }

    template <class T>
    void read(Vector<T>& v, std::istream& f) {
        v.read(f);
    }

    template <class T>
    T maximum_norm_distance(const Vector<T>& a, const Vector<T>& b) {
        if (a.size() == b.size()) {
            T d, max = 0;
            for (std::size_t i = 0; i < a.size(); i++)
                if ((d = std::abs(a[i] - b[i])) > max)
                    max = d;
            return max;
        }
        throw std::invalid_argument(__func__);
    }

    template <class T>
    Vector<T> arithmetic_progression(std::size_t n, const T& a,
                                     const T& d) {
        Vector<T> v(n);
        std::size_t i = 0;
        std::generate_n(v.begin(), n,
                        [a, d, &i] { return a + i++ * d; });
        return v;
    }

}  // namespace SHG

#endif