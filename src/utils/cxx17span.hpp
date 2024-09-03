#pragma once
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <array>

template <typename T>
class cxx17span {
private:
    T* ptr;
    std::size_t sz;
public:
    cxx17span(T* ptr, std::size_t size) :
        ptr(ptr), sz(size) {}

    cxx17span(cxx17span<std::remove_const_t<T>> &vec) :
        ptr(const_cast<T*>(vec.data())), sz(vec.size()) {}

    cxx17span(std::vector<std::remove_const_t<T>> &vec) :
        ptr(vec.data()), sz(vec.size()) {}

    cxx17span(const std::vector<std::remove_const_t<T>> &vec) :
        ptr(const_cast<T*>(vec.data())), sz(vec.size()) {}

    template <std::size_t N>
    cxx17span(std::array<std::remove_const_t<T>, N> &arr) :
        ptr(arr.data()), sz(N) {}

    template <std::size_t N>
    cxx17span(const std::array<std::remove_const_t<T>, N> &arr) :
        ptr(const_cast<T*>(arr.data())), sz(N) {}


    T* data() const { return ptr; }
    std::size_t size() const { return sz; }

    T& operator[](std::size_t index) const {
        if (index >= sz) {
            throw std::out_of_range("Index out of range");
        }
        return ptr[index];
    }
    
    T* begin() { return ptr; }
    T* end() { return ptr + sz; }
    const T* begin() const { return ptr; }
    const T* end() const { return ptr + sz; }

    cxx17span<T> subspan(std::size_t offset, std::size_t count) const {
        if (offset > sz) {
            throw std::out_of_range("Offset out of range");
        }
        if (offset + count > sz) {
            throw std::out_of_range("Subspan exceeds original span size");
        }
        return cxx17span<T>(ptr + offset, count);
    }
};