#ifndef ARRAY_H
#define ARRAY_H

#include <cassert>

template <typename T>
class Array {
private:
    int m_nLength;
    T *m_ptData;

public:
    Array() {
        m_nLength = 0;
        m_ptData = nullptr;
    }

    Array(int nLength) {
        m_nLength = nLength;
        m_ptData = new T[nLength];
    }

    virtual ~Array() {
        delete[] m_ptData;
    }

    void Erase() {
        delete[] m_ptData;
        m_nLength = 0;
        m_ptData = nullptr;
    }

    T& operator[](int nIndex) {
        assert(nIndex >= 0 && nIndex < m_nLength);
        return m_ptData[nIndex];
    }

    void Reallocate(int nNewLength) {
        Erase();
        if (nNewLength <= 0)
            return;

        m_ptData = new T[nNewLength];
        m_nLength = nNewLength;
    }

    void Resize(int nNewLength) {
        if (nNewLength <= 0) {
            Erase();
            return;
        }

        T *ptData = new T[nNewLength];

        if (m_nLength > 0) {
            int nElementsToCopy = (nNewLength > m_nLength) ? m_nLength : nNewLength;
            for (int nIndex = 0; nIndex < nElementsToCopy; ++nIndex) {
                ptData[nIndex] = m_ptData[nIndex];
            }
        }

        delete[] m_ptData;

        m_ptData = ptData;
        m_nLength = nNewLength;
    }

    void InsertBefore(T tValue, int nIndex) {
        assert(nIndex >= 0 && nIndex <= m_nLength);

        T *ptData = new T[m_nLength+1];

        for (int nBefore = 0; nBefore < nIndex; nBefore++) {
            ptData[nBefore] = m_ptData[nBefore];
        }

        ptData[nIndex] = tValue;

        for (int nAfter = nIndex; nAfter < m_nLength; nAfter++) {
            ptData[nAfter+1] = m_ptData[nAfter];
        }

        delete[] m_ptData;

        m_ptData = ptData;
        m_nLength += 1;
    }

    void Remove(int nIndex) {
        assert(nIndex >= 0 && nIndex < m_nLength);

        T *ptData = new T[m_nLength-1];
        
        for (int nBefore = 0; nBefore < nIndex; nBefore++) {
            ptData[nBefore] = m_ptData[nBefore];
        }

        for (int nAfter = nIndex+1; nAfter < m_nLength; nAfter++) {
            ptData[nAfter-1] = m_ptData[nAfter];
        }

        delete[] m_ptData;

        m_ptData = ptData;
        m_nLength -= 1;
    }

    void InsertAtBeginning(T tValue) { InsertBefore(tValue, 0); }
    void InsertAtEnd(T tValue) { InsertBefore(tValue, m_nLength); }

    int GetLength();
};

template <typename T>
int Array<T>::GetLength() { return m_nLength; }

#endif /* end of include guard: ARRAY_H */
