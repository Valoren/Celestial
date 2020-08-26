/*
 * matrix.h
 *
 * Copyright 2019 Miquel Bernat Laporta i Granados
 * <mlaportaigranados@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

#ifndef CELESTIAL_QS_MATRIX_H
#define CELESTIAL_QS_MATRIX_H

#include <vector>

template <typename T> class QSMatrix {
private:
    std::vector<std::vector<T> > mat;
    unsigned rows;
    unsigned cols;

public:
    QSMatrix(unsigned _rows, unsigned _cols, const T& _initial);
    QSMatrix(const QSMatrix<T>& rhs);
    virtual ~QSMatrix();

    // Operator overloading, for "standard" mathematical matrix operations
    QSMatrix<T>& operator=(const QSMatrix<T>& rhs);

    // Matrix mathematical operations
    QSMatrix<T> operator+(const QSMatrix<T>& rhs);
    QSMatrix<T>& operator+=(const QSMatrix<T>& rhs);
    QSMatrix<T> operator-(const QSMatrix<T>& rhs);
    QSMatrix<T>& operator-=(const QSMatrix<T>& rhs);
    QSMatrix<T> operator*(const QSMatrix<T>& rhs);
    QSMatrix<T>& operator*=(const QSMatrix<T>& rhs);
    QSMatrix<T> transpose();

    // Matrix/scalar operations
    QSMatrix<T> operator+(const T& rhs);
    QSMatrix<T> operator-(const T& rhs);
    QSMatrix<T> operator*(const T& rhs);
    QSMatrix<T> operator/(const T& rhs);

    // Matrix/vector operations
    std::vector<T> operator*(const std::vector<T>& rhs);
    std::vector<T> diag_vec();

    // Access the individual elements
    T& operator()(const unsigned& row, const unsigned& col);
    const T& operator()(const unsigned& row, const unsigned& col) const;

// Access the row and column sizes
    unsigned get_rows() const;
    unsigned get_cols() const;

};

#endif