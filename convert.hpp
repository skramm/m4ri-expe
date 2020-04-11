/**
\file convert.hpp
\brief Holds format conversion functions
*/

#ifndef HG_WRAPPER_M4RI_CONVERT_HPP
#define HG_WRAPPER_M4RI_CONVERT_HPP

#include "wrapper_m4ri.hpp"
#include "binary_mat.hpp"



MatM4ri convertToM4ri( const BinaryMatrix& mat_in )
{
	MatM4ri out( mat_in.nbRows(), mat_in.nbCols() );
	size_t row = 0;
	for( const auto& line: mat_in )
	{
		for( size_t col=0; col<line.size(); col++ )
			out.set( row, col, line[col] );
		row++;
	}
	return out;
}

BinaryMatrix
convertFromM4ri( const MatM4ri& mat_in )
{
	BinaryMatrix out( mat_in.nbRows(), mat_in.nbCols() );

	for( size_t row=0; row<mat_in.nbRows(); row++ )
	{
		BinaryVec& vec = *(out.begin() + row);
		for( size_t col=0; col<out.nbCols(); col++ )
			vec[col] = mzd_read_bit( mat_in._data, row, col );
	}
	return out;
}

#endif // HG_WRAPPER_M4RI_CONVERT_HPP
