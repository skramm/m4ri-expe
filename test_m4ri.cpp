/**
\file
\brief test of m4ri
*/

#include <iostream>
#include "wrapper_m4ri.hpp"
#include "binary_mat.hpp"
#include "convert.hpp"

//-------------------------------------------------------------------
int main(int, char*[])
{
	MatM4ri m(10,20);
	m.randomize();
	MatM4ri m2(m);
	std::cout << m;

	std::cout << "mzd_echelonize_naive( m._data, 0 );\n";
	mzd_echelonize_naive( m._data, 0 );
	std::cout << m;
	m = m2;

	std::cout << "mzd_echelonize_naive( m._data, 1 );\n";
	mzd_echelonize_naive( m._data, 1 );
	std::cout << m;

	BinaryMatrix bmat1 = convertFromM4ri( m2 );
	size_t iter=0;
	bmat1.print( std::cout, "bmat1" );
	auto bmat2 = gaussianElim( bmat1, iter );
	bmat2.print( std::cout, "bmat2" );

	MatM4ri m3 = convertToM4ri( bmat2 );
	std::cout << "m3:\n" << m3;
}

//-------------------------------------------------------------------

