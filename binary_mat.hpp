/**
\file
\brief holds class BinaryMatrix that models a matrix in GF(2)
*/
#ifndef HG_BINARY_MATRIX_HPP
#define HG_BINARY_MATRIX_HPP

#include <boost/dynamic_bitset.hpp>


#define COUT std::cout

/// Holds a path as a binary vector
/**
Based on https://www.boost.org/doc/libs/release/libs/dynamic_bitset/dynamic_bitset.html

For a graph of \f$n\f$ vertices, its size needs to be \f$ n.(n-1)/2 \f$

Example: for the path 1-3-4 on a graph of 5 vertices (0 - 4), the vector will have a size of 10 elements:

\verbatim
edge:    0  0  0  0  1  1  1  2  2  3
         1  2  3  4  2  3  4  3  4  4
--------------------------------------
vector:  0  0  0  0  0  1  1  0  0  1
\endverbatim
*/
typedef boost::dynamic_bitset<> BinaryVec;


//-------------------------------------------------------------------------------------------
/// Details on BinaryMatrix
struct BinaryMatInfo
{
    size_t nbLines  = 0;
    size_t nbCols   = 0;
    size_t nbOnes   = 0;
    size_t nb0Cols  = 0;   ///< nb of columns with only 0 values
    size_t nb0Lines = 0;  ///< nb of lines with only 0 values

    void print( std::ostream& f ) const
    {
		f << "BinaryMatInfo:"
			<< "\n-nbLines ="  << nbLines
			<< "\n-nbCols ="   << nbCols
			<< "\n-nbOnes ="   << nbOnes
			<< "\n-nb0Lines =" << nb0Lines
			<< "\n-nb0Cols ="  << nb0Cols
			<< '\n';
    }
};

//-------------------------------------------------------------------------------------------
/// A binary matrix, implemented as a vector of BinaryVec
/**
This type will allow to fetch some relevant information on what the matrix holds
*/
struct BinaryMatrix
{
	std::vector<BinaryVec> _data;

	BinaryMatrix( size_t nbLines, size_t nbCols )
	{
		assert( nbLines>0 );
		assert( nbCols>0 );
		_data.resize( nbLines );
		std::for_each( _data.begin(), _data.end(), [nbCols](BinaryVec& bv){ bv.resize(nbCols); } ); // lambda
//		std::cout << __FUNCTION__ << "(): nb cols=" << _data[0].size() << '\n';
	}
	BinaryMatrix( size_t nbLines )
	{
		assert( nbLines>0 );
		_data.resize( nbLines );
	}
	BinaryMatrix()
	{}

	size_t nbLines() const { return _data.size(); }
	size_t nbCols()  const
	{
		if( 0==nbLines() )
			return 0;
		return _data.at(0).size();
	}

	auto begin() -> decltype(_data.begin()) { return _data.begin(); }
	auto end()   -> decltype(_data.end())   { return _data.end();   }
	const auto begin() const -> decltype(_data.begin()) { return _data.begin(); }
	const auto end()   const -> decltype(_data.end())   { return _data.end();   }

	void addLine( const BinaryVec& bvec )
	{
		_data.push_back(bvec);
	}

	void addCol( const BinaryVec& vin )
	{
        assert( vin.size() == nbLines() );
        for( size_t i=0; i<vin.size(); i++ )
			_data[i].push_back( vin[i] );
	}

	BinaryVec getCol( size_t col ) const
	{
		assert( col<nbCols() );
		BinaryVec out( nbLines() );
		for( size_t i=0; i<nbLines(); i++ )
			out[i] = line(i)[col];
		return out;
	}

	const BinaryVec& line( size_t idx ) const
	{
		assert( idx<nbLines() );
		return _data[idx];
	}
	BinaryVec& line( size_t idx )
	{
		assert( idx<nbLines() );
		return _data[idx];
	}
    BinaryMatInfo getInfo() const
    {
		BinaryMatInfo info;
		info.nbLines = nbLines();
		assert( _data.size() );
		info.nbCols = nbCols();
		std::for_each( _data.begin(), _data.end(), [&info](const BinaryVec& v){ info.nbOnes+= v.count();} ); // count 1

		for( size_t i=0; i<nbCols(); i++ )
		{
			bool foundOne=false;
			for( size_t j=0; j<nbLines(); j++ )
			{
				if( _data[j][i] == 1 )
				{
					foundOne = true;
					break;
				}
			}
			if( !foundOne )
				info.nb0Cols++;
		}
		return info;
    }

	std::vector<size_t> getNonEmptyCols() const
	{
		std::vector<size_t> out;
		for( size_t col=0; col<nbCols(); col++ )
		{
			bool foundOne=false;
			for( size_t row=0; row<nbLines(); row++ )
			{
				if( _data[row][col] == 1 )
				{
					foundOne = true;
					break;
				}
			}
			if( foundOne )
				out.push_back(col);
		}
		return out;
	}

	void print( std::ostream& f, std::string msg=std::string() ) const
	{
		size_t i=0;
		f << "BinaryMatrix: " << msg << ", nbLines=" << nbLines() << " nbCols=" << nbCols() << "\n";
		for( auto line: *this )
		{
			f << std::setw(4) << i++ << ": | ";

			for( size_t i=0; i<line.size(); i++ )
			{
				f << line[i];
				if( !((i+1)%4) && i != line.size()-1 )
					f << '.';
			}
			f << " | #" << line.count() << "\n";
		}
	}

/// Returns a vector having as size the number of columns and holding the number of "1" the column has
	std::vector<size_t> getColumnCount() const
	{
		std::vector<size_t> out( nbCols(), 0 );
		for( size_t col=0; col<nbCols(); col++ )
		{
			size_t n=0;
			for( size_t row=0; row<nbLines(); row++ )
			{
				if( _data[row][col] == 1 )
					n++;
			}
			out[col] = n;
		}
		return out;
	}
};
//-------------------------------------------------------------------------------------------
/// Gaussian binary elimination
/**
- Input: a binary matrix
- Output: a reduced matrix

Assumes no identical rows
*/
//template<typename vertex_t> // TEMP
inline
BinaryMatrix
gaussianElim( BinaryMatrix& m_in, size_t& nbIter ) // /* TEMP */, size_t nbVertices, const std::vector<size_t>& nec )
{
	size_t col = 0;
	size_t nb_rows = m_in.nbLines();
	size_t nb_cols = m_in.nbCols();
	assert( nb_rows > 1 );

	BinaryMatrix m_out;

	nbIter = 0;
	bool done = false;

//	m_in.print( std::cout, "m_in INITIAL" );

// TEMP
//	auto rev_map = buildReverseBinaryMap( nbVertices );

	std::vector<bool> tag(nb_rows,false);
	do
	{
		++nbIter;
		COUT << "\n* start iter " << nbIter << ", current col=" << col
			<< " #tagged lines = " << std::count( tag.begin(),tag.end(), true ) << "\n";

		for( size_t row=0; row<nb_rows; row++ )                // search for first row with a 1 in current column
		{
//			COUT << "considering line " << row << "\n";
			if( tag[row] == false && m_in.line(row)[col] == 1 )    // AND not tagged
			{
				COUT << "row: " << row << ": found 1 in col " << col << "\n"; // found pivot
//				printBitVector( std::cout, m_out.line(row) );
				m_out.addLine( m_in.line(row) );
				COUT << "Adding line " << row << " to OUTMAT at line " << m_out.nbLines()-1 << '\n';

//				printBitVector( std::cout, m_in.line(row) );
//				printBitVector( std::cout, m_out.line(m_out.nbLines()-1) );
				tag[row] = true;
				if( row < nb_rows-1 )
				{
//					for( size_t i=0; i<nb_rows; i++ )      // search for all following rows that have a 1 in that column
//						if( i != row )
					for( size_t i=row+1; i<nb_rows; i++ )      // search for all following rows that have a 1 in that column
					{
						if( tag[i] == false )                  // AND that are not tagged.
							if( m_in.line(i)[col] == 1 )            // it there is, we XOR them with initial line
							{
//								std::cout << " -row " << i << " changes:\nwas: "; printBitVector( std::cout, m_in.line(i) );
								m_in.line(i) = m_in.line(i) ^ m_in.line(row);
//								std::cout << "now: "; printBitVector( std::cout, m_in.line(i) );
#if 0
								auto v_pvertex = buildPairSetFromBinaryVec_v2<vertex_t>( m_in.line(i), rev_map, nec );
								if( false == checkVertexPairSet( v_pvertex, false ) )
									COUT << "Invalid vector!\n";
#endif
							}
					}
				}
				COUT << "BREAK loop\n";
				break;
			}
		}
		COUT << "switch to next col\n";
		col++;
		if( col == nb_cols )
		{
			COUT << "All columns done, end\n";
			done = true;
		}
		if( std::find(tag.begin(),tag.end(), false ) == tag.end() )
		{
			COUT << "All lines tagged, end\n";
			done = true;
		}
#if 0
	COUT << "-m_in: " << m_in.nbLines() << "x" << m_in.nbCols() << '\n';
		m_in.print( std::cout,  "m_in" );
	COUT << "-m_out: " << m_out.nbLines() << "x" << m_out.nbCols() << '\n';
		m_out.print( std::cout, "m_out" );
#endif
	}
	while( !done );
	return m_out;
}
//-------------------------------------------------------------------------------------------

#endif // HG_BINARY_MATRIX_HPP
