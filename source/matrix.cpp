//OE 2018 (based on Cpp Wahlpflichtmodul 2010 by JM Brockmann)
#include "matrix.h"

// ==============================================================
// Standardkonstruktor ohne input, setzt klassenvariablen zu 0
Matrix::Matrix( )
// ==============================================================
{
#ifdef Matrix_DEBUG
   cerr << "DEBUG:: Konstruktor: Matrix::Matrix( ) aufgerufen..." << endl;
#endif
   resize(0, 0);
}

// ==============================================================
// Konstruktor mit Matrizengroesse 
Matrix::Matrix( int r, int c, double val )
// ==============================================================
{
#ifdef Matrix_DEBUG
   cerr << "DEBUG:: Konstruktor: Matrix::Matrix( int r, int c, double val ) aufgerufen..." << endl;
#endif
   resize(r, c);
   fill( begin(), end(), val );
}

// ==============================================================
// Kopierkonstruktor
Matrix::Matrix( const Matrix & other )
// ==============================================================
{
#ifdef Matrix_DEBUG
   cerr << "DEBUG:: Kopierkonstruktor: Matrix::Matrix( const Matrix & other ) aufgerufen..." << endl;
#endif
   // allocate memory
   resize( other.r(), other.c() );
   _rows = other.r();
   _cols = other._cols;
   copy( other.begin(), other.end(), begin() );
}


// ==============================================================
// Destruktor
Matrix::~Matrix( )
// ==============================================================
{
#ifdef Matrix_DEBUG
   cerr << "DEBUG:: Destruktor: Matrix::~Matrix( ) aufgerufen..." << endl;
#endif
}

// ==============================================================
// aendere Dimension der Matrix
void Matrix::resize( int r, int c )
// ==============================================================
{
   _rows = r;
   _cols = c;
   _inc  = 1;
   _ld   = _rows;
   // Alliciere Speicher furr die Daten
   _data.resize(r*c);
}

// ==============================================================
// get element i j
double Matrix::get( int i, int j ) const
// ==============================================================
{
   int idx = i*_inc + j*_ld;
   return( _data.at( idx ) );
}

// ==============================================================
// get element i j
void Matrix::set( int i, int j , double val )
// ==============================================================
{
   int idx = i*_inc + j*_ld;
   _data.at( idx ) = val;
   return;
}

// ==============================================================
// input of a Matrix from keyboard
void Matrix::fromKeyboard( )
// ==============================================================
{
   int rows;
   int cols;
   cerr << "Enter number of rows:";
   cin >> rows;
   cerr << endl << "Enter number of cols:";
   cin >> cols;
   cerr << endl;
   resize( rows, cols );
   double v;
   for( int i = 0; i < r(); i++ )
   {
      for( int j = 0; j < c(); j++ )
      {
         cerr << " Enter element ( " << i << " , " << j << " ) = ";
         cin >> v;
         set(i,j,v);
      }
   }
   return;
}

// ==============================================================
// input of a Matrix from keyboard
void Matrix::disp( )
// ==============================================================
{
   cerr << "[";
   for( int i = 0; i < r(); i++ )
   {
      if ( i > 0 )
          cerr << " ";
      for( int j = 0; j < c(); j++ )
      {
         if( (j == c()-1) )
         {
            cerr << setprecision(5) << setw(7) << setfill( ' ' ) << get(i,j);
         }
         else
         {
            cerr << setprecision(5) << setw(7) << setfill( ' ' ) << get(i,j) << ", ";
         }
      }
      if( i != (r()-1) )
         cerr << endl;
   }
   cerr << " ]_" << size( ) << endl;
   return;
}


// ==============================================================
// formated size
string Matrix::size( ) const
// ==============================================================
{
   stringstream s;
   s << "( " << r() << " x " << c() << " )";
   return( s.str() );
}



// ==============================================================
// Retuen maximal element
double Matrix::max( )
// ==============================================================
{
   return( *max_element( begin(), end() ) );
}

// ==============================================================
// Retuen minimal element
double Matrix::min( )
// ==============================================================
{
   return( *min_element( begin(), end() ) );
}

// ==============================================================
// return transposed copy
Matrix Matrix::transpose( ) const
// ==============================================================
{
   Matrix A;
   A.resize( c(), r() );
   for( int i = 0; i < r(); i++ )
   {
      for( int j = 0; j < c(); j++ )
      {
         A.set(j,i, get(i,j));
      }
   }
   return( A );
}

// ==============================================================
// Kopierkonstruktor
Matrix & Matrix::operator=( const Matrix & other )
// ==============================================================
{
#ifdef Matrix_DEBUG
   cerr << "DEBUG:: Zuweisungsoperator: Matrix & Matrix::operator=( const Matrix & other ) aufgerufen..." << endl;
#endif
   // allocate memory
   resize( other.r(), other.c() );
   _rows = other.r();
   _cols = other._cols;
   copy( other.begin(), other.end(), begin() );
   return *(this);
}

//=========================================================
//Addition
Matrix Matrix::operator+( Matrix &B ) const
//=========================================================
{
   if ( (r()!=B.r()) && (c()!=B.c()))
   {
      cerr<<"Addition unmoeglich, Dimensionen passen nicht."<<endl;
      cerr<<"DIM(A) = "<<  r()<<"x"<<  c()<<" . "<<endl;
      cerr<<"DIM(B) = "<<B.r()<<"x"<<B.c()<<" . "<<endl; exit(-1);
   }
   Matrix Erg = Matrix( r() , c() );
   for(int j = 0 ; j < c() ; j++) 
   {
      for(int i = 0 ; i < r() ; i++)
      {
         Erg.set(i,j, get(i,j)+B.get(i,j) );
      }
   }
   return Erg;
}

//=========================================================
// Subtraktion
Matrix Matrix::operator-( Matrix &B ) const
//=========================================================
{
   if ( (r()!=B.r()) && (c()!=B.c()))
   {
      cerr<<"Addition unmoeglich, Dimensionen passen nicht."<<endl;
      cerr<<"DIM(A) = "<<  r()<<"x"<<  c()<<" . "<<endl;
      cerr<<"DIM(B) = "<<B.r()<<"x"<<B.c()<<" . "<<endl; exit(-1);
   }
   Matrix Erg( r() , c() );
   for(int j = 0 ; j < c() ; j++) 
   {
      for(int i = 0 ; i < r() ; i++)
      {
         Erg(i,j) = operator()( i, j) - B(i,j);
      }
   }
   return Erg;
}


//=========================================================
// Multiplikation
Matrix Matrix::operator*( Matrix &B ) const
//=========================================================
{
   if ( (c()!=B.r()))
   {
      cerr<<"Addition unmoeglich, Dimensionen passen nicht."<<endl;
      cerr<<"DIM(A) = "<<  r()<<"x"<<  c()<<" . "<<endl;
      cerr<<"DIM(B) = "<<B.r()<<"x"<<B.c()<<" . "<<endl; exit(-1);
      exit(-1);
   }
   Matrix Erg( r() , B.c() );
   for(int i = 0 ; i < r() ; i++) 
   {
      for(int j = 0 ; j < c() ; j++)
      {
         for(int k = 0 ; k < B.c() ; k++)
         {
            Erg( i, k ) = this->operator()( i,j ) * B(j,k);
         }
      }
   }
   return Erg;
}

//=========================================================
// Multiplikation mit skalar
Matrix Matrix::operator*( double s ) const
//=========================================================
{
   Matrix erg = *this;
   for(int i = 0 ; i < r()*c() ; i++) 
   {
      erg( i ) *= s;
   }
   return erg;
}

//=========================================================
// ()-operator nur lesend
double Matrix::operator()( int i, int j ) const
//=========================================================
{
   return get(i,j);
}

//=========================================================
// ()-operator lesend/schreibend
double & Matrix::operator()( int i, int j )
//=========================================================
{
	int idx = i*inc() + j*ld();
	return( _data.at( idx ) );
}

//=========================================================
// ()-operator nur lesend
double Matrix::operator()( int i ) const
//=========================================================
{
   return _data.at( i );
}

//=========================================================
// ()-operator lesend/schreibend
double & Matrix::operator()( int i )
//=========================================================
{
	return( _data.at( i ) );
}


//=========================================================
// write binary file
void Matrix::writeBin( const string filename ) const
//=========================================================
{
   ofstream tofile;
   tofile.open( filename.c_str() , ios_base::out | ios_base::binary);
   if (!tofile)
   {
      stringstream error_message;
      error_message << "Open file " << filename << " to write failed...";
      throw runtime_error( error_message.str() );
   }
   int r = _rows;
   int c = _cols;
   double * p = const_cast<double*>(dataPtr());
   tofile.write(reinterpret_cast<char*>(&r), sizeof(int));
   tofile.write(reinterpret_cast<char*>(&c), sizeof(int));
   tofile.write(reinterpret_cast<char*>(p), r*c*sizeof(double));
   tofile.close();
   return;
}

//=============================================================================
//write Matrix data to binary file
void Matrix::readBin(string filename)
//=============================================================================
{ 
   ifstream infile;
   infile.open( filename.c_str() , ios_base::in | ios_base::binary );
   if (!infile)
   {
      stringstream error_message;
      error_message << "Open file " << filename << " to read failed...";
      throw runtime_error( error_message.str() );
   }
	
   infile.read((char *)&_rows, sizeof(int));
   infile.read((char *)&_cols, sizeof(int));
   resize(_rows,_cols);
   infile.read((char *)&_data[0], _rows*_cols*sizeof(double));
   infile.close();
   return;
}

//===========================================================
// sum of all elements
double Matrix::sum() const
//===========================================================
{ 
   double sum = 0.0;
   for( int j = 0; j < _cols ; j++ )
      for( int i = 0; i < _rows; i++ )
         sum+=this->operator()(i,j);
   return(sum);
}

