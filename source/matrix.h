//OE 2018 (based on Cpp Wahlpflichtmodul 2010 by JM Brockmann)
#ifndef Matrix_H
#define Matrix_H
#include<iterator>
#include<algorithm>
#include<stdexcept>
#include "common.h"


using namespace std;
// Schluesselwort class leitet die Klasse ein
class Matrix
{
      //================================================================
      // Deklaration friend Funktionen
      //================================================================
      // Ausgabe auf Standardausgabe
      friend ostream & operator<<( ostream & out, const Matrix & A);
      // Schreiben in Datei
      friend void operator<<( string filename, const Matrix & A);
      // Lesen aus Datei
      friend void operator>>( string filename, Matrix & A);
   public: 
      //================================================================
      // Konstruktoren
      //================================================================
      // default Konstruktor
      Matrix( );
      // Konstruktor mit Matrizengroesse (Fordert direkt speicherplatz an)
      Matrix( int r, int c, double val = 0.0 );
      // Kopierkonstruktor
      Matrix( const Matrix & other );
      // Destruktor
      ~Matrix( );
      //================================================================
      // Public Memberfunktion
      //================================================================
      // aendert Dimension der Matrix
      void resize( int r, int c );
      // Zugriff auf die Daten
      double r( )   const { return _rows;} // Nur lesender Zugriff
      double c( )   const { return _cols;} // Nur lesender Zugriff
      double ld( )  const { return _ld;}   // Nur lesender Zugriff
      double inc( ) const { return _inc;}  // Nur lesender Zugriff
      // get and set elements
      double get( int i, int j ) const;
      void set( int i, int j, double val );
      // Matrix input
      void fromKeyboard();
      // Matrix output
      void disp( );
      string size( ) const;
      // Matrix computations
      double max( );
      double min( );
      Matrix transpose( ) const;
      // Iterator auf den Speicherbereich der Matrix
      vector<double>::iterator begin( ) { return _data.begin();}
      vector<double>::iterator end( )   { return _data.end();}
      vector<double>::const_iterator begin( ) const { return _data.begin();}
      vector<double>::const_iterator end( )   const { return _data.end();}
      // lesend und schreibend
      double * dataPtr() { return &_data.at(0); }
      // nur lesend
      const double * dataPtr() const { return &_data.at(0); }
      //================================================================
      // Operatoren
      //================================================================
      // Rechnen
      Matrix operator+( Matrix &B ) const;
      Matrix operator-( Matrix &B ) const;
      Matrix operator*( Matrix &B ) const;
      Matrix operator*( double s ) const;

      // Elementweiser Zugriff
      double   operator()( int i, int j ) const; // lesend
      double & operator()( int i, int j );       // lesend und schreibend
      double   operator()( int i ) const; // lesend
      double & operator()( int i );       // lesend und schreibend
      // Zuweisung
      Matrix & operator=( const Matrix & other );
      void writeBin( const string filename ) const;
      void readBin(string filename);
      double sum( ) const;
   private:
      // Spaltenweise Speicherung
      int _inc;
      int _ld;
      int _rows;
      int _cols;
      // Daten der Matrix in Spaltenweiser Speicherung
      vector<double> _data;
}; // Semikolen nicht vergessen !!!

#endif // Matrix_H
