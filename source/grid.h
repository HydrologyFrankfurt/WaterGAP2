#ifndef __grid__h
#define __grid__h

#include "def.h"
#include <array>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include "grid_io_adapters.h"


/*

Usage of the Grid class:

This class should be used whenever you have data over all grid cells.

The Grid<Celltype, AlwaysInitialize, CellCount>-class holds CellCount cells of CellType.
The count default to the number of gridcells, celltype defaults to double.

You can write or read an element of this class using the square brackets, as you would on an array, eg.

Grid GWstorage;
GWstorage[456] = 0.0; 

The AlwaysInitialize-Flag means the needed memory will be allocated as the object is constructed,
if it is false, you need to call the initialize()-Function to allocate the memory.
This is used thoughout WaterGAP if the usage of an Grid is dependent on some option.

If you want to set the number of cells at runtime, not at compile time, just set
AlwaysIntiialize to false, and pass a number of cells to initialize(CellCount). The template
parameter for the number of cells will be ignored.

You can use the VariableChannelGrid class, which has an additional template parameter for the
number of channels. This represents a matrix with size ChannelCount x CellCount.

You can access the elements using brackets, the first element is the cell id and the second the
channel id, e.g.

DailyGrid G_daily;
G_daily(6782, 364) = 0;

Additionally MonthlyGrid, Daily31Grid and DailyGrid are defined as Grids with 12, 31 and 365 channels respectively.

 */


/**
 * @brief A class for storing a Grid array.
 * 
 * This class stores values defined in an array.
 *  
 * @tparam CellType the type for each grid cell 
 * @tparam AlwaysInitialize wether this grid should be initialized on construction 
 * @tparam CellCount the number of grid cells
 */
template < typename CellType  = double,const bool AlwaysInitialize = true, int const CellCount = ng>
class Grid
{
    protected:
        std::vector<CellType> dataArray; 
        bool wasInitialized = false;
        int actualCellCount = CellCount;
    public:          
        Grid()
        {
            if(AlwaysInitialize)
            {
                this->initialize();
            }            
        }

        /**
         * @brief Construct a new Grid object from another Grid object.
         * 
         * The values will be static_casted to this Grid's CellType.
         * 
         * @tparam OtherType the type of the other grid's cells
         * @param other the other grid
         */
        template<typename OtherType,bool init>
        Grid(Grid<OtherType,init,CellCount>& other) {
            this->initialize();

            #pragma omp parallel for
            for(size_t i = 0; i < CellCount; i++)
            {                
                dataArray[i] = static_cast<CellType>(other[i]);                
            }
        }

        /**
         * @brief manually triggers the initialization of this grid instance
         * 
         * a different number of grid cells can be specified, to allow allocation of memory
         * on runtime. the template parameter will be ignored in this case
         * I strongly recommend setting the template parameter to -1 in this case,
         * as there will be an error on missed initialization. The needed constraints
         * on template parameters are available in C++20 and should be added here, such
         * that a negative number of cells is possible if AlwaysInitialize is false. 
         * 
         * @param the number of cells needed, default to template parameter
         */
        virtual void initialize(int actualCellCount = CellCount)
        {          
            dataArray.resize(actualCellCount);
            this->actualCellCount = actualCellCount;
            wasInitialized = true;
        }

        void ensureInitialized()
        {
            if(!wasInitialized)
            {
                initialize();
            }
        }

        /**
         * @brief Get a pointer to the underlying data array. May not be initialized!
         * 
         * @return CellType* pointer to the data of the underlying STL vector
         */
        CellType* getDataPointer()
        {
            return dataArray.data();
        }

        /**
         * @brief access a grid cell
         * 
         * @param cellIndex the index to access
         * @return CellType& a reference to the content of this grid cell
         */
        inline CellType& operator[](const int cellIndex)
        {
            return dataArray[cellIndex];       
        }
        
        /**
         * @brief const access of a grid cell
         * 
         * Allows passing const Grids around. Grids need to be initialized beforehand
         * as there is no way to initialize the grid here...
         * 
         * @throws std::logic_error if the grid is not initialized
         * 
         * @param cellIndex the index to access to
         * @return const CellType& 
         */
        const CellType& operator[](const int cellIndex) const
        {
            return dataArray[cellIndex];            
        }

        /**
         * @brief Read this grid from file using the default/currently selected IO Adapter
         * 
         * @param filename the file to read from
         */
		void read(std::string filename)
		{
            ensureInitialized();
			GridUNFAdapter<CellType,CellCount,AlwaysInitialize> adapter;
			this->read(filename, &adapter);			
		}
    
        /**
         * @brief Read this grid from file using the given IO Adapter
         * 
         * @param filename the file to read from
         * @param io_adapter the adapter to use
         */
        void read(std::string filename, IGridAdapter<CellType,CellCount,AlwaysInitialize>* io_adapter)
        {
            ensureInitialized();
            io_adapter->readGrid(*this, filename);
        }

        /**
         * @brief Writes this grid instance to file using the default/currently selected IO adapter.
         * 
         * @param filename the file to write the grid to
         */
		void write(std::string filename)
		{
            ensureInitialized();
			GridUNFAdapter<CellType,CellCount,AlwaysInitialize> adapter;
			this->write(filename, &adapter);		
		}

        /**
         * @brief Writes this Grid instance (to file) using the given IO Adapter.
         * 
         * @param filename the file to write the grid to
         * @param io_adapter the adapter to use
         */
        void write(std::string filename, IGridAdapter<CellType,CellCount,AlwaysInitialize>* io_adapter)
        {           
            ensureInitialized();
            io_adapter->writeGrid(*this, filename);
        }

        /**
         * @brief Get the number of grid cells in this grid.
         * 
         * This is the value given as the template parameter.
         * 
         * @return int the number of cells in this grid.
         */
        const int getCellCount() const
        {
            return actualCellCount;
        }

        /**
         * @brief sets each cell in the grid to value
         * 
         * @param value the value to set each grid cell to
         */
        void fill(CellType value)
        {            
            ensureInitialized();
            #pragma omp parallel for
            for (size_t i = 0; i < actualCellCount; i++)
            {
                dataArray[i] = value;
            }          
        }

        /**
         * @brief adds each cell value of the right hand side to the cells in this grid
         * 
         * Each value of the grid RHS is added to the corresponding cell in this grid.
         * The resulting grid is also is returned.
         * 
         * @param rhs the grid to add to this grid
         * @return Grid<CellType, CellCount>& the result of the operation
         */
        Grid<CellType, AlwaysInitialize, CellCount>& operator+=(const Grid<CellType, AlwaysInitialize, CellCount>& rhs)
        {
            ensureInitialized();
            #pragma omp parallel for
            for(size_t i = 0; i < actualCellCount; i++)
            {
                dataArray[i] += rhs[i];
            }
            return *this;
        }

        /**
         * @brief sets each cell in the grid to value and return the result
         * 
         * @param value the value to set each grid cell to
         * @return Grid<CellType, CellCount>& the filled grid
         */
        Grid<CellType, AlwaysInitialize, CellCount>& operator=(CellType value)
        {
            ensureInitialized();
            this->fill(value);
            return *this;
        }

        /**
         * @brief sets each cell in the grid to the value of the same gridcell in other
         * 
         * @param other the grid to copy
         * @return Grid<CellType, CellCount>& the filled grid
         */
        template<typename OtherType>
        Grid<CellType, AlwaysInitialize, CellCount> operator=(const Grid<OtherType, AlwaysInitialize, CellCount>& other)
        {
            ensureInitialized();
            if(std::is_same<CellType, OtherType>::value) 
            {
                #pragma omp parallel for
                for(size_t i = 0; i < actualCellCount; i++)
                {                
                    dataArray[i] = other[i];                
                }
            }
            else
            {
                #pragma omp parallel for
                for(size_t i = 0; i < actualCellCount; i++)
                {                
                    dataArray[i] = static_cast<CellType>(other[i]);                
                }
            }
            
            return *this;
        }

        /**
         * @brief Check the grid for negative values
         * 
         * @return true The grid contains at least one negative value
         * @return false The grid contains no values below 0
         */
        bool containsNegative()
        {
            if (!wasInitialized)
            {                
                initialize();
                return false;
            }
            for(size_t i = 0; i < actualCellCount; i++)
            {
                if(dataArray[i] < 0)
                    return true;
            }
            return false;
        }

        /**
         * @brief ensures every grid cell value is below or equal maximumValue.
         * 
         * Every grid cell above maximumValue will be set to maximumValue
         * 
         * @param maximumValue the maximum possible value of the resulting grid
         */
        void ensureBelow(const CellType maximumValue)
        {
            ensureInitialized();
            #pragma omp parallel for
            for(size_t i = 0; i < actualCellCount; i++)
            {
                if(dataArray[i] > maximumValue)
                    dataArray[i] = maximumValue;
            }
        }

        /**
         * @brief ensures every grid cell value is above or equal minimumValue.
         * 
         * Every grid cell below minimumValue will be set to minimumValue
         * 
         * @param minimumValue the minimum possible value of the resulting grid
         */
        void ensureAbove(const CellType minimumValue)
        {
            ensureInitialized();
            #pragma omp parallel for
            for(size_t i = 0; i < actualCellCount; i++)
            {
                if(dataArray[i] < minimumValue)
                    dataArray[i] = minimumValue;
            }
        }

        /**
         * @brief swaps the byte order of each grid cell
         * 
         * This function swaps the byte order of each grid cell value.
         * This is needed for UNF file output. 
         * Depending on the CellType of the Grid, gcc/clangs __builtin_bwap
         * is applied to the value. This can not be compiled on MSVC!
         * 
         */
        void doByteSwap()
        {
            ensureInitialized();
            if(sizeof(CellType) == 2)
            {
                #pragma omp parallel for
                for(size_t i = 0; i < actualCellCount; i++)
                {
                    uint16_t val = * reinterpret_cast<uint16_t*>(&(dataArray[i]));
                    val = __builtin_bswap16(val);
                    dataArray[i] = * reinterpret_cast<CellType*>(&val);                    
                }
                
            }
            else if(sizeof(CellType) == 4)
            {
                #pragma omp parallel for
                for(size_t i = 0; i < actualCellCount; i++)
                {
                    uint32_t val = * reinterpret_cast<uint32_t*>(&(dataArray[i]));
                    val = __builtin_bswap32(val);
                    dataArray[i] = * reinterpret_cast<CellType*>(&val);
                }
            }
            else if(sizeof(CellType) == 8)
            {
                #pragma omp parallel for
                for(size_t i = 0; i < actualCellCount; i++)
                {
                    uint64_t val = * reinterpret_cast<uint64_t*>(&(dataArray[i]));
                    val = __builtin_bswap64(val);
                    dataArray[i] = * reinterpret_cast<CellType*>(&val);
                }
            }
            else
            {
                #pragma omp parallel for
                for(size_t i = 0; i < actualCellCount; i++)
                {
                    char tmp;
                    char* val = (char*) & (dataArray[i]);
                    for(int j = 0; j < sizeof(CellType)/2;j++)
                    {
                        tmp = val[j];
                        val[j] = val[sizeof(CellType)-1-j];
                        val[sizeof(CellType)-1-j] = tmp;                    
                    }    
                }                
            }
        }
};

/**
 * @brief A Grid with multiple channels for each cell.
 * 
 * This represents a 2D-array, in other words a CellCount x ChannelCount matrix.
 * 
 * @tparam ChannelCount
 * @tparam AlwaysInitialize wether to initialize this grid on construction or not
 * @tparam CellType the type of each grid cell, defaults to double
 * @tparam CellCount the number of grid cells, defaults to "ng" defined in "def.h"
 */
template<int const ChannelCount,  typename CellType = double, const bool AlwaysInitialize = true, const int CellCount = ng>
class VariableChannelGrid : public Grid<CellType,AlwaysInitialize,ChannelCount*CellCount>
{
    public:

        VariableChannelGrid() : Grid<CellType,AlwaysInitialize,ChannelCount*CellCount>() {}

        virtual void initialize(int actualCellCount = CellCount)
        {        
            // this does not work?
            // Grid<double, AlwaysInitialize, CellCount>::initialize(actualCellCount*ChannelCount);

            this->dataArray.resize(actualCellCount*ChannelCount);
            this->actualCellCount = actualCellCount*ChannelCount;
            this->wasInitialized = true;
        }

        /**
         * @brief Get the number of channels
         * 
         * This is the template parameter ChannelCount.
         * 
         * @return const int the number of channels of this grid
         */
        const int getChannelCount() const
		{
			return ChannelCount;
		}

        /**
         * @brief Get the number of grid cells
         * 
         * This is the template parameter CellCount
         * 
         * @return const int the number of grid cells in this grid
         */
        const int getCellCount() const
        {
            return CellCount;
        }

        /**
         * @brief Access an element of this grid
         * 
         * @param cell the cell index of the element to be accessed
         * @param channel the channel index of the element to be accessed
         * @return CellType& a reference to the element
         */
        inline CellType& operator()(const int cell, const int channel)
        {
            return this->dataArray[cell*ChannelCount+channel];            
        }

        /**
         * @brief const access an element of this grid
         *  
         * Allows passing const VariableChannelGrids around. Grids need to be initialized beforehand
         * as there is no way to initialize the grid here...
         * 
         * @throws std::logic_error if the grid is not initialized
         * 
         * @param cell the cell index of the element to be accessed
         * @param channel the channel index of the element to be accessed
         * @return CellType& a reference to the element
         */
        const CellType& operator()(const int cell, const int channel) const
        {
            return this->dataArray[cell*ChannelCount+channel];
        }

};

/**
 * @brief A VariableChannelGrid with 12 channels for each cell, one for each month in a year.
 * 
 * In other words, a CellCount x 12 matrix. 
 * 
 * @tparam AlwaysInitialize wether to initialize this grid on construction or not
 * @tparam CellType the type of each grid cell, defaults to double
 * @tparam CellCount the number of grid cells, default to "ng" set in def.h
 */
template<const bool AlwaysInitialize = true, typename CellType = double, const int CellCount = ng> using MonthlyGrid = VariableChannelGrid<12, CellType, AlwaysInitialize, CellCount>;

/**
 * @brief A VariableChannelGrid with 365 channels for each cell, one for each day in a year.
 * 
 * In other words, a CellCount x 365 matrix. 
 * 
 * 
 * @tparam AlwaysInitialize wether to initialize this grid on construction or not
 * @tparam CellType the type of each grid cell, defaults to double
 * @tparam CellCount the number of grid cells, default to "ng" set in def.h
 */
template<const bool AlwaysInitialize = true, typename CellType = double, const int CellCount = ng> using DailyGrid = VariableChannelGrid<365, CellType, AlwaysInitialize, CellCount>;

/**
 * @brief A VariableChannelGrid with 31 channels for each cell, one for each day in a month.
 * 
 * In other words, a CellCount x 31 matrix. 
 * 
 * @tparam AlwaysInitialize wether to initialize this grid on construction or not
 * @tparam CellType the type of each grid cell, defaults to double
 * @tparam CellCount the number of grid cells, default to "ng" set in def.h
 */
template<const bool AlwaysInitialize = true, typename CellType = double, const int CellCount = ng> using Daily31Grid = VariableChannelGrid< 31, CellType, AlwaysInitialize, CellCount>;



#endif