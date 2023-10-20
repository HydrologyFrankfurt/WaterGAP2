#pragma once

#include "def.h"
#include <string>
#include <type_traits>
#include <chrono>

// forward declaration
template < typename CellType , const bool AlwaysInitialize , int const CellCount >
class Grid;

/**
 * @brief interface for all GridAdapters to come in the future and present...
 * 
 * @tparam CellType the CellType of the grid to be written and read with this adapter
 * @tparam CellCount the CellCount of the grid to be written and read with this adapter
 */
template<typename CellType = double, const int CellCount = ng, const bool init = true>
class IGridAdapter{
	public:		

		/**
		 * @brief writes the grid to file
		 * 
		 * The grid cant be passed as const, as we need to swap the byte order
		 * in place for maximum speed
		 * 
		 * @param grid the grid to write
		 * @param filename the file to write to
		 */
		virtual void writeGrid(Grid<CellType,init,CellCount>& grid, std::string filename) {};

		/**
		 * @brief reads the grid from file
		 * 
		 * @param grid reference to the grid to read in. will be filled by this function.
		 * @param filename the file to write to
		 */
		virtual void readGrid(Grid<CellType,init,CellCount>& grid, std::string filename) {};
};

template<typename CellType = double, const int CellCount = ng, const bool init = true>
class GridUNFAdapter : public IGridAdapter<CellType,CellCount,init>
{
	public:

		/**
		 * @brief writes the grid to file
		 * 
		 * The grid cant be passed as const, as we need to swap the byte order
		 * in place for maximum speed
		 * 
		 * @param grid the grid to write
		 * @param filename the file to write to
		 */
		void writeGrid(Grid<CellType,init,CellCount>& grid, std::string filename) override
		{
			std::ofstream file(filename, std::ios::binary);
	
			// store double grids as float grids
			if(std::is_same<CellType,double>::value)
			{
				// use a static function member so we dont need to allocate memory everytime
				Grid<float,false,CellCount> floatgrid;
				floatgrid.initialize(grid.getCellCount());

				floatgrid = grid;

				floatgrid.doByteSwap();

				file.write((char*)floatgrid.getDataPointer(),sizeof(float)*floatgrid.getCellCount());
			}
			else
			{	
				grid.doByteSwap();

				file.write((char*)grid.getDataPointer(),sizeof(CellType)*grid.getCellCount());

				// swap the bytes again, so the original byte order is restored.
				grid.doByteSwap();
			}

			file.close();
		};

		/**
		 * @brief reads the grid from file
		 * 
		 * @param grid reference to the grid to read in. will be filled by this function.
		 * @param filename the file to write to
		 */
		void readGrid(Grid<CellType,init,CellCount>& grid, std::string filename) override
		{		
			std::ifstream file(filename, std::ios::binary);
            if(!file) {
                std::cout << filename << " not found." << std::endl;
                exit(-1);
            }
			// double grids are stored as float grids!
			if(std::is_same<CellType,double>::value)
			{
				// use a static function member so we dont need to allocate memory everytime
				Grid<float,false,CellCount> floatgrid;
				floatgrid.initialize(grid.getCellCount());

				file.read((char*)floatgrid.getDataPointer(),sizeof(float)*floatgrid.getCellCount());

				floatgrid.doByteSwap();
				
				#pragma omp parallel for
				for(size_t i = 0; i < grid.getCellCount(); i++)
				{
					grid[i] = static_cast<double>(floatgrid[i]);
				}

			}
			else
			{				
				file.read((char*)grid.getDataPointer(),sizeof(CellType)*grid.getCellCount());

				grid.doByteSwap();

			}
			
			file.close();

		};
		
};
