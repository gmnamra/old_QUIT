//
//  fdf.h
//  Agilent
//
//  Created by Tobias Wood on 23/08/2013.
//
//

#ifndef AGILENT_FDF
#define AGILENT_FDF

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>

#include <dirent.h>

#include "Eigen/Geometry"

#include "fdfFile.h"
#include "procpar.h"

using namespace std;
using namespace Eigen;

namespace Agilent {

class fdfImage {
	public:
		enum class OpenMode : size_t {
			Read
		};
		
	protected:
		string m_folderPath, m_filePrefix, m_dtype;
		ProcPar m_pp;
		size_t m_rank, m_dim[5]; //x y z(slices) t e(echoes)
		map<string, fdfFile> m_files;
		Matrix4f orientation;
		const string filePath(const size_t slice, const size_t image, const size_t echo) const;
		
	public:
		fdfImage();
		fdfImage(const string &path, const OpenMode &m = OpenMode::Read);
		
		void open(const string &path, const OpenMode &m = OpenMode::Read);
		void close();
		
		const size_t voxelsPerSlice() const;
		const size_t voxelsPerVolume() const;
		
		template<typename T>
		vector<T> readVolume(const size_t vol, const size_t echo = 0) {
			vector<T> Tbuffer(voxelsPerVolume());
			for (size_t slice = 0; slice < m_dim[2]; slice++) {
				size_t offset = slice * voxelsPerSlice();
				vector<T> Tslice = m_files.find(filePath(slice, vol, echo))->second.readData<T>();
				for (size_t i = 0; i < voxelsPerSlice(); i++)
					Tbuffer[offset + i] = Tslice[i];
			}
			return Tbuffer;
		}
};

} // End namespace Agilent

#endif
