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
#include <memory>

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
		size_t m_rank, m_sls, m_sl_size, // Number of slabs or slices, number of voxels per slab or slice
		       m_dim[5]; //x y z t e(echoes)
		double m_voxdim[3];
		map<string, shared_ptr<fdfFile>> m_files;
		Matrix4d m_transform;
		const string filePath(const size_t slice, const size_t image, const size_t echo) const;
		
	public:
		fdfImage();
		fdfImage(const string &path, const OpenMode &m = OpenMode::Read);
		
		void open(const string &path, const OpenMode &m = OpenMode::Read);
		void close();
		
		const size_t dim(const size_t d) const;
		const double voxdim(const size_t d) const;
		const size_t voxelsPerSlice() const;
		const size_t voxelsPerVolume() const;
		const Matrix4d &ijk_to_xyz() const;
		
		template<typename T>
		vector<T> readVolume(const size_t vol, const size_t echo = 0) {
			vector<T> Tbuffer(voxelsPerVolume());
			for (size_t sl = 0; sl < m_sls; sl++) { // Slice or slab
				size_t offset = sl * m_sl_size;
				string name = filePath(sl, vol, echo);
				auto file = m_files.find(name);
				if (file == m_files.end())
					throw(runtime_error("Could not find file: " + name));
				vector<T> Tsl = file->second->readData<T>();
				for (size_t i = 0; i < Tsl.size(); i++)
					Tbuffer[offset + i] = Tsl[i];
			}
			return Tbuffer;
		}
};

} // End namespace Agilent

#endif
