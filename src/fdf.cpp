//
//  fdf.cpp
//  Agilent
//
//  Created by Tobias Wood on 23/08/2013.
//
//

#include "fdf.h"

namespace Agilent {

fdfImage::fdfImage() { }
fdfImage::fdfImage(const string &path, const OpenMode &mode) { open(path, mode); }
void fdfImage::open(const string &path, const OpenMode &mode) {
	// Remove any trailing / to prevent paths being messed up later, check extension
	if (path.back() == '/')
		m_folderPath = path.substr(0, path.size() - 1);
	else
		m_folderPath = path;
	string ext = m_folderPath.substr(m_folderPath.rfind("."));
	if (ext != ".img")
		throw(invalid_argument("Invalid fdf folder extension " + ext + ", must be .img"));
	
	struct dirent *dp;
	DIR *dfd;
	dfd = opendir(m_folderPath.c_str());
	if (dfd == NULL) {
		throw(runtime_error("Could not open fdf folder: " + m_folderPath));
	}
	ifstream pp_file(m_folderPath + "/procpar");
	if (!pp_file)
		throw(runtime_error("Could not open propcar in folder: " + m_folderPath));
	pp_file >> m_pp;
	
	while ((dp = readdir(dfd))) {
		if (dp->d_name[0] == '.') // Ignore ., .., and hidden files
			continue;
		if (strcmp(dp->d_name, "procpar") == 0) {
			continue;
		} else if (strstr(dp->d_name, ".fdf")) { // Ignore any other files that VnmrJ has saved
			// Grab interesting information from first file
			string fname(dp->d_name);
			size_t prefixEnd = fname.find_first_not_of("0123456789");
			string prefix = fname.substr(0, prefixEnd);
			if (m_filePrefix == "") { // Haven't set a prefix yet
				m_filePrefix = prefix;
			} else if (m_filePrefix != prefix) {
				throw(runtime_error("Detected multiple file prefixes in: " + m_folderPath));
			}
			m_files.insert(pair<string, fdfFile>(fname, fdfFile(m_folderPath + "/" + fname)));
		}
	}
	auto f = m_files.begin();
	m_rank = f->second.rank();
	m_dim[0] = f->second.dim(0);
	m_dim[1] = f->second.dim(1);
	if (m_rank == 2) {
		m_dim[2] = static_cast<size_t>(m_pp.realValue("ns"));
	} else {
		m_dim[2] = f->second.dim(2);
	}
	m_dim[3] = m_files.size() / static_cast<size_t>(m_pp.realValue("ns") * m_pp.realValue("ne"));
	m_dim[4] = static_cast<size_t>(m_pp.realValue("ne"));
	closedir(dfd);
}

void fdfImage::close() {
	for (auto &f : m_files)
		f.second.close();
}

const size_t fdfImage::voxelsPerSlice() const { return m_dim[0] * m_dim[1]; }
const size_t fdfImage::voxelsPerVolume() const { return voxelsPerSlice() * m_dim[2]; }
const string fdfImage::filePath(const size_t slice, const size_t image, const size_t echo) const {
	const size_t w = 3;
	stringstream p;
	p << m_folderPath << "/" << m_filePrefix << setw(w) << setfill('0') << slice
	  << "image" << image << "echo" << echo << ".fdf";
	return p.str();
}

} // End namespace Agilent