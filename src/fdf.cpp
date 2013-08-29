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
			size_t prefixEnd = fname.find_first_of("0123456789");
			string prefix = fname.substr(0, prefixEnd);
			if (m_filePrefix == "") { // Haven't set a prefix yet
				m_filePrefix = prefix;
			} else if (m_filePrefix != prefix) {
				throw(runtime_error("Detected multiple file prefixes in: " + m_folderPath));
			}
			auto temp = make_shared<fdfFile>(m_folderPath + "/" + fname);
			m_files.insert(pair<string, shared_ptr<fdfFile>>(fname, temp));
		}
	}
	closedir(dfd);
	
	auto f = m_files.begin();
	size_t ne;
	if (m_pp.contains("ns"))
		m_sls = static_cast<size_t>(m_pp.realValue("ns"));
	else
		m_sls = 1;
	if (m_pp.contains("ne"))
		ne = static_cast<size_t>(m_pp.realValue("ne"));
	else
		ne = 1;
	m_rank = f->second->rank();
	m_dim[0] = f->second->dim(0);
	m_dim[1] = f->second->dim(1);
	if (m_rank == 2) {
		m_dim[2] = m_sls;
	} else {
		m_dim[2] = f->second->dim(2);
	}
	m_dim[3] = m_files.size() / (m_sls * ne);
	m_dim[4] = static_cast<size_t>(m_pp.realValue("ne"));
	m_sl_size = f->second->dataSize();
	// Now we have the joy of calculating a correct orientation field
	// Get "Euler" angles. These describe how to get to the user frame
	// from the magnet frame.
	double psi = m_pp.realValue("psi"), phi = m_pp.realValue("phi"),
	       tht = m_pp.realValue("theta");
	
	double sinphi = sin(phi*M_PI/180.), cosphi = cos(phi*M_PI/180.);
	double sintht = sin(tht*M_PI/180.), costht = cos(tht*M_PI/180.);
	double sinpsi = sin(psi*M_PI/180.), cospsi = cos(psi*M_PI/180.);
	
	// Now for the vox dimensions and offsets in the user frame
	Array3d offset;
	m_voxdim[0] = m_pp.realValue("lro")/m_dim[0];
	offset(0)   = m_pp.realValue("pro") - (m_pp.realValue("lro") - m_voxdim[0])/2.;
	m_voxdim[1] = m_pp.realValue("lpe")/m_dim[1];
	offset(1)   = m_pp.realValue("ppe") - (m_pp.realValue("lpe") - m_voxdim[1])/2.;
	if (m_rank == 2) {
		m_voxdim[2] = m_pp.realValue("thk") + m_pp.realValue("gap");
		offset[2]   = m_pp.realValue("pss", 0); // Find the most negative slice center
		for (size_t i = 1; i < m_sls; i++) {
			if (m_pp.realValue("pss", i) < offset[2])
				offset[2] = m_pp.realValue("pss", i);
		}
	} else {
		m_voxdim[2] = m_pp.realValue("lpe2")/m_dim[2];
		offset[2]   = m_pp.realValue("ppe2") - (m_pp.realValue("lpe2") - m_voxdim[2])/2.;
	}
	// Now build the transform matrix - the 10 is to convert from cm to mm
	Affine3d S; S = Scaling(m_voxdim[0], m_voxdim[1], m_voxdim[2]) * 10.;
	Affine3d T; T = Translation3d(offset[0], offset[1], offset[2]);
	// From Michael Gyngell
	Matrix3d Rd;
	Rd << -cospsi*sinphi + sinpsi*costht*cosphi, -cospsi*cosphi - sinpsi*costht*sinphi, sinpsi*sintht,
	       sinpsi*sinphi + cospsi*costht*cosphi,  sinpsi*cosphi - cospsi*costht*sinphi, cospsi*sintht,
		  -sintht*cosphi, sintht*sinphi, costht;
	Affine3d R(Rd);
	m_transform = (R*T*S).matrix();
}

void fdfImage::close() {
	for (auto &f : m_files)
		f.second->close();
}

const size_t fdfImage::dim(const size_t d) const {
	if (d > 4)
		throw(invalid_argument("Tried to access dimension " + std::to_string(d) + " of image: " + m_folderPath));
	else
		return m_dim[d];
}
const double fdfImage::voxdim(const size_t d) const {
	if (d > 2)
		throw(invalid_argument("Tried to access voxdim " + std::to_string(d) + " of image: " + m_folderPath));
	else
		return m_voxdim[d];
}
const size_t fdfImage::voxelsPerSlice() const { return m_dim[0] * m_dim[1]; }
const size_t fdfImage::voxelsPerVolume() const { return voxelsPerSlice() * m_dim[2]; }
const string fdfImage::filePath(const size_t slice, const size_t image, const size_t echo) const {
	const size_t w = 3;
	stringstream p;
	// fdf Indices are base 1
	p << setfill('0') << m_filePrefix << setw(w) << slice + 1
	  << "image" << setw(w) << image + 1 << "echo" << setw(w) << echo + 1 << ".fdf";
	return p.str();
}
const Matrix4d &fdfImage::ijk_to_xyz() const { return m_transform; }

} // End namespace Agilent