//
//  FID.cpp
//  Nrecon
//
//  Created by Tobias Wood on 11/12/2012.
//
//

#include "FID.h"

namespace Recon {

FID::FID(const string &cpath) {
	
	string path(cpath);
	// Remove any trailing directory delimiters
	if (path.back() == '/')
		path.resize(path.size() - 1);
	std::string ext = path.substr(path.find_last_of(".") + 1);
	if (ext != ".fid")
		NRECON_FAIL("Invalid extension for FID Bundle: " + path);
	_fid.open(path + "/fid");
	if (!ReadProcpar(path + "/procpar", _procpar))
		NRECON_FAIL("Failed to read procpar within FID Bundle: " + path);
}

} // End namespace Nrecon
