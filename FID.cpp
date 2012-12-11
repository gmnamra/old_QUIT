//
//  FID.cpp
//  Nrecon
//
//  Created by Tobias Wood on 11/12/2012.
//
//

#include "FID.h"

namespace Nrecon {

FID::FID(const string &path) {
	
	// Remove any trailing directory delimiters
	if (path.back() == "/")
		path.erase(path.end() - 1);
	std::string ext = path.substr(fname.find_last_of(".") + 1);
	if (ext != ".fid")
		NRECON_FAIL("Invalid extension for FID Bundle: " + path);
	_fid.open(path + "/fid");
	_procpar = ReadProcpar(path + "/procpar");
}

} // End namespace Nrecon
