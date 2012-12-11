//
//  FID.cpp
//  Nrecon
//
//  Created by Tobias Wood on 11/12/2012.
//
//

#include "FID.h"

namespace Recon {
const map<string, AppTypes> FID::AppTypeMap
{
	{"im2D", im2D},
	{"im3D", im3D}
};

FID::FID(const string &path) {
	
	_bundlePath = path;
	// Remove any trailing directory delimiters
	if (_bundlePath.back() == '/')
		_bundlePath.resize(path.size() - 1);
	std::string ext = _bundlePath.substr(_bundlePath.find_last_of(".") + 1);
	if (ext != "fid")
		NRECON_FAIL("Invalid extension for FID Bundle: " + _bundlePath);
	_fid.open(_bundlePath + "/fid");
	if (!ReadProcpar(_bundlePath + "/procpar", _procpar))
		NRECON_FAIL("Failed to read procpar within FID Bundle: " + _bundlePath);
	if (!ParExists(_procpar, "seqcon") || !ParExists(_procpar, "apptype"))
		NRECON_FAIL("No apptype or seqcon found in FID Bundle: " + _bundlePath);
	_appType = AppTypeMap.find(StringValue(_procpar, "apptype"))->second;
}

const string FID::print_info() const {
	stringstream ss;
	
	ss << "FID Bundle: " << _bundlePath << endl
	   << _fid.print_header() << endl
	   << "Procpar contains " << _procpar.size() << " parameters." << endl;
	return ss.str();
}

const complex<double> *FID::readKSpace() {
	complex<double> *kSpace = new complex<double>[_fid.nComplexPerBlock() * _fid.nBlocks()];
	
	int blockOffset = 0;
	for (int b = 0; b < _fid.nBlocks(); b++) {
		const complex<double> *thisBlock = _fid.readBlock(b);
		for (int k = 0; k < _fid.nComplexPerBlock(); k++)
			kSpace[blockOffset + k] = thisBlock[k];
		blockOffset += _fid.nComplexPerBlock();
		delete[] thisBlock;
	}
	return kSpace;
}

const int FID::nVolumes() const { return 1; }
const int FID::nDim0() const { return static_cast<int>(RealValue(_procpar, "np") / 2); }
const int FID::nDim1() const { return static_cast<int>(RealValue(_procpar, "nv")); }
const int FID::nDim2() const {
	if (_appType == im2D)
		return static_cast<int>(RealValue(_procpar, "ns"));
	else if (_appType == im3D)
		return static_cast<int>(RealValue(_procpar, "nv2"));
	else
		NRECON_FAIL("Unknown application type.");
}

} // End namespace Nrecon
