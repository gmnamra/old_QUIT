//
//  FID.h
//  Nrecon
//
//  Created by Tobias Wood on 11/12/2012.
//
//

#ifndef AGILENT_FID
#define AGILENT_FID

#include <string>
#include <iostream>
#include <complex>
#include <map>

#include "recon_util.h"
#include "FIDFile.h"
#include "procpar.h"

namespace Recon {

enum AppTypes {
	im2D = 0,
	im3D
};

class FID {
	private:
		string _bundlePath;
		FIDFile _fid;
		ParameterList _procpar;
		AppTypes _appType;
		
		static const map<string, AppTypes> AppTypeMap;
	public:
		FID(const string &path);
		//~FID();
		
		const string print_info() const;
		complex<double> *readKSpace();
		
		const int nVolumes() const;
		const int nDim0() const;    //!< First data size (usually read-out)
		const int nDim1() const;    //!< Second data size (usually phase-encode)
		const int nDim2() const;    //!< Third data size (usually phase-encode 2 or slices)
};

} // End namespace Nrecon

#endif /* defined(__Nrecon__FID__) */
