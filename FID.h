//
//  FID.h
//  Nrecon
//
//  Created by Tobias Wood on 11/12/2012.
//
//

#ifndef __Nrecon__FID__
#define __Nrecon__FID__

#include <iostream>

#include <recon_util.h>
#include <FIDFile.h>
#include <procpar.h>

namespace Nrecon {

class FID {

	private:
		FIDFile _fid;
		ParameterList _procpar;
	
	public:
		FID(const string &path);
		~FID();
};

} // End namespace Nrecon

#endif /* defined(__Nrecon__FID__) */
