//
//  main.cpp
//  Nrecon
//
//  Created by Tobias Wood on 21/11/2012.
//
//

#include <string>
#include <iostream>
using namespace std;
#include "fid.h"
#include "NiftiImage.h"

int main(int argc, const char * argv[])
{
	string path(argv[1]);
	
	Recon::FID thisFid(path);
	cout << thisFid.print_info() << endl;
	const complex<double> *kSpace = thisFid.readKSpace();
	cout << kSpace[0] << endl;
	cout << thisFid.nDim0() << " " << thisFid.nDim1() << " " << thisFid.nDim2() << endl;
	
	NiftiImage output;
	output.setDims(thisFid.nDim0(), thisFid.nDim1(), thisFid.nDim2(), thisFid.nVolumes());
	output.setDatatype(NIFTI_TYPE_FLOAT32);
	output.open("output.nii.gz", NiftiImage::NIFTI_WRITE);
	output.writeAllVolumes(kSpace);
	output.close();
	
    return 0;
}

