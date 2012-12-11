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

int main(int argc, const char * argv[])
{
	string path(argv[1]);
	
	Recon::FID_File thisFid(path);
	cout << thisFid.print_header() << endl;
	const complex<double> *one = thisFid.readBlock(0);
	cout << one[0] << endl;
	const complex<double> *two = thisFid.readBlock(22);
	cout << two[0] << endl;
    return 0;
}

