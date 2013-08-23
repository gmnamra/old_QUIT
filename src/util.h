//
//  recon_util.h
//  Nrecon
//
//  Created by Tobias Wood on 11/12/2012.
//
//

#ifndef AGILENT_AGILENT
#define AGILENT_AGILENT

namespace Agilent {

#define NRECON_ERROR( err ) do { std::cerr << __PRETTY_FUNCTION__ << ": " << ( err ) << std::flush << std::endl; } while(0)
#define NRECON_FAIL( err ) do { NRECON_ERROR( err ); exit(EXIT_FAILURE); } while(0)

enum Endianness {
	LittleEndian = 0,
	BigEndian
};

Endianness HostEndianness();

template <typename T>
void SwapEndianness(T *ptr, int n = 1) {
	for (int i = 0; i < n; i++) {
		char swap;
		char *lo = reinterpret_cast<char *>(ptr);
		char *hi = lo + (sizeof(T) - 1);
		while (hi > lo) {
			swap = *lo; *lo = *hi; *hi = swap;
			lo++; hi--;
		}
		ptr++;
	}
}

} // End namespace Agilent
#endif
