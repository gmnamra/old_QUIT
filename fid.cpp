//
//  fid.cpp
//  Nrecon
//
//  Created by Tobias Wood on 21/11/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include "fid.h"

/*! Swap siz bytes at a time from the given array of n sets of bytes
 *
 *  Declared void * so that the fields from the headers can be passed through
 *  without casting.
 */
void FID::SwapBytes(size_t n, int size, void *bytes) {
	size_t i;
	char *cp0 = (char *)bytes, *cp1, *cp2;
	char swap;
	
	for(i=0; i < n; i++) {
		cp1 = cp0;
		cp2 = cp0 + (size-1);
		while (cp2 > cp1)
		{
			swap = *cp1; *cp1 = *cp2; *cp2 = swap;
			cp1++; cp2--;
		}
		cp0 += size;
	}
}

FID::SwapFileHeader(FileHeader *hdr) {
	SwapBytes(1, 4, &hdr->nblocks);
	SwapBytes(1, 4, &hdr->ntraces);
	SwapBytes(1, 4, &hdr->np);
	SwapBytes(1, 4, &hdr->ebytes);
	SwapBytes(1, 4, &hdr->tbytes);
	SwapBytes(1, 4, &hdr->bbytes);
	SwapBytes(1, 2, &hdr->vers_id);
	SwapBytes(1, 2, &hdr->status);
	SwapBytes(1, 4, &hdr->nbheaders);
}

FID::SwapBlockHeader(BlockHeader *hdr) {
	SwapBytes(1, 2, &hdr->scale);
	SwapBytes(1, 2, &hdr->status);
	SwapBytes(1, 2, &hdr->index);
	SwapBytes(1, 2, &hdr->mode);
	SwapBytes(1, 4, &hdr->ctcount);
	SwapBytes(1, 4, &hdr->lpval);
	SwapBytes(1, 4, &hdr->rpval);
	SwapBytes(1, 4, &hdr->lvl);
	SwapBytes(1, 4, &hdr->tlt);
}

bool FID::CPUisBigEndian() {
	union {
		unsigned char ch[sizeof(short)];
		short         ss;
	} byteOrder;
	// Little endian will put this in byteOrder.ch[0],
	// big endian will put this in byteOrder.ch[1]
	byteOrder.ss = 1;
	if (byteOrder.ch[1] == 1) // We have a big endian CPU, need to flip bytes in the FID
		return true;
	else
		return false;
}

FID::FID(const string& path) {
	
	FileHeader hdr;
	ifstream FID_File(path, ios::in | ios::binary);
	if (FID_File.read(reinterpret_cast<char *>(&hdr), sizeof(hdr))) {
		_swap = CPUisBigEndian();
		if (_swap)
			SwapFileHeader(&hdr);
		
		_numBlocks = hdr.nblocks;
		_numTraces = hdr.ntraces;
		_numPoints = hdr.np;
		_bytesPerPoint = hdr.ebytes;
		_bytesPerTrace = hdr.tbytes;
		_bytesPerBlock = hdr.bbytes;
		_status = bitset(hdr.status);
		_version_id = bitset(hdr.ver_id);
		_numBlockHeaders = hdr.nbheaders;
}

const string FID::print_info() const {
	stringstream ss;
	
	ss << "Number of blocks: " << _numBlocks << endl
	   << "Number of traces per block: " << _numTraces << endl
	   << "Number of points per trace: " << _numPoints << endl
	   << "Number of bytes per point: " << _bytesPerPoint << endl
	   << "Number of bytes per trace: " << _bytesPerTrace << endl
	   << "Number of bytes per block: " << _bytesPerBlock << endl
	   << "Status bits: " << status << " Version/ID bits: " << _version_id << endl
	   << "Number of block headers per block: " << _numBlockHeaders << endl;
	return ss.str();
}

