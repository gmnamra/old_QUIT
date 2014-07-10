//
//  Internal.cpp
//  NiftiImage
//
//  Created by Tobias Wood on 18/09/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#include "Internal.h"

namespace Nifti {

/*! Swap size bytes at a time from the given array of n sets of bytes
 *
 *  Declared void * so that the fields from the headers can be passed through
 *  without casting.
 */
void SwapBytes(size_t n, size_t size, void *bytes)
{
	size_t i;
	char *cp0 = (char *)bytes, *cp1, *cp2;
	char swap;
	
	for(i=0; i < n; i++) {
		cp1 = cp0;
		cp2 = cp0 + (size-1);
		while (cp2 > cp1) {
			swap = *cp1; *cp1 = *cp2; *cp2 = swap;
			cp1++; cp2--;
		}
		cp0 += size;
	}
}

/*
 *  Byte swap the individual fields of a NIFTI-1 header struct.
 */
void SwapNiftiHeader(struct nifti_1_header *h)
{
	SwapBytes(1, 4, &h->sizeof_hdr);
	SwapBytes(1, 4, &h->extents);
	SwapBytes(1, 2, &h->session_error);
	
	SwapBytes(8, 2, h->dim);
	SwapBytes(1, 4, &h->intent_p1);
	SwapBytes(1, 4, &h->intent_p2);
	SwapBytes(1, 4, &h->intent_p3);
	
	SwapBytes(1, 2, &h->intent_code);
	SwapBytes(1, 2, &h->datatype);
	SwapBytes(1, 2, &h->bitpix);
	SwapBytes(1, 2, &h->slice_start);
	
	SwapBytes(8, 4, h->pixdim);
	
	SwapBytes(1, 4, &h->vox_offset);
	SwapBytes(1, 4, &h->scl_slope);
	SwapBytes(1, 4, &h->scl_inter);
	SwapBytes(1, 2, &h->slice_end);
	
	SwapBytes(1, 4, &h->cal_max);
	SwapBytes(1, 4, &h->cal_min);
	SwapBytes(1, 4, &h->slice_duration);
	SwapBytes(1, 4, &h->toffset);
	SwapBytes(1, 4, &h->glmax);
	SwapBytes(1, 4, &h->glmin);
	
	SwapBytes(1, 2, &h->qform_code);
	SwapBytes(1, 2, &h->sform_code);
	
	SwapBytes(1, 4, &h->quatern_b);
	SwapBytes(1, 4, &h->quatern_c);
	SwapBytes(1, 4, &h->quatern_d);
	SwapBytes(1, 4, &h->qoffset_x);
	SwapBytes(1, 4, &h->qoffset_y);
	SwapBytes(1, 4, &h->qoffset_z);
	
	SwapBytes(4, 4, h->srow_x);
	SwapBytes(4, 4, h->srow_y);
	SwapBytes(4, 4, h->srow_z);
	
	return ;
}

/*
 *! Byte swap as an ANALYZE 7.5 header
 */
void SwapAnalyzeHeader(nifti_analyze75 * h)
{
	SwapBytes(1, 4, &h->sizeof_hdr);
	SwapBytes(1, 4, &h->extents);
	SwapBytes(1, 2, &h->session_error);
	
	SwapBytes(8, 2, h->dim);
	SwapBytes(1, 2, &h->unused8);
	SwapBytes(1, 2, &h->unused9);
	SwapBytes(1, 2, &h->unused10);
	SwapBytes(1, 2, &h->unused11);
	SwapBytes(1, 2, &h->unused12);
	SwapBytes(1, 2, &h->unused13);
	SwapBytes(1, 2, &h->unused14);
	
	SwapBytes(1, 2, &h->datatype);
	SwapBytes(1, 2, &h->bitpix);
	SwapBytes(1, 2, &h->dim_un0);
	
	SwapBytes(8, 4, h->pixdim);
	
	SwapBytes(1, 4, &h->vox_offset);
	SwapBytes(1, 4, &h->funused1);
	SwapBytes(1, 4, &h->funused2);
	SwapBytes(1, 4, &h->funused3);
	
	SwapBytes(1, 4, &h->cal_max);
	SwapBytes(1, 4, &h->cal_min);
	SwapBytes(1, 4, &h->compressed);
	SwapBytes(1, 4, &h->verified);
	
	SwapBytes(1, 4, &h->glmax);
	SwapBytes(1, 4, &h->glmin);
	
	SwapBytes(1, 4, &h->views);
	SwapBytes(1, 4, &h->vols_added);
	SwapBytes(1, 4, &h->start_field);
	SwapBytes(1, 4, &h->field_skip);
	
	SwapBytes(1, 4, &h->omax);
	SwapBytes(1, 4, &h->omin);
	SwapBytes(1, 4, &h->smax);
	SwapBytes(1, 4, &h->smin);
}

} // End namespace Nifti
