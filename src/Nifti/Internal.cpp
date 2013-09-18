//
//  Internal.cpp
//  NiftiImage
//
//  Created by Tobias Wood on 18/09/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#include "Internal.h"

/*! Swap size bytes at a time from the given array of n sets of bytes
 *
 *  Declared void * so that the fields from the headers can be passed through
 *  without casting.
 */
void swapBytes(size_t n, size_t size, void *bytes)
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
void swapNiftiHeader(struct nifti_1_header *h)
{
	swapBytes(1, 4, &h->sizeof_hdr);
	swapBytes(1, 4, &h->extents);
	swapBytes(1, 2, &h->session_error);
	
	swapBytes(8, 2, h->dim);
	swapBytes(1, 4, &h->intent_p1);
	swapBytes(1, 4, &h->intent_p2);
	swapBytes(1, 4, &h->intent_p3);
	
	swapBytes(1, 2, &h->intent_code);
	swapBytes(1, 2, &h->datatype);
	swapBytes(1, 2, &h->bitpix);
	swapBytes(1, 2, &h->slice_start);
	
	swapBytes(8, 4, h->pixdim);
	
	swapBytes(1, 4, &h->vox_offset);
	swapBytes(1, 4, &h->scl_slope);
	swapBytes(1, 4, &h->scl_inter);
	swapBytes(1, 2, &h->slice_end);
	
	swapBytes(1, 4, &h->cal_max);
	swapBytes(1, 4, &h->cal_min);
	swapBytes(1, 4, &h->slice_duration);
	swapBytes(1, 4, &h->toffset);
	swapBytes(1, 4, &h->glmax);
	swapBytes(1, 4, &h->glmin);
	
	swapBytes(1, 2, &h->qform_code);
	swapBytes(1, 2, &h->sform_code);
	
	swapBytes(1, 4, &h->quatern_b);
	swapBytes(1, 4, &h->quatern_c);
	swapBytes(1, 4, &h->quatern_d);
	swapBytes(1, 4, &h->qoffset_x);
	swapBytes(1, 4, &h->qoffset_y);
	swapBytes(1, 4, &h->qoffset_z);
	
	swapBytes(4, 4, h->srow_x);
	swapBytes(4, 4, h->srow_y);
	swapBytes(4, 4, h->srow_z);
	
	return ;
}

/*
 *! Byte swap as an ANALYZE 7.5 header
 */
void swapAnalyzeHeader(nifti_analyze75 * h)
{
	swapBytes(1, 4, &h->sizeof_hdr);
	swapBytes(1, 4, &h->extents);
	swapBytes(1, 2, &h->session_error);
	
	swapBytes(8, 2, h->dim);
	swapBytes(1, 2, &h->unused8);
	swapBytes(1, 2, &h->unused9);
	swapBytes(1, 2, &h->unused10);
	swapBytes(1, 2, &h->unused11);
	swapBytes(1, 2, &h->unused12);
	swapBytes(1, 2, &h->unused13);
	swapBytes(1, 2, &h->unused14);
	
	swapBytes(1, 2, &h->datatype);
	swapBytes(1, 2, &h->bitpix);
	swapBytes(1, 2, &h->dim_un0);
	
	swapBytes(8, 4, h->pixdim);
	
	swapBytes(1, 4, &h->vox_offset);
	swapBytes(1, 4, &h->funused1);
	swapBytes(1, 4, &h->funused2);
	swapBytes(1, 4, &h->funused3);
	
	swapBytes(1, 4, &h->cal_max);
	swapBytes(1, 4, &h->cal_min);
	swapBytes(1, 4, &h->compressed);
	swapBytes(1, 4, &h->verified);
	
	swapBytes(1, 4, &h->glmax);
	swapBytes(1, 4, &h->glmin);
	
	swapBytes(1, 4, &h->views);
	swapBytes(1, 4, &h->vols_added);
	swapBytes(1, 4, &h->start_field);
	swapBytes(1, 4, &h->field_skip);
	
	swapBytes(1, 4, &h->omax);
	swapBytes(1, 4, &h->omin);
	swapBytes(1, 4, &h->smax);
	swapBytes(1, 4, &h->smin);
}
