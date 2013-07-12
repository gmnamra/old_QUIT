//
//  Extension.cpp
//  NiftiImage
//
//  Created by Tobias Wood on 08/07/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#include "Extension.h"

namespace Nifti {

///*----------------------------------------------------------------------*/
///*! check whether the extension code is valid
// 
// \return 1 if valid, 0 otherwise
// *//*--------------------------------------------------------------------*/
//int nifti_is_valid_ecode( int ecode )
//{
//	if( ecode < NIFTI_ECODE_IGNORE  ||   /* minimum code number (0) */
//       ecode > NIFTI_MAX_ECODE     ||   /* maximum code number     */
//       ecode & 1 )                      /* cannot be odd           */
//		return 0;
//	
//	return 1;
//}
//
//
///*----------------------------------------------------------------------
// * check for valid size and code, as well as can be done
// *----------------------------------------------------------------------*/
//static int nifti_check_extension(nifti_image *nim, int size, int code, int rem)
//{
//	/* check for bad code before bad size */
//	if( ! nifti_is_valid_ecode(code) ) {
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d invalid extension code %d\n",code);
//		return 0;
//	}
//	
//	if( size < 16 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d ext size %d, no extension\n",size);
//		return 0;
//	}
//	
//	if( size > rem ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d ext size %d, space %d, no extension\n", size, rem);
//		return 0;
//	}
//	
//	if( size & 0xf ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d nifti extension size %d not multiple of 16\n",size);
//		return 0;
//	}
//	
//	if( nim->nifti_type == NIFTI_FTYPE_ASCII && size > LNI_MAX_NIA_EXT_LEN ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d NVE, bad nifti_type 3 size %d\n", size);
//		return 0;
//	}
//	
//	return 1;
//}
//
//
//
///*----------------------------------------------------------------------*/
///*! compute the total size of all extensions
// 
// \return the total of all esize fields
// 
// Note that each esize includes 4 bytes for ecode, 4 bytes for esize,
// and the bytes used for the data.  Each esize also needs to be a
// multiple of 16, so it may be greater than the sum of its 3 parts.
// *//*--------------------------------------------------------------------*/
//int nifti_extension_size(nifti_image *nim)
//{
//	int c, size = 0;
//	
//	if( !nim || nim->num_ext <= 0 ) return 0;
//	
//	if( g_opts.debug > 2 ) fprintf(stderr,"-d ext sizes:");
//	
//	for ( c = 0; c < nim->num_ext; c++ ){
//		size += nim->ext_list[c].esize;
//		if( g_opts.debug > 2 ) fprintf(stderr,"  %d",nim->ext_list[c].esize);
//	}
//	
//	if( g_opts.debug > 2 ) fprintf(stderr," (total = %d)\n",size);
//	
//	return size;
//}




}; // End namespace Nifti