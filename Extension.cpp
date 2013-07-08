//
//  Extension.cpp
//  NiftiImage
//
//  Created by Tobias Wood on 08/07/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#include "Extension.h"

namespace Nifti {

/*----------------------------------------------------------------------
 * Read the extensions into the nifti_image struct   08 Dec 2004 [rickr]
 *
 * This function is called just after the header struct is read in, and
 * it is assumed the file pointer has not moved.  The value in remain
 * is assumed to be accurate, reflecting the bytes of space for potential
 * extensions.
 *
 * return the number of extensions read in, or < 0 on error
 *----------------------------------------------------------------------*/
//static int nifti_read_extensions( nifti_image *nim, znzFile fp, int remain )
//{
//	nifti1_extender    extdr;      /* defines extension existence  */
//	nifti1_extension   extn;       /* single extension to process  */
//	nifti1_extension * Elist;      /* list of processed extensions */
//	long               posn, count;
//	
//	if( !nim || znz_isnull(fp) ) {
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"** nifti_read_extensions: bad inputs (%p,%p)\n",
//					(void *)nim, (void *)fp);
//		return -1;
//	}
//	
//	posn = znztell(fp);
//	
//	if( (posn != sizeof(nifti_1_header)) &&
//       (nim->nifti_type != NIFTI_FTYPE_ASCII) )
//		fprintf(stderr,"** WARNING: posn not header size (%ld, %d)\n",
//				posn, (int)sizeof(nifti_1_header));
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"-d nre: posn = %ld, offset = %d, type = %d, remain = %d\n",
//				posn, nim->iname_offset, nim->nifti_type, remain);
//	
//	if( remain < 16 ){
//		if( g_opts.debug > 2 ){
//			if( g_opts.skip_blank_ext )
//				fprintf(stderr,"-d no extender in '%s' is okay, as "
//						"skip_blank_ext is set\n",nim->fname);
//			else
//				fprintf(stderr,"-d remain=%d, no space for extensions\n",remain);
//		}
//		return 0;
//	}
//	
//	count = (int)znzread( extdr.extension, 1, 4, fp ); /* get extender */
//	
//	if( count < 4 ){
//		if( g_opts.debug > 1 )
//			fprintf(stderr,"-d file '%s' is too short for an extender\n",
//					nim->fname);
//		return 0;
//	}
//	
//	if( extdr.extension[0] != 1 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d extender[0] (%d) shows no extensions for '%s'\n",
//					extdr.extension[0], nim->fname);
//		return 0;
//	}
//	
//	remain -= 4;
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"-d found valid 4-byte extender, remain = %d\n", remain);
//	
//	/* so we expect extensions, but have no idea of how many there may be */
//	
//	count = 0;
//	Elist = NULL;
//	while (nifti_read_next_extension(&extn, nim, remain, fp) > 0)
//	{
//		if( nifti_add_exten_to_list(&extn, &Elist, count+1) < 0 ){
//			if( g_opts.debug > 0 )
//				fprintf(stderr,"** failed adding ext %ld to list\n", count);
//			return -1;
//		}
//		
//		/* we have a new extension */
//		if( g_opts.debug > 1 ){
//			fprintf(stderr,"+d found extension #%ld, code = 0x%x, size = %d\n",
//					count, extn.ecode, extn.esize);
//			if( extn.ecode == NIFTI_ECODE_AFNI && g_opts.debug > 2 ) /* ~XML */
//				fprintf(stderr,"   AFNI extension: %.*s\n",
//						extn.esize-8,extn.edata);
//			else if( extn.ecode == NIFTI_ECODE_COMMENT && g_opts.debug > 2 )
//				fprintf(stderr,"   COMMENT extension: %.*s\n",        /* TEXT */
//						extn.esize-8,extn.edata);
//		}
//		remain -= extn.esize;
//		count++;
//	}
//	
//	if( g_opts.debug > 2 ) fprintf(stderr,"+d found %ld extension(s)\n", count);
//	
//	nim->num_ext = (int)count;
//	nim->ext_list = Elist;
//	
//	return (int)count;
//}
//
//
///*----------------------------------------------------------------------*/
///*! nifti_add_extension - add an extension, with a copy of the data
// 
// Add an extension to the nim->ext_list array.
// Fill this extension with a copy of the data, noting the
// length and extension code.
// 
// \param nim    - nifti_image to add extension to
// \param data   - raw extension data
// \param length - length of raw extension data
// \param ecode  - extension code
// 
// \sa extension codes NIFTI_ECODE_* in nifti1_io.h
// \sa nifti_free_extensions, valid_nifti_extensions, nifti_copy_extensions
// 
// \return 0 on success, -1 on error (and free the entire list)
// *//*--------------------------------------------------------------------*/
//int nifti_add_extension(nifti_image *nim, const char * data, int len, int ecode)
//{
//	nifti1_extension ext;
//	
//	/* error are printed in functions */
//	if( nifti_fill_extension(&ext, data, len, ecode) )                 return -1;
//	if( nifti_add_exten_to_list(&ext, &nim->ext_list, nim->num_ext+1)) return -1;
//	
//	nim->num_ext++;  /* success, so increment */
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------*/
///* nifti_add_exten_to_list     - add a new nifti1_extension to the list
// 
// We will append via "malloc, copy and free", because on an error,
// the list will revert to the previous one (sorry realloc(), only
// quality dolphins get to become part of St@rk!st brand tunafish).
// 
// return 0 on success, -1 on error (and free the entire list)
// *//*--------------------------------------------------------------------*/
//static int nifti_add_exten_to_list( nifti1_extension *  new_ext,
//								   nifti1_extension ** list, long new_length )
//{
//	nifti1_extension * tmplist;
//	
//	tmplist = *list;
//	*list = (nifti1_extension *)malloc(new_length * sizeof(nifti1_extension));
//	
//	/* check for failure first */
//	if( ! *list ){
//		fprintf(stderr,"** failed to alloc %ld extension structs (%ld bytes)\n",
//				new_length, new_length*(int)sizeof(nifti1_extension));
//		if( !tmplist ) return -1;  /* no old list to lose */
//		
//		*list = tmplist;  /* reset list to old one */
//		return -1;
//	}
//	
//	/* if an old list exists, copy the pointers and free the list */
//	if( tmplist ){
//		memcpy(*list, tmplist, (new_length-1)*sizeof(nifti1_extension));
//		free(tmplist);
//	}
//	
//	/* for some reason, I just don't like struct copy... */
//	(*list)[new_length-1].esize = new_ext->esize;
//	(*list)[new_length-1].ecode = new_ext->ecode;
//	(*list)[new_length-1].edata = new_ext->edata;
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d allocated and appended extension #%ld to list\n",
//				new_length);
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------*/
///* nifti_fill_extension  - given data and length, fill an extension struct
// 
// Allocate memory for data, copy data, set the size and code.
// 
// return 0 on success, -1 on error (and free the entire list)
// *//*--------------------------------------------------------------------*/
//static int nifti_fill_extension( nifti1_extension *ext, const char * data,
//                                int len, int ecode)
//{
//	int esize;
//	
//	if( !ext || !data || len < 0 ){
//		fprintf(stderr,"** fill_ext: bad params (%p,%p,%d)\n",
//				(void *)ext, data, len);
//		return -1;
//	} else if( ! nifti_is_valid_ecode(ecode) ){
//		fprintf(stderr,"** fill_ext: invalid ecode %d\n", ecode);
//		return -1;
//	}
//	
//	/* compute esize, first : len+8, and take ceiling up to a mult of 16 */
//	esize = len+8;
//	if( esize & 0xf ) esize = (esize + 0xf) & ~0xf;
//	ext->esize = esize;
//	
//	/* allocate esize-8 (maybe more than len), using calloc for fill */
//	ext->edata = (char *)calloc(esize-8, sizeof(char));
//	if( !ext->edata ){
//		fprintf(stderr,"** NFE: failed to alloc %d bytes for extension\n",len);
//		return -1;
//	}
//	
//	memcpy(ext->edata, data, len);  /* copy the data, using len */
//	ext->ecode = ecode;             /* set the ecode */
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d alloc %d bytes for ext len %d, ecode %d, esize %d\n",
//				esize-8, len, ecode, esize);
//	
//	return 0;
//}
//
//
///*----------------------------------------------------------------------
// * nifti_read_next_extension  - read a single extension from the file
// *
// * return (>= 0 is okay):
// *
// *     success      : esize
// *     no extension : 0
// *     error        : -1
// *----------------------------------------------------------------------*/
//static int nifti_read_next_extension( nifti1_extension * nex, nifti_image *nim,
//									 int remain, znzFile fp )
//{
//	int swap = nim->byteorder != nifti_short_order();
//	int count, size, code;
//	
//	/* first clear nex */
//	nex->esize = nex->ecode = 0;
//	nex->edata = NULL;
//	
//	if( remain < 16 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d only %d bytes remain, so no extension\n", remain);
//		return 0;
//	}
//	
//	/* must start with 4-byte size and code */
//	count = (int)znzread( &size, 4, 1, fp );
//	if( count == 1 ) count += (int)znzread( &code, 4, 1, fp );
//	
//	if( count != 2 ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d current extension read failed\n");
//		znzseek(fp, -4*count, SEEK_CUR); /* back up past any read */
//		return 0;                        /* no extension, no error condition */
//	}
//	
//	if( swap ){
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"-d pre-swap exts: code %d, size %d\n", code, size);
//		
//		nifti_swap_4bytes(1, &size);
//		nifti_swap_4bytes(1, &code);
//	}
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"-d potential extension: code %d, size %d\n", code, size);
//	
//	if( !nifti_check_extension(nim, size, code, remain) ){
//		if( znzseek(fp, -8, SEEK_CUR) < 0 ){      /* back up past any read */
//			fprintf(stderr,"** failure to back out of extension read!\n");
//			return -1;
//		}
//		return 0;
//	}
//	
//	/* now get the actual data */
//	nex->esize = size;
//	nex->ecode = code;
//	
//	size -= 8;  /* subtract space for size and code in extension */
//	nex->edata = (char *)malloc(size * sizeof(char));
//	if( !nex->edata ){
//		fprintf(stderr,"** failed to allocate %d bytes for extension\n",size);
//		return -1;
//	}
//	
//	count = (int)znzread(nex->edata, 1, size, fp);
//	if( count < size ){
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"-d read only %d (of %d) bytes for extension\n",
//					count, size);
//		free(nex->edata);
//		nex->edata = NULL;
//		return -1;
//	}
//	
//	/* success! */
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d successfully read extension, code %d, size %d\n",
//				nex->ecode, nex->esize);
//	
//	return nex->esize;
//}
//
//
///*----------------------------------------------------------------------*/
///*! for each extension, check code, size and data pointer
// *//*--------------------------------------------------------------------*/
//int valid_nifti_extensions(const nifti_image * nim)
//{
//	nifti1_extension * ext;
//	int                c, errs;
//	
//	if( nim->num_ext <= 0 || nim->ext_list == NULL ){
//		if( g_opts.debug > 2 ) fprintf(stderr,"-d empty extension list\n");
//		return 0;
//	}
//	
//	/* for each extension, check code, size and data pointer */
//	ext = nim->ext_list;
//	errs = 0;
//	for ( c = 0; c < nim->num_ext; c++ ){
//		if( ! nifti_is_valid_ecode(ext->ecode) ) {
//			if( g_opts.debug > 1 )
//				fprintf(stderr,"-d ext %d, invalid code %d\n", c, ext->ecode);
//			errs++;
//		}
//		
//		if( ext->esize <= 0 ){
//			if( g_opts.debug > 1 )
//				fprintf(stderr,"-d ext %d, bad size = %d\n", c, ext->esize);
//			errs++;
//		} else if( ext->esize & 0xf ){
//			if( g_opts.debug > 1 )
//				fprintf(stderr,"-d ext %d, size %d not multiple of 16\n",
//						c, ext->esize);
//			errs++;
//		}
//		
//		if( ext->edata == NULL ){
//			if( g_opts.debug > 1 ) fprintf(stderr,"-d ext %d, missing data\n", c);
//			errs++;
//		}
//		
//		ext++;
//	}
//	
//	if( errs > 0 ){
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"-d had %d extension errors, none will be written\n",
//					errs);
//		return 0;
//	}
//	
//	/* if we're here, we're good */
//	return 1;
//}
//
//
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
///*--------------------------------------------------------------------------*/
///*! free the nifti extensions
// 
// - If any edata pointer is set in the extension list, free() it.
// - Free ext_list, if it is set.
// - Clear num_ext and ext_list from nim.
// 
// \return 0 on success, -1 on error
// 
// \sa nifti_add_extension, nifti_copy_extensions
// *//*------------------------------------------------------------------------*/
//int nifti_free_extensions( nifti_image *nim )
//{
//	int c ;
//	if( nim == NULL ) return -1;
//	if( nim->num_ext > 0 && nim->ext_list ){
//		for( c = 0; c < nim->num_ext; c++ )
//			if ( nim->ext_list[c].edata ) free(nim->ext_list[c].edata);
//		free(nim->ext_list);
//	}
//	/* or if it is inconsistent, warn the user (if we are not in quiet mode) */
//	else if ( (nim->num_ext > 0 || nim->ext_list != NULL) && (g_opts.debug > 0) )
//		fprintf(stderr,"** warning: nifti extension num/ptr mismatch (%d,%p)\n",
//				nim->num_ext, (void *)nim->ext_list);
//	
//	if( g_opts.debug > 2 )
//		fprintf(stderr,"+d free'd %d extension(s)\n", nim->num_ext);
//	
//	nim->num_ext = 0;
//	nim->ext_list = NULL;
//	
//	return 0;
//}
///* return number of extensions written, or -1 on error */
//static int nifti_write_extensions(znzFile fp, nifti_image *nim)
//{
//	nifti1_extension * list;
//	char               extdr[4] = { 0, 0, 0, 0 };
//	int                c, size, ok = 1;
//	
//	if( znz_isnull(fp) || !nim || nim->num_ext < 0 ){
//		if( g_opts.debug > 0 )
//			fprintf(stderr,"** nifti_write_extensions, bad params\n");
//		return -1;
//	}
//	
//	/* if no extensions and user requests it, skip extender */
//	if( g_opts.skip_blank_ext && (nim->num_ext == 0 || ! nim->ext_list ) ){
//		if( g_opts.debug > 1 )
//			fprintf(stderr,"-d no exts and skip_blank_ext set, "
//					"so skipping 4-byte extender\n");
//		return 0;
//	}
//	
//	/* if invalid extension list, clear num_ext */
//	if( ! valid_nifti_extensions(nim) ) nim->num_ext = 0;
//	
//	/* write out extender block */
//	if( nim->num_ext > 0 ) extdr[0] = 1;
//	if( nifti_write_buffer(fp, extdr, 4) != 4 ){
//		fprintf(stderr,"** failed to write extender\n");
//		return -1;
//	}
//	
//	list = nim->ext_list;
//	for ( c = 0; c < nim->num_ext; c++ ){
//		size = (int)nifti_write_buffer(fp, &list->esize, sizeof(int));
//		ok = (size == (int)sizeof(int));
//		if( ok ){
//			size = (int)nifti_write_buffer(fp, &list->ecode, sizeof(int));
//			ok = (size == (int)sizeof(int));
//		}
//		if( ok ){
//			size = (int)nifti_write_buffer(fp, list->edata, list->esize - 8);
//			ok = (size == list->esize - 8);
//		}
//		
//		if( !ok ){
//			fprintf(stderr,"** failed while writing extension #%d\n",c);
//			return -1;
//		} else if ( g_opts.debug > 2 )
//			fprintf(stderr,"+d wrote extension %d of %d bytes\n", c, size);
//		
//		list++;
//	}
//	
//	if( g_opts.debug > 1 )
//		fprintf(stderr,"+d wrote out %d extension(s)\n", nim->num_ext);
//	
//	return nim->num_ext;
//}
///*----------------------------------------------------------------------*/
///*! \fn int nifti_copy_extensions(nifti_image * nim_dest, nifti_image * nim_src)
// \brief copy the nifti1_extension list from src to dest
// 
// Duplicate the list of nifti1_extensions.  The dest structure must
// be clear of extensions.
// \return 0 on success, -1 on failure
// 
// \sa nifti_add_extension, nifti_free_extensions
// */
//int nifti_copy_extensions(nifti_image * nim_dest, const nifti_image * nim_src)
//{
//	char   * data;
//	size_t   bytes;
//	int      c, size, old_size;
//	
//	if( nim_dest->num_ext > 0 || nim_dest->ext_list != NULL ){
//		fprintf(stderr,"** will not copy extensions over existing ones\n");
//		return -1;
//	}
//	
//	if( g_opts.debug > 1 )
//		fprintf(stderr,"+d duplicating %d extension(s)\n", nim_src->num_ext);
//	
//	if( nim_src->num_ext <= 0 ) return 0;
//	
//	bytes = nim_src->num_ext * sizeof(nifti1_extension);  /* I'm lazy */
//	nim_dest->ext_list = (nifti1_extension *)malloc(bytes);
//	if( !nim_dest->ext_list ){
//		fprintf(stderr,"** failed to allocate %d nifti1_extension structs\n",
//				nim_src->num_ext);
//		return -1;
//	}
//	
//	/* copy the extension data */
//	nim_dest->num_ext = 0;
//	for( c = 0; c < nim_src->num_ext; c++ ){
//		size = old_size = nim_src->ext_list[c].esize;
//		if( size & 0xf ) size = (size + 0xf) & ~0xf; /* make multiple of 16 */
//		if( g_opts.debug > 2 )
//			fprintf(stderr,"+d dup'ing ext #%d of size %d (from size %d)\n",
//					c, size, old_size);
//		/* data length is size-8, as esize includes space for esize and ecode */
//		data = (char *)calloc(size-8,sizeof(char));      /* maybe size > old */
//		if( !data ){
//			fprintf(stderr,"** failed to alloc %d bytes for extention\n", size);
//			if( c == 0 ) { free(nim_dest->ext_list); nim_dest->ext_list = NULL; }
//			/* otherwise, keep what we have (a.o.t. deleting them all) */
//			return -1;
//		}
//		/* finally, fill the new structure */
//		nim_dest->ext_list[c].esize = size;
//		nim_dest->ext_list[c].ecode = nim_src->ext_list[c].ecode;
//		nim_dest->ext_list[c].edata = data;
//		memcpy(data, nim_src->ext_list[c].edata, old_size-8);
//		
//		nim_dest->num_ext++;
//	}
//	
//	return 0;
//}
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