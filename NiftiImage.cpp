#include "NiftiImage.h"

/*******************************
#pragma mark Methods for ZipFile
*******************************/
bool ZipFile::open(const string &path, const string &mode, const bool zip) {
	_zip = zip;
	if (_zip) {
		_zipped = gzopen(path.c_str(), mode.c_str());
	} else {
		_unzipped = fopen(path.c_str(), mode.c_str());
	}
	if (!(_zipped || _unzipped)) {
		NIFTI_ERROR("Failed to open file: " + path);
		return false;
	}
	return true;
}

void ZipFile::close() {
	if (_zip)
		gzclose(_zipped);
	else
		fclose(_unzipped);
	_zipped = _unzipped = NULL;
}

int ZipFile::read(void *buff, int size)
{
	if (_zip) {
		int remaining = size, totalRead = 0, nread, chunkSize;
		char *cbuff = (char *)buff;
		while (remaining > 0) {
			chunkSize = (remaining < MaxZippedBytes) ? remaining : MaxZippedBytes;
			nread = gzread(_zipped, cbuff, chunkSize);
			remaining -= nread;
			if (nread < chunkSize)
				NIFTI_ERROR("Zipped read short by " << remaining << " bytes.");
			cbuff += nread;
			totalRead += nread;
		}
		return totalRead;
	} else {
		return (int)fread(buff, size, 1, _unzipped) * size;
	}
}

int ZipFile::write(const void *buff, int size)
{
	if (_zip) {
		int remaining = size, chunkSize, totalWritten = 0, nwritten;
		char *chunk = (char *)buff;
		while(remaining > 0 ) {
			chunkSize = (remaining < MaxZippedBytes) ? remaining : MaxZippedBytes;
			nwritten = gzwrite(_zipped, chunk, chunkSize);
			if (nwritten == 0) {
				NIFTI_ERROR("Error during write to file.");
				return 0;
			}
			remaining -= nwritten;
			if(nwritten < (int)chunkSize) {
				NIFTI_ERROR("Zipped write short by " << (chunkSize - nwritten) << " bytes.");
				return totalWritten;
			}
			chunk += nwritten;
			totalWritten += nwritten;
		}
		return totalWritten;
	} else {
		return (int)fwrite(buff, size, 1, _unzipped) * size;
	}
}

long ZipFile::seek(long offset, int whence) {
	long error;
	if (_zip)
		error = gzseek(_zipped, offset, whence);
	else
		error = fseek(_unzipped, offset, whence);
	return error;
}

long ZipFile::tell() {
	if (_zip)
		return gztell(_zipped);
	else
		return ftell(_unzipped);
}

void ZipFile::flush() {
	if (_zip)
		gzflush(_zipped, Z_FINISH);
	else
		fflush(_unzipped);
}

/**********************************
#pragma mark Methods for NiftiImage
**********************************/
/*
 * Map for string representations of NIfTI unit codes.
 *
 *\sa NIFTI1_UNITS group in nifti1.h
 */
const string &NiftiImage::spaceUnits() const
{
	static const StringMap Units
	{
		{ NIFTI_UNITS_METER,  "m" },
		{ NIFTI_UNITS_MM,     "mm" },
		{ NIFTI_UNITS_MICRON, "um" }
	};
	
	static string unknown("Unknown space units code");
	StringMap::const_iterator it = Units.find(xyz_units);
	if (it == Units.end())
		return unknown;
	else
		return it->second;
}
const string &NiftiImage::timeUnits() const
{
	static const StringMap Units
	{
		{ NIFTI_UNITS_SEC,    "s" },
		{ NIFTI_UNITS_MSEC,   "ms" },
		{ NIFTI_UNITS_USEC,   "us" },
		{ NIFTI_UNITS_HZ,     "Hz" },
		{ NIFTI_UNITS_PPM,    "ppm" },
		{ NIFTI_UNITS_RADS,   "rad/s" }
	};
	static string unknown("Unknown time units code");
	StringMap::const_iterator it = Units.find(time_units);
	if (it == Units.end())
		return unknown;
	else
		return it->second;
}

/*
 * Map for string representations of NIfTI intent types.
 *
 *\sa NIFTI1_INTENT_CODES group in nifti1.h
 */
const string &NiftiImage::intentName() const
{
	const StringMap Intents
	{
		{ NIFTI_INTENT_CORREL,     "Correlation statistic" },
		{ NIFTI_INTENT_TTEST,      "T-statistic" },
		{ NIFTI_INTENT_FTEST,      "F-statistic" },
		{ NIFTI_INTENT_ZSCORE,     "Z-score"     },
		{ NIFTI_INTENT_CHISQ,      "Chi-squared distribution" },
		{ NIFTI_INTENT_BETA,       "Beta distribution" },
		{ NIFTI_INTENT_BINOM,      "Binomial distribution" },
		{ NIFTI_INTENT_GAMMA,      "Gamma distribution" },
		{ NIFTI_INTENT_POISSON,    "Poisson distribution" },
		{ NIFTI_INTENT_NORMAL,     "Normal distribution" },
		{ NIFTI_INTENT_FTEST_NONC, "F-statistic noncentral" },
		{ NIFTI_INTENT_CHISQ_NONC, "Chi-squared noncentral" },
		{ NIFTI_INTENT_LOGISTIC,   "Logistic distribution" },
		{ NIFTI_INTENT_LAPLACE,    "Laplace distribution" },
		{ NIFTI_INTENT_UNIFORM,    "Uniform distribition" },
		{ NIFTI_INTENT_TTEST_NONC, "T-statistic noncentral" },
		{ NIFTI_INTENT_WEIBULL,    "Weibull distribution" },
		{ NIFTI_INTENT_CHI,        "Chi distribution" },
		{ NIFTI_INTENT_INVGAUSS,   "Inverse Gaussian distribution" },
		{ NIFTI_INTENT_EXTVAL,     "Extreme Value distribution" },
		{ NIFTI_INTENT_PVAL,       "P-value" },
				
		{ NIFTI_INTENT_LOGPVAL,    "Log P-value" },
		{ NIFTI_INTENT_LOG10PVAL,  "Log10 P-value" },
				
		{ NIFTI_INTENT_ESTIMATE,   "Estimate" },
		{ NIFTI_INTENT_LABEL,      "Label index" },
		{ NIFTI_INTENT_NEURONAME,  "NeuroNames index" },
		{ NIFTI_INTENT_GENMATRIX,  "General matrix" },
		{ NIFTI_INTENT_SYMMATRIX,  "Symmetric matrix" },
		{ NIFTI_INTENT_DISPVECT,   "Displacement vector" },
		{ NIFTI_INTENT_VECTOR,     "Vector" },
		{ NIFTI_INTENT_POINTSET,   "Pointset" },
		{ NIFTI_INTENT_TRIANGLE,   "Triangle" },
		{ NIFTI_INTENT_QUATERNION, "Quaternion" },
				
		{ NIFTI_INTENT_DIMLESS,    "Dimensionless number" }
	};
	static string unknown("Unknown intent code");
	StringMap::const_iterator it = Intents.find(intent_code);
	if (it == Intents.end())
		return unknown;
	else
		return it->second;
}

/*
 * Map for string representations of NIfTI transform codes.
 *
 *\sa NIFTI1_XFORM_CODES group in nifti1.h
 */
const string &NiftiImage::transformName() const
{
	static const StringMap Transforms
	{
		{ NIFTI_XFORM_SCANNER_ANAT, "Scanner Anat" },
		{ NIFTI_XFORM_ALIGNED_ANAT, "Aligned Anat" },
		{ NIFTI_XFORM_TALAIRACH,    "Talairach" },
		{ NIFTI_XFORM_MNI_152,      "MNI_152" }
	};
	static string unknown("Unknown transform code");
	StringMap::const_iterator it = Transforms.find(sform_code);
	if (it == Transforms.end())
		return unknown;
	else
		return it->second;
}

/*
 * Map for string representations of NIfTI slice_codes
 *
 *\sa NIFTI1_SLICE_ORDER group in nifti1.h
 */
const string &NiftiImage::sliceName() const
{
	static const StringMap SliceOrders
	{
		{ NIFTI_SLICE_SEQ_INC,  "sequential_increasing"    },
		{ NIFTI_SLICE_SEQ_DEC,  "sequential_decreasing"    },
		{ NIFTI_SLICE_ALT_INC,  "alternating_increasing"   },
		{ NIFTI_SLICE_ALT_DEC,  "alternating_decreasing"   },
		{ NIFTI_SLICE_ALT_INC2, "alternating_increasing_2" },
		{ NIFTI_SLICE_ALT_DEC2, "alternating_decreasing_2" }
	};
	static string unknown("Unknown slice order code");
	StringMap::const_iterator it = SliceOrders.find(slice_code);
	if (it == SliceOrders.end())
		return unknown;
	else
		return it->second;
}

/*! Swap size bytes at a time from the given array of n sets of bytes
 *
 *  Declared void * so that the fields from the headers can be passed through
 *  without casting.
 */
void NiftiImage::SwapBytes(size_t n, int size, void *bytes)
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
void NiftiImage::SwapNiftiHeader(struct nifti_1_header *h)
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
void NiftiImage::SwapAnalyzeHeader(nifti_analyze75 * h)
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

NiftiImage::~NiftiImage()
{
	if (_mode != CLOSED)
		close();
}

NiftiImage::NiftiImage() :
	_mode(CLOSED), _gz(false), _nii(false), _swap(false), _voxoffset(0),
	_dim(Matrix<int, 7, 1>::Ones()), _voxdim(Matrix<float, 7, 1>::Ones()),
	scaling_slope(1.), scaling_inter(0.), calibration_min(0.), calibration_max(0.),
	freq_dim(0), phase_dim(0), slice_dim(0),
	slice_code(0), slice_start(0), slice_end(0), slice_duration(0),
	toffset(0), xyz_units(NIFTI_UNITS_MM), time_units(NIFTI_UNITS_SEC),
	intent_code(NIFTI_INTENT_DIMLESS), intent_p1(0), intent_p2(0), intent_p3(0),
	intent_name(""), description(""), aux_file(""),
	qform_code(NIFTI_XFORM_SCANNER_ANAT), sform_code(NIFTI_XFORM_SCANNER_ANAT)
{
	_qform.setIdentity(); _sform.setIdentity();
	setDatatype(NIFTI_TYPE_FLOAT32);
}

NiftiImage::NiftiImage(const string &filename, const char &mode) :
	NiftiImage()
{
	if (!open(filename, mode)) {
		NIFTI_FAIL("Could not open file in constructor.");
	}
}

NiftiImage::NiftiImage(const int nx, const int ny, const int nz, const int nt,
		               const float dx, const float dy, const float dz, const float dt,
		               const int datatype) :
	NiftiImage()
{
	setDatatype(datatype);
	_qform.setIdentity(); _sform.setIdentity();
	_dim[0] = nx < 1 ? 1 : nx;
	_dim[1] = ny < 1 ? 1 : ny;
	_dim[2] = nz < 1 ? 1 : nz;
	_dim[3] = nt < 1 ? 1 : nt;
	_voxdim[0] = dx; _voxdim[1] = ny; _voxdim[2] = dz; _voxdim[3] = dt;
}

NiftiImage::NiftiImage(const ArrayXi &dim, const ArrayXf &voxdim, const int &datatype,
                       const Matrix4f &qform, const Matrix4f &sform) :
	NiftiImage()
{
	assert(dim.rows() < 8);
	assert(dim.rows() == voxdim.rows());
	
	_dim.head(dim.rows()) = dim;
	_voxdim.head(voxdim.rows()) = voxdim;
	_qform = qform;
	_sform = sform;
	setDatatype(datatype);
}

NiftiImage &NiftiImage::operator=(const NiftiImage &other)
{
	if (this == &other)
		return *this;
	else if (_mode != CLOSED)
		close();
	
	_dim = other._dim;
	_voxdim = other._voxdim;
	_qform = other._qform;
	_sform = other._sform;
	_basepath = other._basepath;
	_gz = other._gz;
	_nii = other._nii;
	_mode = CLOSED;
	_voxoffset = 0;
	setDatatype(other.datatype());
	scaling_slope = other.scaling_slope;
	scaling_inter = other.scaling_inter;
	calibration_min = other.calibration_min;
	calibration_max = other.calibration_max;
	qform_code = other.qform_code;
	sform_code = other.sform_code;
	freq_dim = other.freq_dim;
	phase_dim = other.phase_dim;
	slice_dim = other.slice_dim;
	slice_code = other.slice_code;
	slice_start = other.slice_start;
	slice_end = other.slice_end;
	slice_duration = other.slice_duration;
	toffset = other.toffset;
	xyz_units = other.xyz_units;
	time_units = other.time_units;
	intent_code = other.intent_code;
	intent_p1 = other.intent_p1;
	intent_p2 = other.intent_p2;
	intent_p3 = other.intent_p3;
	intent_name = other.intent_name;
	description = other.description;
	aux_file = other.aux_file;
	return *this;
}

bool isGZippedFile(const string &fname)
{
	if (fname.find_last_of(".") != string::npos) {
		string ext = fname.substr(fname.find_last_of(".") + 1);
		transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
		if (ext == "gz")
			return true;
	}
	return false;
}

const string NiftiImage::basePath() const {
	return _basepath;
}

const string NiftiImage::imagePath() const {
	string path(_basepath);
	if (_nii) {
		path += ".nii";
	} else {
		path += ".img";
	}
	if (_gz)
		path += ".gz";
	
	return path;
}

const string NiftiImage::headerPath() const {
	string path(_basepath);
	if (_nii) {
		path += ".nii";
	} else {
		path += ".hdr";
	}
	if (_gz)
		path += ".gz";
	
	return path;
}

inline float NiftiImage::fixFloat(const float f)
{
	if (isfinite(f))
		return f;
	else
		return 0.;
}

bool NiftiImage::readHeader() {
	struct nifti_1_header nhdr;
	
	if (_file.read(&nhdr, sizeof(nhdr)) < sizeof(nhdr)) {
		NIFTI_ERROR("Could not read header structure from " + headerPath());
		return false;
	}
	
	// Check if disk and CPU byte order match.
	// The sizeof_hdr field should always be 352, as this is the size of a
	// NIfTI-1 header
	if (nhdr.sizeof_hdr != sizeof(nhdr)) {
		SwapBytes(1, 4, &nhdr.sizeof_hdr);
		if (nhdr.sizeof_hdr != sizeof(nhdr)) {
			NIFTI_ERROR("Could not determine byte order of header " + headerPath());
			return false;
		}
		// If we didn't fail, then we need to swap the header (first swap sizeof back)
		_swap = true;
		SwapBytes(1, 4, &nhdr.sizeof_hdr);
	}
	
	// Check the magic string is set to one of the possible NIfTI values,
	// otherwise process as an ANALYZE file
	int is_nifti = ((nhdr.magic[0]=='n' && nhdr.magic[3]=='\0') &&
                    (nhdr.magic[1]=='i' || nhdr.magic[1]=='+') &&
                    (nhdr.magic[2]>='1' && nhdr.magic[2]<='9')) ? true : false;
	
	if (_swap && is_nifti)
		SwapNiftiHeader(&nhdr);
	else if (_swap)
		SwapAnalyzeHeader((nifti_analyze75 *)&nhdr);
	
	if(nhdr.datatype == DT_BINARY || nhdr.datatype == DT_UNKNOWN  ) {
		NIFTI_ERROR("Bad datatype in header " << headerPath());
		return false;
	}
	if(nhdr.dim[1] <= 0) {
		NIFTI_ERROR("Bad first dimension in header " << headerPath());
		return false;
	}
	
	setDatatype(nhdr.datatype);
	// Set up dimensions. For now, trust dim[0] to hold the number of dimensions
	// and set trailing dims to sensible numbers
	for (int i = 0; i < nhdr.dim[0]; i++) {
		_dim[i] = nhdr.dim[i + 1];
		_voxdim[i] = nhdr.pixdim[i + 1];
	}
	for (int i = nhdr.dim[0]; i < 7; i++) {
		_dim[i] = 1;
		_voxdim[i] = 1.;
	}
	// Compute Q-Form
	Affine3f S; S = Scaling(_voxdim[0], _voxdim[1], _voxdim[2]);
	if( !is_nifti || nhdr.qform_code <= 0 ) {
		// If Q-Form not set or ANALYZE then just use voxel scaling
		_qform = S.matrix();
		qform_code = NIFTI_XFORM_UNKNOWN ;
	} else {
		float b = fixFloat(nhdr.quatern_b);
		float c = fixFloat(nhdr.quatern_c);
		float d = fixFloat(nhdr.quatern_d);
		float a = sqrt(1 - (b*b + c*c + d*d));
				
		float x = fixFloat(nhdr.qoffset_x);
		float y = fixFloat(nhdr.qoffset_y);
		float z = fixFloat(nhdr.qoffset_z);

		Quaternionf Q(a, b, c, d);
		Affine3f T; T = Translation3f(x, y, z);
		_qform = T*Q*S;
		
		// Fix left-handed co-ords in a very dumb way (see writeHeader())
		if (nhdr.pixdim[0] < 0.)
			_qform.matrix().block(0, 2, 3, 1) *= -1.;
		qform_code = nhdr.qform_code;
	}
	// Load S-Form
	if( !is_nifti || nhdr.sform_code <= 0 ) {
		sform_code = NIFTI_XFORM_UNKNOWN ;
	} else {
		_sform.setIdentity();
		for (int i = 0; i < 4; i++) {
			_sform(0, i) = nhdr.srow_x[i];
			_sform(1, i) = nhdr.srow_y[i];
			_sform(2, i) = nhdr.srow_z[i];
		}
		sform_code = nhdr.sform_code ;
	}
	
	if (is_nifti) {
		scaling_slope   = fixFloat(nhdr.scl_slope);
		if (scaling_slope == 0.)
			scaling_slope = 1.;
		scaling_inter   = fixFloat(nhdr.scl_inter);
		intent_code = nhdr.intent_code;
		intent_p1 = fixFloat(nhdr.intent_p1);
		intent_p2 = fixFloat(nhdr.intent_p2);
		intent_p3 = fixFloat(nhdr.intent_p3);
		toffset   = fixFloat(nhdr.toffset);
		
		intent_name = string(nhdr.intent_name);
		xyz_units = XYZT_TO_SPACE(nhdr.xyzt_units);
		time_units = XYZT_TO_TIME(nhdr.xyzt_units);
		freq_dim  = DIM_INFO_TO_FREQ_DIM (nhdr.dim_info);
		phase_dim = DIM_INFO_TO_PHASE_DIM(nhdr.dim_info);
		slice_dim = DIM_INFO_TO_SLICE_DIM(nhdr.dim_info);
		slice_code     = nhdr.slice_code;
		slice_start    = nhdr.slice_start;
		slice_end      = nhdr.slice_end;
		slice_duration = fixFloat(nhdr.slice_duration);
	}
	calibration_min = fixFloat(nhdr.cal_min);
	calibration_max = fixFloat(nhdr.cal_max);
	
	description = string(nhdr.descrip);
	aux_file    = string(nhdr.aux_file);
	
	if (_nii) {
		_voxoffset = (int)nhdr.vox_offset;
		if (_voxoffset < (int)sizeof(nhdr)) _voxoffset = (int)sizeof(nhdr);
	} else {
		_voxoffset = (int)nhdr.vox_offset ;
	}
		
	return true;
}

bool NiftiImage::writeHeader() {
	struct nifti_1_header nhdr;
	memset(&nhdr,0,sizeof(nhdr)) ;  /* zero out header, to be safe */
	/**- load the ANALYZE-7.5 generic parts of the header struct */
	nhdr.sizeof_hdr = sizeof(nhdr);
	nhdr.regular    = 'r';             /* for some stupid reason */
	
	nhdr.dim[0] = dimensions(); //pixdim[0] is set later with qform
	for (int i = 0; i < 7; i++) {	// Copy this way so types can be changed
		nhdr.dim[i + 1] = _dim[i];
		nhdr.pixdim[i + 1] = _voxdim[i];
	}
	
	nhdr.datatype = _datatype.code;
	nhdr.bitpix   = 8 * _datatype.size;
	
	if(calibration_max > calibration_min) {
		nhdr.cal_max = calibration_max;
		nhdr.cal_min = calibration_min;
	}
	
	if(scaling_slope != 0.0) {
		nhdr.scl_slope = scaling_slope;
		nhdr.scl_inter = scaling_inter;
	}
	
	strncpy(nhdr.descrip, description.c_str(), 80);
	strncpy(nhdr.aux_file, aux_file.c_str(), 24);
	

	if(_nii)
		strcpy(nhdr.magic,"n+1");
	else
		strcpy(nhdr.magic,"ni1");
	
	nhdr.intent_code = intent_code;
	nhdr.intent_p1   = intent_p1;
	nhdr.intent_p2   = intent_p2;
	nhdr.intent_p3   = intent_p3;
	strncpy(nhdr.intent_name, intent_name.c_str(), 16);
	
	// Check that _voxoffset is sensible
	if (_nii && (_voxoffset < nhdr.sizeof_hdr))
		_voxoffset = 352;
	nhdr.vox_offset = _voxoffset ;
	nhdr.xyzt_units = SPACE_TIME_TO_XYZT(xyz_units, time_units);
	nhdr.toffset    = toffset ;
	
	nhdr.qform_code = qform_code;
	Quaternionf Q(_qform.rotation());
	Translation3f T(_qform.translation());
	// Fix left-handed co-ord systems in an incredibly dumb manner.
	// First - NIFTI stores this information in pixdim[0], with both inconsistent
	// documentation and a reference implementation that hides pixdim[0] on reading
	// Second - Eigen .rotation() simultaneously calculates a scaling, and so may
	// hide axes flips. Hence we need to use .linear() to get the determinant
	if (_qform.linear().determinant() < 0)
		nhdr.pixdim[0] = -1.;
	else
		nhdr.pixdim[0] = 1.;
	// NIfTI REQUIRES a (or w) >= 0. Because Q and -Q represent the same
	// rotation, if w < 0 simply store -Q
	if (Q.w() < 0) {
		nhdr.quatern_b  = -Q.x();
		nhdr.quatern_c  = -Q.y();
		nhdr.quatern_d  = -Q.z();
	} else {
		nhdr.quatern_b = Q.x();
		nhdr.quatern_c = Q.y();
		nhdr.quatern_d = Q.z();
	}
	nhdr.qoffset_x  = T.x();
	nhdr.qoffset_y  = T.y();
	nhdr.qoffset_z  = T.z();
	
	nhdr.sform_code = sform_code;
	for (int i = 0; i < 4; i++) {
		nhdr.srow_x[i]  = _sform(0, i);
		nhdr.srow_y[i]  = _sform(1, i);
		nhdr.srow_z[i]  = _sform(2, i);
	}
	
	nhdr.dim_info = FPS_INTO_DIM_INFO(freq_dim, phase_dim, slice_dim);
	nhdr.slice_code     = slice_code;
	nhdr.slice_start    = slice_start;
	nhdr.slice_end      = slice_end;
	nhdr.slice_duration = slice_duration;
	
	if(_file.write(&nhdr, sizeof(nhdr)) < sizeof(nhdr)) {
		NIFTI_ERROR("Could not write header to file " + headerPath() + ".");
		return false;
	}
	return true;
}

/**
  *   Reads a sequence of bytes from the open NIfTI image.
  *
  *   Internal function to actually read bytes from an image file.
  *   @param start Location in file to start the read
  *   @param length Number of bytes to read
  *   @param buffer Location to read bytes to. If none, or NULL, specified,
  *          memory will be allocated internally and a pointer returned.
  *   @return On success a pointer to the read bytes (if buffer was specified
  *           then will be the same). NULL on fail.
  */
char *NiftiImage::readBytes(size_t start, size_t length, char *buffer)
{
	if (_mode == CLOSED) {
		NIFTI_ERROR("Cannot read from a closed file.");
		return NULL;
	}
	if (_mode == WRITE) {
		NIFTI_ERROR("Cannot read from a file opened for writing.");
		return NULL;
	}
	if (length == 0) {
		NIFTI_ERROR("Asked to read a buffer of 0 bytes.");
		return NULL;
	}
	bool didAllocate = false;
	if (!buffer) {
		buffer = new char[length];
		didAllocate = true;
	}
	_file.seek(_voxoffset + start, SEEK_SET);
	if (_file.read(buffer, static_cast<unsigned int>(length)) != length)
		NIFTI_ERROR("Read buffer returned wrong number of bytes.");
	
	if (_datatype.swapsize > 1 && _swap)
		SwapBytes(length / _datatype.swapsize, _datatype.swapsize, buffer);
	return buffer;
}

/**
  *   Writes a sequence of bytes to the open NIfTI image.
  *
  *   Internal function to actually write bytes to an image file.
  *   @param buffer Buffer to write bytes from.
  *   @param start Location in file to start the write
  *   @param length Number of bytes to write
  */
void NiftiImage::writeBytes(char *buffer, size_t start, size_t length)
{
	if (_mode == CLOSED) {
		NIFTI_ERROR("Cannot write to a closed file.");
		return;
	}
	if (_mode == READ) {
		NIFTI_ERROR("Cannot write to a file opened for writing.");
		return;
	}
	if (length < 1) {
		NIFTI_ERROR("Invalid write of length " << length << " bytes attempted. None written.");
		return;
	}
	_file.seek(_voxoffset + start, SEEK_SET);
	if (_file.write(buffer, static_cast<unsigned int>(length)) != length)
		NIFTI_ERROR("Write buffer failed.");
}

char *NiftiImage::readRawVolume(const int vol)
{
	size_t bytesPerVolume = voxelsPerVolume() * _datatype.size;
	char *raw = readBytes(vol * bytesPerVolume, bytesPerVolume);
	return raw;
}
char *NiftiImage::readRawAllVolumes()
{
	char *raw =	readBytes(0, voxelsTotal() * _datatype.size);
	return raw;
}
		
bool NiftiImage::open(const string &path, const char &mode) {
	_basepath = path.substr(0, path.find_first_of("."));
	if (path.substr(path.find_last_of(".") + 1) == "gz") {
		_gz = true;
	} else {
		_gz = false;
	}
	string ext = path.substr(path.find_first_of(".") + 1, 3);
	if (ext == "hdr" || ext == "img") {
		_nii = false;
	} else if (ext == "nii") {
		_nii = true;
	} else {
		NIFTI_ERROR("Invalid NIfTI extension: " + ext);
		return false;
	}
	
	if (_mode != CLOSED) {
		NIFTI_ERROR("Attempted to open file " << path <<
		           " using NiftiImage that is already open with file " + imagePath());
		return false;
	} else {
		_mode = mode;
		if ((mode == READ) || (mode == READ_HEADER)) {
			if(!_file.open(headerPath().c_str(), "rb", _gz)) {
				NIFTI_ERROR("Could not open header file " + headerPath() + " for reading.");
				return false;
			}
			if (!readHeader())
				return false;
			if (_mode == READ_HEADER) {
				close();
				return true;
			}
		} else if (mode == WRITE) {
			if(!_file.open(headerPath().c_str(), "wb", _gz)) {
				NIFTI_ERROR("Could not open header file " + headerPath() + " for writing.");
				return false;
			}
			if (!writeHeader())
				return false;
		} else {
			NIFTI_FAIL(string("Invalid NiftImage mode '") + mode + "'.");
		}
		_file.seek(_voxoffset, SEEK_SET);
		if (!_nii) {
			// Need to close the header and open the image
			_file.close();
			bool result;
			if (mode == READ)
				result = _file.open(imagePath().c_str(), "rb", _gz);
			else
				result = _file.open(imagePath().c_str(), "wb", _gz);
			if (!result) {
				NIFTI_ERROR("Could not open image file: " << imagePath());
				return false;
			}
		}
	}
	return true;
}

void NiftiImage::close()
{
	if (_mode == CLOSED) {
		NIFTI_ERROR("File is already closed: " << imagePath());
	} else if ((_mode == READ) || (_mode == READ_HEADER)) {
		_file.close();
		_mode = CLOSED;
	} else if (_mode == WRITE) {
		// If we've been writing subvolumes then we may not have written a complete file
		// Write a single zero-byte at the end to persuade the OS to write a file of the
		// correct size.
		_file.seek(0, SEEK_END);
		long correctEnd = (voxelsTotal() * _datatype.size + _voxoffset);
		char zero = 0;
		long pos = _file.tell();
		if (pos < correctEnd) {
			_file.seek(correctEnd - 1, SEEK_SET);
			_file.write(&zero, 1);
		}
		_file.flush();
		_file.close();
		_mode = CLOSED;
	}
}

int NiftiImage::dimensions() const {
	for (int d = 7; d > 0; d--) {
		if (_dim[d - 1] > 1) {
			return d;
		}
	}
	return 1;
}
void NiftiImage::setDimensions(const ArrayXi &dims, const ArrayXf &voxDims) {
	if (dims.rows() != voxDims.rows())
		NIFTI_FAIL("Length of dimension and voxel size arrays did not match.");
	if ((dims.rows() < 1) || (dims.rows() > 7))
		NIFTI_FAIL("Tried to set an invalid number of dimensions (" << dims.rows() << "). Valid numbers are 1 to 7.");
	_dim.head(dims.rows()) = dims;
	_voxdim.head(voxDims.rows()) = voxDims;
	for (int i = 0; i < _dim.rows(); i++) {
		if (_dim[i] < 1)
			_dim[i] = 1;
	}
}
	
int NiftiImage::dim(const int d) const {
	if ((d > 0) && (d <= _dim.rows())) {
		return _dim[d - 1];
	} else {
		NIFTI_ERROR("Tried to read invalid dimension: " << d << " from file: " << imagePath());
		return -1;
	}
}
void NiftiImage::setDim(const int d, const int n) {
	if (_mode == CLOSED) {
		if ((d > 0) && (d < 8)) {
			_dim[d - 1] = n;
		} else {
			NIFTI_ERROR("Tried to write invalid dimension: " << d << " for file: " << imagePath());
		}
	} else {
		NIFTI_FAIL("Cannot change image dimensions for open file: " << imagePath());
	}
}
const ArrayXi NiftiImage::dims() const { return _dim.head(dimensions()); }
void NiftiImage::setDims(const ArrayXi &n) {
	if (_mode == CLOSED) {
		if (n.rows() <= _dim.rows())
			_dim.head(n.rows()) = n;
		else
			NIFTI_ERROR("Tried to set an invalid number of dimensions for file: " << imagePath());
	} else {
		NIFTI_FAIL("Cannot change image dimensions for open file: " << imagePath());
	}
}

int NiftiImage::voxelsPerSlice() const  { return _dim[0]*_dim[1]; };
int NiftiImage::voxelsPerVolume() const { return _dim[0]*_dim[1]*_dim[2]; };
int NiftiImage::voxelsTotal() const     { return _dim.prod(); }

float NiftiImage::voxDim(const int d) const {
	if ((d > 0) && (d <= _voxdim.rows()))
		return _voxdim[d - 1];
	else
		NIFTI_FAIL("Tried to read voxel size for invalid dimension: " << d << " from file: " << imagePath());
}
void NiftiImage::setVoxDim(const int d, const float f) {
	if (_mode == CLOSED) {
		if ((d > 0) && (d <= _voxdim.rows()))
			_voxdim[d] = f;
		else
			NIFTI_FAIL("Tried to write voxel szie for invalid dimension: " << d << "from file: " << imagePath());
	} else
		NIFTI_FAIL("Cannot change voxel sizes for open file.");
}
const ArrayXf NiftiImage::voxDims() const { return _voxdim; }
void NiftiImage::setVoxDims(const ArrayXf &n) {
	if (_mode == CLOSED) {
		if (n.rows() <= _voxdim.rows())
			_voxdim.head(n.rows()) = n;
		else
			NIFTI_FAIL("Tried to set an invalid number of voxel sizes for file: " << imagePath());
	} else
		NIFTI_FAIL("Cannot change voxel sizes for open file: " << imagePath());
}

const int &NiftiImage::datatype() const { return _datatype.code; }
const string &NiftiImage::dtypeName() const { return _datatype.name; }
const int &NiftiImage::bytesPerVoxel() const { return _datatype.size; }
void NiftiImage::setDatatype(const int dt)
{
	/*  The map is declared here because making it a static member of NiftiImage was
	 *  causing problems with looking up the datatype in close() when called by 
	 *  ~NiftiImage. It's possible for C++ to destruct static members even when
	 *  objects still exist in another translation unit.
	 */
	static const DTMap DataTypes{
		{NIFTI_TYPE_UINT8,    {NIFTI_TYPE_UINT8, 1, 0, "NIFTI_TYPE_UINT8"} },
		{NIFTI_TYPE_INT16,    {NIFTI_TYPE_INT16, 2, 2, "NIFTI_TYPE_INT16"} },
		{NIFTI_TYPE_INT32,    {NIFTI_TYPE_INT32, 4, 4, "NIFTI_TYPE_INT32"} },
		{NIFTI_TYPE_FLOAT32,   {NIFTI_TYPE_FLOAT32, 4, 4, "NIFTI_TYPE_FLOAT32"} },
		{NIFTI_TYPE_COMPLEX64,   {NIFTI_TYPE_COMPLEX64, 8, 4, "NIFTI_TYPE_COMPLEX64"} },
		{NIFTI_TYPE_FLOAT64,   {NIFTI_TYPE_FLOAT64, 8, 8, "NIFTI_TYPE_FLOAT64"} },
		{NIFTI_TYPE_RGB24,  {NIFTI_TYPE_RGB24, 3, 0, "NIFTI_TYPE_RGB24"} },
		{NIFTI_TYPE_INT8,  {NIFTI_TYPE_INT8, 1, 0, "NIFTI_TYPE_INT8"} },
		{NIFTI_TYPE_UINT16,  {NIFTI_TYPE_UINT16, 2, 2, "NIFTI_TYPE_UINT16"} },
		{NIFTI_TYPE_UINT32,  {NIFTI_TYPE_UINT32, 4, 4, "NIFTI_TYPE_UINT32"} },
		{NIFTI_TYPE_INT64, {NIFTI_TYPE_INT64, 8, 8, "NIFTI_TYPE_INT64"} },
		{NIFTI_TYPE_UINT64, {NIFTI_TYPE_UINT64, 8, 8, "NIFTI_TYPE_UINT64"} },
		{NIFTI_TYPE_FLOAT128, {NIFTI_TYPE_FLOAT128, 16, 16, "NIFTI_TYPE_FLOAT128"} },
		{NIFTI_TYPE_COMPLEX128, {NIFTI_TYPE_COMPLEX128, 16,  8, "NIFTI_TYPE_COMPLEX128"} },
		{NIFTI_TYPE_COMPLEX256, {NIFTI_TYPE_COMPLEX256, 32, 16, "NIFTI_TYPE_COMPLEX256"} },
		{NIFTI_TYPE_RGBA32, {NIFTI_TYPE_RGBA32, 4,   0, "NIFTI_TYPE_RGBA32"} }
	};

	if (_mode == READ) {
		NIFTI_ERROR("Cannot set the datatype of a file opened for reading.");
		return;
	}
    DTMap::const_iterator it = DataTypes.find(dt);
	if (it == DataTypes.end())
		NIFTI_ERROR("Attempted to set invalid datatype.");
	else
		_datatype = it->second;
}

bool NiftiImage::matchesVoxels(const NiftiImage &other) const
{
	// Only check the first 3 dimensions
	if ((_dim.head(3) == other._dim.head(3)).all() && (_voxdim.head(3).isApprox(other._voxdim.head(3))))
		return true;
	else
		return false;
}

bool NiftiImage::matchesSpace(const NiftiImage &other) const
{
	if (matchesVoxels(other) && ijk_to_xyz().isApprox(other.ijk_to_xyz()))
		return true;
	else
		return false;	
}

const Matrix4f &NiftiImage::qform() const { return _qform.matrix(); }
const Matrix4f &NiftiImage::sform() const { return _sform.matrix(); }
const Matrix4f &NiftiImage::ijk_to_xyz() const
{
	if ((sform_code > 0) && (sform_code >= qform_code))
		return _sform.matrix();
	else // There is always a _qform matrix
		return _qform.matrix();
}
const Matrix4f &NiftiImage::xyz_to_ijk() const
{
	static Matrix4f inverse;
	if ((sform_code > 0) && (sform_code >= qform_code))
		inverse = _sform.matrix().inverse();
	else // There is always a _qform matrix
		inverse = _qform.matrix().inverse();
	return inverse;
}
