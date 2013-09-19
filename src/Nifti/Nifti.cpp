#include "Nifti/Nifti.h"
#include "Nifti/Internal.h"

const Nifti::DataType Nifti::DataTypeForCode(const int code) {
	static const map<int, DataType> c2dt{
		{NIFTI_TYPE_UINT8, DataType::UINT8},
		{NIFTI_TYPE_UINT16, DataType::UINT16},
		{NIFTI_TYPE_UINT32, DataType::UINT32},
		{NIFTI_TYPE_UINT64, DataType::UINT64},
		{NIFTI_TYPE_INT8, DataType::INT8},
		{NIFTI_TYPE_INT16, DataType::INT16},
		{NIFTI_TYPE_INT32, DataType::INT32},
		{NIFTI_TYPE_INT64, DataType::INT64},
		{NIFTI_TYPE_FLOAT32, DataType::FLOAT32},
		{NIFTI_TYPE_FLOAT64, DataType::FLOAT64},
		{NIFTI_TYPE_FLOAT128, DataType::FLOAT128},
		{NIFTI_TYPE_COMPLEX64, DataType::COMPLEX64},
		{NIFTI_TYPE_COMPLEX128, DataType::COMPLEX128},
		{NIFTI_TYPE_COMPLEX256, DataType::COMPLEX256},
		{NIFTI_TYPE_RGB24, DataType::RGB24},
		{NIFTI_TYPE_RGBA32, DataType::RGBA32},
	};
	auto dt = c2dt.find(code);
	if (dt != c2dt.end())
		return dt->second;
	else
		throw(std::invalid_argument("Unsupported data format code: " + to_string(code)));
}
	

/*	Internal map of datatype properties
 *
 *  The map is declared here because making it a static member of Nifti was
 *  causing problems with looking up the datatype in close() when called by 
 *  ~Nifti. It's possible for C++ to destruct static members even when
 *  objects still exist in another translation unit.
 */
const Nifti::DataTypeInfo &Nifti::TypeInfo(const DataType dt) {
	static const map<DataType, DataTypeInfo> DTInfo{
		{DataType::UINT8,      {DataType::UINT8, NIFTI_TYPE_UINT8, 1, 0, "NIFTI_TYPE_UINT8"} },
		{DataType::UINT16,     {DataType::UINT16, NIFTI_TYPE_UINT16, 2, 2, "NIFTI_TYPE_UINT16"} },
		{DataType::UINT32,     {DataType::UINT32, NIFTI_TYPE_UINT32, 4, 4, "NIFTI_TYPE_UINT32"} },
		{DataType::UINT64,     {DataType::UINT64, NIFTI_TYPE_UINT64, 8, 8, "NIFTI_TYPE_UINT64"} },
		{DataType::INT8,       {DataType::INT8, NIFTI_TYPE_INT8, 1, 0, "NIFTI_TYPE_INT8"} },
		{DataType::INT16,      {DataType::INT16, NIFTI_TYPE_INT16, 2, 2, "NIFTI_TYPE_INT16"} },
		{DataType::INT32,      {DataType::INT32, NIFTI_TYPE_INT32, 4, 4, "NIFTI_TYPE_INT32"} },
		{DataType::INT64,      {DataType::INT64, NIFTI_TYPE_INT64, 8, 8, "NIFTI_TYPE_INT64"} },
		{DataType::FLOAT32,    {DataType::FLOAT32, NIFTI_TYPE_FLOAT32, 4, 4, "NIFTI_TYPE_FLOAT32"} },
		{DataType::FLOAT64,    {DataType::FLOAT64, NIFTI_TYPE_FLOAT64, 8, 8, "NIFTI_TYPE_FLOAT64"} },
		{DataType::FLOAT128,   {DataType::FLOAT128, NIFTI_TYPE_FLOAT128, 16, 16, "NIFTI_TYPE_FLOAT128"} },
		{DataType::COMPLEX64,  {DataType::COMPLEX64, NIFTI_TYPE_COMPLEX64, 8, 4, "NIFTI_TYPE_COMPLEX64"} },
		{DataType::COMPLEX128, {DataType::COMPLEX128, NIFTI_TYPE_COMPLEX128, 16,  8, "NIFTI_TYPE_COMPLEX128"} },
		{DataType::COMPLEX256, {DataType::COMPLEX256, NIFTI_TYPE_COMPLEX256, 32, 16, "NIFTI_TYPE_COMPLEX256"} },
		{DataType::RGB24,      {DataType::RGB24, NIFTI_TYPE_RGB24, 3, 0, "NIFTI_TYPE_RGB24"} },
		{DataType::RGBA32,     {DataType::RGBA32, NIFTI_TYPE_RGBA32, 4,   0, "NIFTI_TYPE_RGBA32"} }
	};
	auto info = DTInfo.find(dt);
	if (info != DTInfo.end())
		return info->second;
	else
		throw(std::invalid_argument("Missing type information, contact libNifti author."));
}

/*
 * Returns the string representation of a NIfTI space unit code.
 *
 *\sa NIFTI1_UNITS group in nifti1.h
 */
const string &Nifti::spaceUnits() const {
	static const map<int, string> Units {
		{ NIFTI_UNITS_METER,  "m" },
		{ NIFTI_UNITS_MM,     "mm" },
		{ NIFTI_UNITS_MICRON, "um" }
	};
	static const string unknown("Unknown space units code");
	auto it = Units.find(xyz_units);
	if (it == Units.end())
		return unknown;
	else
		return it->second;
}
/*
 * Returns the string description of a NIfTI time unit code.
 *
 *\sa NIFTI1_UNITS group in nifti1.h
 */
const string &Nifti::timeUnits() const {
	static const map<int, string> Units {
		{ NIFTI_UNITS_SEC,    "s" },
		{ NIFTI_UNITS_MSEC,   "ms" },
		{ NIFTI_UNITS_USEC,   "us" },
		{ NIFTI_UNITS_HZ,     "Hz" },
		{ NIFTI_UNITS_PPM,    "ppm" },
		{ NIFTI_UNITS_RADS,   "rad/s" }
	};
	static const string unknown("Unknown time units code");
	auto it = Units.find(time_units);
	if (it == Units.end())
		return unknown;
	else
		return it->second;
}

/*
 * Returns the string representation of NIfTI intent code.
 *
 *\sa NIFTI1_INTENT_CODES group in nifti1.h
 */
const string &Nifti::intentName() const {
	static const map<int, string> Intents {
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
	static const string unknown("Unknown intent code");
	auto it = Intents.find(intent_code);
	if (it == Intents.end())
		return unknown;
	else
		return it->second;
}

/*
 * Returns the string representation of a NIfTI transform code.
 *
 *\sa NIFTI1_XFORM_CODES group in nifti1.h
 */
const string Nifti::XFormName(const Nifti::XForm c) {
	switch (c) {
		case XForm::Unknown: return "Unknown"; break;
		case XForm::ScannerAnatomy: return "Scanner Anatomy"; break;
		case XForm::AlignedAnatomy: return "Aligned Anatomy"; break;
		case XForm::Talairach: return "Talairach"; break;
		case XForm::MNI_152: return "MNI 152"; break;
	}
}
const int Nifti::XFormCode(const Nifti::XForm c) {
	switch (c) {
		case XForm::Unknown: return NIFTI_XFORM_UNKNOWN; break;
		case XForm::ScannerAnatomy: return NIFTI_XFORM_SCANNER_ANAT; break;
		case XForm::AlignedAnatomy: return NIFTI_XFORM_ALIGNED_ANAT; break;
		case XForm::Talairach: return NIFTI_XFORM_TALAIRACH; break;
		case XForm::MNI_152: return NIFTI_XFORM_MNI_152; break;
	}
}
const Nifti::XForm Nifti::XFormForCode(const int c) {
	switch (c) {
		case NIFTI_XFORM_UNKNOWN: return XForm::Unknown; break;
		case NIFTI_XFORM_SCANNER_ANAT: return XForm::ScannerAnatomy; break;
		case NIFTI_XFORM_ALIGNED_ANAT: return XForm::AlignedAnatomy; break;
		case NIFTI_XFORM_TALAIRACH: return XForm::Talairach; break;
		case NIFTI_XFORM_MNI_152: return XForm::MNI_152; break;
		default:
			throw(std::invalid_argument("Invalid transform code: " + to_string(c)));
	}
}

/*
 * Returns the string representation of a NIfTI slice_code
 *
 *\sa NIFTI1_SLICE_ORDER group in nifti1.h
 */
const string &Nifti::sliceName() const {
	static const map<int, string> SliceOrders {
		{ NIFTI_SLICE_SEQ_INC,  "sequential_increasing"    },
		{ NIFTI_SLICE_SEQ_DEC,  "sequential_decreasing"    },
		{ NIFTI_SLICE_ALT_INC,  "alternating_increasing"   },
		{ NIFTI_SLICE_ALT_DEC,  "alternating_decreasing"   },
		{ NIFTI_SLICE_ALT_INC2, "alternating_increasing_2" },
		{ NIFTI_SLICE_ALT_DEC2, "alternating_decreasing_2" }
	};
	static const string unknown("Unknown slice order code");
	auto it = SliceOrders.find(slice_code);
	if (it == SliceOrders.end())
		return unknown;
	else
		return it->second;
}

Nifti::~Nifti()
{
	if (m_mode != Mode::Closed)
		close();
}

Nifti::Nifti() :
	m_mode(Mode::Closed), m_gz(false), m_nii(false), m_swap(false), m_voxoffset(0),
	m_dim(Array<size_t, 7, 1>::Ones()), m_voxdim(Array<float, 7, 1>::Ones()),
	m_basepath(""), m_typeinfo(TypeInfo(DataType::FLOAT32)),
	m_qcode(XForm::Unknown), m_scode(XForm::Unknown),
	scaling_slope(1.), scaling_inter(0.), calibration_min(0.), calibration_max(0.),
	freq_dim(0), phase_dim(0), slice_dim(0),
	slice_code(0), slice_start(0), slice_end(0), slice_duration(0),
	toffset(0), xyz_units(NIFTI_UNITS_MM), time_units(NIFTI_UNITS_SEC),
	intent_code(NIFTI_INTENT_NONE), intent_p1(0), intent_p2(0), intent_p3(0),
	intent_name(""), description(""), aux_file("")
{
	m_qform.setIdentity(); m_sform.setIdentity();
}

Nifti::Nifti(const Nifti &other) :
	m_mode(other.m_mode), m_gz(other.m_gz), m_nii(other.m_nii),
	m_swap(other.m_swap), m_voxoffset(other.m_voxoffset),
	m_dim(other.m_dim), m_voxdim(other.m_voxdim),
	m_qform(other.m_qform), m_sform(other.m_sform),
	m_qcode(other.m_qcode), m_scode(other.m_scode),
	m_typeinfo(other.m_typeinfo), m_basepath(other.m_basepath), m_file(),
	m_extensions(other.m_extensions),
	scaling_slope(other.scaling_slope), scaling_inter(other.scaling_inter),
	calibration_min(other.calibration_min), calibration_max(other.calibration_max),
	freq_dim(other.freq_dim), phase_dim(other.phase_dim), slice_dim(other.slice_dim),
	slice_code(other.slice_code), slice_start(other.slice_start),
	slice_end(other.slice_end), slice_duration(other.slice_duration),
	toffset(other.toffset), xyz_units(other.xyz_units), time_units(other.time_units),
	intent_code(other.intent_code), intent_p1(other.intent_p1), intent_p2(other.intent_p2), intent_p3(other.intent_p3),
	intent_name(other.intent_name), description(other.description), aux_file(other.aux_file)
{
	if ((m_mode == Mode::Read) || (m_mode == Mode::ReadSkipExt)) {
		m_file.open(imagePath(), "rb", m_gz);
		m_file.seek(other.m_file.tell(), SEEK_SET);
	} else if ((m_mode == Mode::Write) || (m_mode == Mode::WriteSkipExt)) {
		m_file.open(imagePath(), "wb", m_gz);
		m_file.seek(other.m_file.tell(), SEEK_SET);
	}
}

Nifti::Nifti(const Nifti &other, const size_t nt, const DataType dtype) : Nifti() {
	m_dim.head(3) = other.m_dim.head(3);
	setDim(4, nt);
	m_voxdim.head(3) = other.m_voxdim.head(3);
	m_qform = other.m_qform; m_sform = other.m_sform;
	m_qcode = other.m_qcode; m_scode = other.m_scode;
	xyz_units = other.xyz_units;
	m_typeinfo = TypeInfo(dtype);
}

Nifti::Nifti(Nifti &&other) noexcept :
	m_mode(other.m_mode), m_gz(other.m_gz), m_nii(other.m_nii),
	m_swap(other.m_swap), m_voxoffset(other.m_voxoffset),
	m_dim(other.m_dim), m_voxdim(other.m_voxdim),
	m_qform(other.m_qform), m_sform(other.m_sform),
	m_qcode(other.m_qcode), m_scode(other.m_scode),
	m_typeinfo(other.m_typeinfo), m_basepath(other.m_basepath), m_file(other.m_file),
	m_extensions(other.m_extensions),
	scaling_slope(other.scaling_slope), scaling_inter(other.scaling_inter),
	calibration_min(other.calibration_min), calibration_max(other.calibration_max),
	freq_dim(other.freq_dim), phase_dim(other.phase_dim), slice_dim(other.slice_dim),
	slice_code(other.slice_code), slice_start(other.slice_start),
	slice_end(other.slice_end), slice_duration(other.slice_duration),
	toffset(other.toffset), xyz_units(other.xyz_units), time_units(other.time_units),
	intent_code(other.intent_code), intent_p1(other.intent_p1), intent_p2(other.intent_p2), intent_p3(other.intent_p3),
	intent_name(other.intent_name), description(other.description), aux_file(other.aux_file)
{
	other.m_mode = Mode::Closed;
}

Nifti::Nifti(const string &filename, const Mode &mode) :
	Nifti()
{
	open(filename, mode);
}

Nifti::Nifti(const int nx, const int ny, const int nz, const int nt,
		     const float dx, const float dy, const float dz, const float dt,
			 const DataType dtype, const Affine3f &xform) :
	Nifti()
{
	m_typeinfo = TypeInfo(dtype);
	m_dim[0] = nx < 1 ? 1 : nx;
	m_dim[1] = ny < 1 ? 1 : ny;
	m_dim[2] = nz < 1 ? 1 : nz;
	m_dim[3] = nt < 1 ? 1 : nt;
	m_voxdim[0] = dx; m_voxdim[1] = dy; m_voxdim[2] = dz; m_voxdim[3] = dt;
	setTransform(xform);
}

Nifti::Nifti(const ArrayXs &dim, const ArrayXf &voxdim,
             const DataType dtype, const Affine3f &xform) :
	Nifti()
{
	assert(dim.rows() < 8);
	assert(dim.rows() == voxdim.rows());
	
	m_dim.head(dim.rows()) = dim;
	m_voxdim.head(voxdim.rows()) = voxdim;
	m_typeinfo = TypeInfo(dtype);
	setTransform(xform);
}

Nifti &Nifti::operator=(const Nifti &other)
{
	if (this == &other)
		return *this;
	else if (m_mode != Mode::Closed)
		close();
	
	m_dim = other.m_dim;
	m_voxdim = other.m_voxdim;
	m_qform = other.m_qform;
	m_sform = other.m_sform;
	m_qcode = other.m_qcode;
	m_scode = other.m_scode;
	m_basepath = other.m_basepath;
	m_gz = other.m_gz;
	m_nii = other.m_nii;
	m_mode = Mode::Closed;
	m_voxoffset = 0;
	m_typeinfo = other.m_typeinfo;
	scaling_slope = other.scaling_slope;
	scaling_inter = other.scaling_inter;
	calibration_min = other.calibration_min;
	calibration_max = other.calibration_max;
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

const string Nifti::basePath() const { return m_basepath; }
const string Nifti::imagePath() const {
	string path(m_basepath);
	if (m_nii) {
		path += ".nii";
	} else {
		path += ".img";
	}
	if (m_gz)
		path += ".gz";
	
	return path;
}
const string Nifti::headerPath() const {
	string path(m_basepath);
	if (m_nii) {
		path += ".nii";
	} else {
		path += ".hdr";
	}
	if (m_gz)
		path += ".gz";
	
	return path;
}

void Nifti::readHeader() {
	struct nifti_1_header nhdr;
	
	if (m_file.read(&nhdr, sizeof(nhdr)) < sizeof(nhdr)) {
		throw(std::runtime_error("Could not read header structure from " + headerPath()));
	}
	
	// Check if disk and CPU byte order match.
	// The sizeof_hdr field should always be 352, as this is the size of a
	// NIfTI-1 header
	if (nhdr.sizeof_hdr != sizeof(nhdr)) {
		swapBytes(1, 4, &nhdr.sizeof_hdr);
		if (nhdr.sizeof_hdr != sizeof(nhdr)) {
			throw(std::runtime_error("Could not determine byte order of header " + headerPath()));
		}
		// If we didn't fail, then we need to swap the header (so swap sizeof back)
		m_swap = true;
		swapBytes(1, 4, &nhdr.sizeof_hdr);
	}
	
	// Check the magic string is set to one of the possible NIfTI values,
	// otherwise process as an ANALYZE file
	int is_nifti = ((nhdr.magic[0]=='n' && nhdr.magic[3]=='\0') &&
                    (nhdr.magic[1]=='i' || nhdr.magic[1]=='+') &&
                    (nhdr.magic[2]>='1' && nhdr.magic[2]<='9')) ? true : false;
	
	if (m_swap && is_nifti)
		swapNiftiHeader(&nhdr);
	else if (m_swap)
		swapAnalyzeHeader((nifti_analyze75 *)&nhdr);
	
	if(nhdr.datatype == DT_BINARY || nhdr.datatype == DT_UNKNOWN  ) {
		throw(std::runtime_error("Bad datatype in header: " + headerPath()));
	}
	m_typeinfo = TypeInfo(DataTypeForCode(nhdr.datatype));
	
	if(nhdr.dim[1] <= 0) {
		throw(std::runtime_error("Bad first dimension in header: " + headerPath()));
	}
	for (int i = 0; i < nhdr.dim[0]; i++) {
		m_dim[i] = nhdr.dim[i + 1];
		m_voxdim[i] = nhdr.pixdim[i + 1];
	}
	for (int i = nhdr.dim[0]; i < 7; i++) {
		m_dim[i] = 1;
		m_voxdim[i] = 1.;
	}
	calcStrides();
	// Compute Q-Form
	Affine3f S; S = Scaling(m_voxdim[0], m_voxdim[1], m_voxdim[2]);
	if( !is_nifti || nhdr.qform_code <= 0 ) {
		// If Q-Form not set or ANALYZE then just use voxel scaling
		m_qform = S.matrix();
		m_qcode = XForm::Unknown;
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
		m_qform = T*Q*S;
		
		// Fix left-handed co-ords in a very dumb way (see writeHeader())
		if (nhdr.pixdim[0] < 0.)
			m_qform.matrix().block(0, 2, 3, 1) *= -1.;
		m_qcode = XFormForCode(nhdr.qform_code);
	}
	// Load S-Form
	if( !is_nifti || nhdr.sform_code <= 0 ) {
		m_sform = S.matrix();
		m_scode = XForm::Unknown;
	} else {
		m_sform.setIdentity();
		for (int i = 0; i < 4; i++) {
			m_sform(0, i) = nhdr.srow_x[i];
			m_sform(1, i) = nhdr.srow_y[i];
			m_sform(2, i) = nhdr.srow_z[i];
		}
		m_scode = XFormForCode(nhdr.sform_code);
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
	
	if (m_nii) {
		m_voxoffset = (int)nhdr.vox_offset;
		if (m_voxoffset < (int)sizeof(nhdr)) m_voxoffset = (int)sizeof(nhdr);
	} else {
		m_voxoffset = (int)nhdr.vox_offset ;
	}
}

void Nifti::readExtensions()
{
	long target = m_voxoffset;
	if (!m_nii) {
		m_file.seek(0, SEEK_END);
		target = m_file.tell();
	}
	m_file.seek(sizeof(nifti_1_header), SEEK_SET);
	char extender[4];
	if (m_file.read(extender, 4) != 4) {
		throw(std::runtime_error("While checking for extensions hit end of file: " + headerPath()));
	}
	if (extender[0] != 1) // There are no extensions
		return;
	
	while (m_file.tell() < target) {		
		if(m_file.tell() > target - 16 ){
			throw(std::runtime_error("Insufficient space for remaining extensions in file: " + headerPath()));
		}
		
		int size, code;
		long bytesRead = m_file.read(&size, 4);
		bytesRead += m_file.read(&code, 4);
		if (bytesRead != 8) {
			throw(std::runtime_error("Error while reading extension size and code in file: " + headerPath()));
		}
		
		if (m_swap) {
			swapBytes(1, 4, &size);
			swapBytes(1, 4, &code);
		}
		
		vector<char> dataBytes(size - 8);
		if (m_file.read(dataBytes.data(), size - 8) < (size - 8)) {
			throw(std::runtime_error("Could not read extension in file: " + headerPath()));
		}
		m_extensions.emplace_back(code, dataBytes);

		if (m_nii && (m_file.tell() > m_voxoffset)) {
			throw(std::runtime_error("Went past start of voxel data while reading extensions in file: " + headerPath()));
		}
	}
}

void Nifti::addExtension(const int code, const vector<char> &data) {
	m_extensions.emplace_back(code, data);
}

void Nifti::addExtension(const Nifti::Extension &e) {
	m_extensions.push_back(e);
}

const list<Nifti::Extension> &Nifti::extensions() const {
	return m_extensions;
}

void Nifti::writeHeader() {
	struct nifti_1_header nhdr;
	memset(&nhdr,0,sizeof(nhdr)) ;  /* zero out header, to be safe */
	/**- load the ANALYZE-7.5 generic parts of the header struct */
	nhdr.sizeof_hdr = sizeof(nhdr);
	nhdr.regular    = 'r';             /* for some stupid reason */
	
	nhdr.dim[0] = dimensions(); //pixdim[0] is set later with qform
	for (size_t i = 0; i < 7; i++) { // Convert back to int/float
		if (m_dim[i] > numeric_limits<short>::max()) {
			throw(std::runtime_error("Nifti does not support dimensions greater than " + to_string(numeric_limits<short>::max())));
		}
		nhdr.dim[i + 1] = m_dim[i];
		nhdr.pixdim[i + 1] = m_voxdim[i];
	}
	
	nhdr.datatype = m_typeinfo.code;
	nhdr.bitpix   = 8 * m_typeinfo.size;
	
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
	

	if(m_nii)
		strcpy(nhdr.magic,"n+1");
	else
		strcpy(nhdr.magic,"ni1");
	
	nhdr.intent_code = intent_code;
	nhdr.intent_p1   = intent_p1;
	nhdr.intent_p2   = intent_p2;
	nhdr.intent_p3   = intent_p3;
	strncpy(nhdr.intent_name, intent_name.c_str(), 16);
	
	// Check that m_voxoffset is sensible
	m_voxoffset = 352;
	if (m_nii && (m_mode != Mode::WriteSkipExt))
		m_voxoffset += totalExtensionSize();
	nhdr.vox_offset = m_voxoffset ;
	nhdr.xyzt_units = SPACE_TIME_TO_XYZT(xyz_units, time_units);
	nhdr.toffset    = toffset ;
	
	nhdr.qform_code = XFormCode(m_qcode);
	Quaternionf Q(m_qform.rotation());
	Translation3f T(m_qform.translation());
	// Fix left-handed co-ord systems in an incredibly dumb manner.
	// First - NIFTI stores this information in pixdim[0], with both inconsistent
	// documentation and a reference implementation that hides pixdim[0] on reading
	// Second - Eigen .rotation() simultaneously calculates a scaling, and so may
	// hide axes flips. Hence we need to use .linear() to get the determinant
	if (m_qform.linear().determinant() < 0)
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
	
	nhdr.sform_code = XFormCode(m_scode);
	for (int i = 0; i < 4; i++) {
		nhdr.srow_x[i]  = m_sform(0, i);
		nhdr.srow_y[i]  = m_sform(1, i);
		nhdr.srow_z[i]  = m_sform(2, i);
	}
	
	nhdr.dim_info = FPS_INTO_DIM_INFO(freq_dim, phase_dim, slice_dim);
	nhdr.slice_code     = slice_code;
	nhdr.slice_start    = slice_start;
	nhdr.slice_end      = slice_end;
	nhdr.slice_duration = slice_duration;
	
	if(m_file.write(&nhdr, sizeof(nhdr)) < sizeof(nhdr)) {
		throw(std::runtime_error("Could not write header to file: " + headerPath()));
	}
}

int Nifti::totalExtensionSize() {
	int total = 0;
	for (auto ext: m_extensions) {
		total += ext.size();
	}
	return total;
}

void Nifti::writeExtensions() {
	m_file.seek(sizeof(nifti_1_header), SEEK_SET);
	char extender[4] = {0, 0, 0, 0};
	if (m_extensions.size() > 0)
		extender[0] = 1;
	if (m_file.write(extender, 4) < 4) {
		throw(std::runtime_error("Could not write extender block to file: " + headerPath()));
	}
	
	for (auto ext : m_extensions) {
		int size = ext.size();
		int padding = ext.padding();
		long bytesWritten = m_file.write(&size, sizeof(int));
		int code = ext.code();
		bytesWritten += m_file.write(&code, sizeof(int));
		if (bytesWritten != (2*sizeof(int))) {
			throw(std::runtime_error("Could not write extension size and code to file: " + headerPath()));
		}
		if (m_file.write(ext.data().data(), ext.rawSize()) != (ext.rawSize())) {
			throw(std::runtime_error("Could not write extension data to file: " + headerPath()));
		}
		if (padding) {
			vector<char> pad(ext.padding(), 0);
			if (m_file.write(pad.data(), ext.padding()) != ext.padding()) {
				throw(std::runtime_error("Could not write extension padding to file: " + headerPath()));
			}
		}
	}
	if ((m_file.tell() - totalExtensionSize() - 4) != sizeof(nifti_1_header)) {
		throw(std::runtime_error("Wrote wrong number of bytes for extensions to file: " + headerPath()));
	}
}

/**
  *   Simple function to calculate the strides into the data on disk. Used for
  *   subvolume/voxel-wise reads.
  */
void Nifti::calcStrides() {
	m_strides = Array<size_t, 7, 1>::Ones();
	for (size_t i = 1; i < dimensions(); i++) {
		m_strides(i) = m_strides(i - 1) * m_dim(i - 1);
	}
}

/**
  * Seeks to a particular voxel on the disk.
  *
  * @param target Desired voxel to seek to on disk.
  *
  * @throws std::out_of_range if the target is outside the image dimensions.
  * @throws std::runtime_error if the seek fails.
  */
void Nifti::seekToVoxel(const ArrayXs &target) {
	if (target.rows() > dimensions()) {
		throw(std::out_of_range("Too many dimensions for seeking."));
	}
	if ((target > m_dim.head(target.rows())).any()) {
		throw(std::out_of_range("Target voxel is outside image dimensions."));
	}
	size_t index = (target * m_strides.head(target.rows())).sum() * m_typeinfo.size + m_voxoffset;
	size_t current = m_file.tell();
	if (!m_file.seek(current - index, SEEK_CUR)) {
		throw(std::runtime_error("Failed to seek to index: " + to_string(index) + " in file: " + imagePath()));
	}
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
char *Nifti::readBytes(size_t start, size_t length, char *buffer) {
	if (m_mode == Mode::Closed) {
		throw(std::logic_error("Cannot read from closed file: " + imagePath()));
	}
	if (m_mode == Mode::Write) {
		throw(std::logic_error("Cannot read from a file opened for reading: " + imagePath()));
	}
	if (length == 0) {
		throw(std::invalid_argument("Asked to read 0 bytes from file: " + imagePath()));
	}
	if (!buffer) {
		buffer = new char[length];
	}
	if (!m_file.seek(m_voxoffset + start, SEEK_SET)) {
		throw(std::runtime_error("Failed seek in file: " + imagePath()));
	}
	if (m_file.read(buffer, static_cast<unsigned int>(length)) != length) {
		throw(std::runtime_error("Read wrong number of bytes from file: " + imagePath()));
	}
	if (m_typeinfo.swapsize > 1 && m_swap)
		swapBytes(length / m_typeinfo.swapsize, m_typeinfo.swapsize, buffer);
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
void Nifti::writeBytes(size_t start, size_t length, char *buffer) {
	if (m_mode == Mode::Closed) {
		throw(std::logic_error("Cannot write to closed file: " + imagePath()));
	}
	if (m_mode == Mode::Read) {
		throw(std::logic_error("Cannot write to file opened for writing: " + imagePath()));
	}
	if (length == 0) {
		throw(std::invalid_argument("Asked to write 0 bytes to file: " + imagePath()));
	}
	if (!m_file.seek(m_voxoffset + start, SEEK_SET)) {
		throw(std::runtime_error("Failed seek in file: " + imagePath()));
	}
	if (m_file.write(buffer, static_cast<unsigned int>(length)) != length) {
		throw(std::runtime_error("Wrote wrong number of bytes to file: " + imagePath()));
	}
}

void Nifti::open(const string &path, const Mode &mode) {
	size_t lastDot = path.find_last_of(".");
	string ext;
	if (path.substr(lastDot + 1) == "gz") {
		m_gz = true;
		size_t extDot = path.find_last_of(".", lastDot - 1);
		ext = path.substr(extDot + 1, lastDot - extDot - 1);
		m_basepath = path.substr(0, extDot);
	} else {
		m_gz = false;
		ext = path.substr(lastDot + 1);
		m_basepath = path.substr(0, lastDot);
	}
	if (ext == "hdr" || ext == "img") {
		m_nii = false;
	} else if (ext == "nii") {
		m_nii = true;
	} else {
		throw(std::invalid_argument("Invalid NIfTI extension for file: " + path));
	}
	
	if (m_mode != Mode::Closed) {
		throw(std::logic_error("Attempted to open file: " + path +
		           " when file: " + imagePath() + " is already open."));
	} else {
		if ((mode == Mode::Read) || (mode == Mode::ReadHeader) || (mode == Mode::ReadSkipExt)) {
			if(!m_file.open(headerPath(), "rb", m_gz)) {
				throw(std::runtime_error("Failed to open file: " + headerPath()));
			}
			readHeader();
			if (mode != Mode::ReadSkipExt ) {
				readExtensions();
			}
		} else if (mode == Mode::Write) {
			if(!m_file.open(headerPath(), "wb", m_gz)) {
				throw(std::runtime_error("Failed to open file: " + headerPath()));
			}
			writeHeader();
			if (mode == Mode::Write) {
				writeExtensions();
			}
		} else {
			throw(std::invalid_argument("Invalid opening mode for file: " + path));
		}
		
		if (mode == Mode::ReadHeader) {
			// Don't do anything in this case
		} else {
			if (!m_nii) {
				// Need to close the header and open the image
				m_file.close();
				bool result;
				if (mode == Mode::Read)
					result = m_file.open(imagePath(), "rb", m_gz);
				else
					result = m_file.open(imagePath(), "wb", m_gz);
				if (!result) {
					throw(std::runtime_error("Could not open image file: " + imagePath()));
				}
			}
			if (!m_file.seek(m_voxoffset, SEEK_SET)) {
				throw(std::runtime_error("Could not seek to voxel offset in file: " + imagePath()));
			}
			// Only set the mode here when we have successfully opened the file and
			// not thrown any errors. Throwing an error triggers the destructor and
			// we don't want to be in the wrong state there.
			m_mode = mode;
		}
	}
}

bool Nifti::isOpen() {
	if (m_mode == Mode::Closed)
		return false;
	else
		return true;
}
void Nifti::close()
{
	if (m_mode == Mode::Closed) {
		throw(std::logic_error("Cannot close already closed file: " + imagePath()));
	} else if ((m_mode == Mode::Read) || (m_mode == Mode::ReadHeader)) {
		m_file.close();
		m_mode = Mode::Closed;
	} else if (m_mode == Mode::Write) {
		// If we've been writing subvolumes then we may not have written a complete file
		// Write a single zero-byte at the end to persuade the OS to write a file of the
		// correct size.
		m_file.seek(0, SEEK_END);
		long correctEnd = (voxelsTotal() * m_typeinfo.size + m_voxoffset);
		char zero{0};
		long pos = m_file.tell();
		if (pos < correctEnd) {
			m_file.seek(correctEnd - 1, SEEK_SET);
			m_file.write(&zero, 1);
		}
		m_file.flush();
		m_file.close();
		m_mode = Mode::Closed;
	}
}

size_t Nifti::dimensions() const {
	for (size_t d = m_voxdim.rows(); d > 0; d--) {
		if (m_dim[d - 1] > 1) {
			return d;
		}
	}
	return 1;
}
	
size_t Nifti::dim(const size_t d) const {
	assert((d > 0) && (d <= m_voxdim.rows()));
	return m_dim[d - 1];
}
void Nifti::setDim(const size_t d, const size_t n) {
	if (m_mode == Mode::Closed) {
		assert((d > 0) && (d < 8));
		m_dim[d - 1] = n;
		calcStrides();
	} else {
		throw(std::logic_error("Cannot change image dimensions for open file: " + imagePath()));
	}
}
const ArrayXs Nifti::dims() const { return m_dim.head(dimensions()); }
void Nifti::setDims(const ArrayXs &n) {
	if (m_mode == Mode::Closed) {
		assert(n.rows() <= m_voxdim.rows());
		m_dim.head(n.rows()) = n;
		calcStrides();
	} else {
		throw(std::logic_error("Cannot change image dimensions for open file: " + imagePath()));
	}
}

size_t Nifti::voxelsPerSlice() const  { return m_dim[0]*m_dim[1]; };
size_t Nifti::voxelsPerVolume() const { return m_dim[0]*m_dim[1]*m_dim[2]; };
size_t Nifti::voxelsTotal() const     { return m_dim.prod(); }

float Nifti::voxDim(const size_t d) const {
	assert((d > 0) && (d <= m_voxdim.rows()));
	return m_voxdim[d - 1];
}
void Nifti::setVoxDim(const size_t d, const float f) {
	if (m_mode == Mode::Closed) {
		assert((d > 0) && (d <= m_voxdim.rows()));
		m_voxdim[d] = f;
	} else
		throw(std::logic_error("Cannot change voxel sizes for open file: " + imagePath()));
}
const ArrayXf Nifti::voxDims() const { return m_voxdim; }
void Nifti::setVoxDims(const ArrayXf &n) {
	if (m_mode == Mode::Closed) {
		assert(n.rows() <= m_voxdim.rows());
		m_voxdim.head(n.rows()) = n;
	} else
		throw(std::logic_error("Cannot change voxel sizes for open file: " + imagePath()));
}

const Nifti::DataType &Nifti::datatype() const { return m_typeinfo.type; }
void Nifti::setDatatype(const Nifti::DataType dt) {
	if (m_mode != Mode::Closed) {
		throw(std::logic_error("Cannot set the datatype of open file: " + imagePath()));
		return;
	}
    m_typeinfo = TypeInfo(dt);
}

bool Nifti::matchesVoxels(const Nifti &other) const {
	// Only check the first 3 dimensions
	if ((m_dim.head(3) == other.m_dim.head(3)).all() && (m_voxdim.head(3).isApprox(other.m_voxdim.head(3))))
		return true;
	else
		return false;
}

bool Nifti::matchesSpace(const Nifti &other) const {
	if (matchesVoxels(other) && transform().isApprox(other.transform()))
		return true;
	else
		return false;	
}

const Affine3f &Nifti::qform() const { return m_qform; }
const Affine3f &Nifti::sform() const { return m_sform; }
const Nifti::XForm Nifti::qcode() const { return m_qcode; }
const Nifti::XForm Nifti::scode() const { return m_scode; }
const Affine3f &Nifti::transform() const {
	if ((m_scode > XForm::Unknown) && (m_scode >= m_qcode))
		return m_sform;
	else // There is always a m_qform matrix
		return m_qform;
}
void Nifti::setTransform(const Affine3f &t, const XForm tc) {
	m_qform = t.linear();
	m_sform = t;
	m_qcode = tc;
	m_scode = tc;
}
