/** \file Extension.h
 \brief Declaration for Nifti Extension class
 - Written by Tobias Wood, IoP KCL
 - Based on nifti1_io.h (Thanks to Robert Cox et al)
 - This code is released to the public domain. Do with it what you will.
 */

#ifndef NIFTI_EXTENSION
#define NIFTI_EXTENSION

#include <vector>
#include <string>
#include <map>

using namespace std;

namespace Nifti {

#pragma mark Extension Codes
/* NIfTI-1.1 extension codes: see http://nifti.nimh.nih.gov/nifti-1/documentation/faq#Q21 */

#define NIFTI_ECODE_IGNORE           0  /* changed from UNKNOWN, 29 June 2005 */
#define NIFTI_ECODE_DICOM            2  /* intended for raw DICOM attributes  */
#define NIFTI_ECODE_AFNI             4  /* Robert W Cox: rwcox@nih.gov http://afni.nimh.nih.gov/afni */
#define NIFTI_ECODE_COMMENT          6  /* plain ASCII text only              */
#define NIFTI_ECODE_XCEDE            8  /* David B Keator: dbkeator@uci.edu http://www.nbirn.net/Resources/Users/Applications//xcede/index.htm */
#define NIFTI_ECODE_JIMDIMINFO      10  /* Mark A Horsfield:mah5@leicester.ac.uk http://someplace/something         */
#define NIFTI_ECODE_WORKFLOW_FWDS   12  /* Kate Fissell: fissell@pitt.edu http://kraepelin.wpic.pitt.edu/~fissell/NIFTI_ECODE_WORKFLOW_FWDS/NIFTI_ECODE_WORKFLOW_FWDS.html   */
#define NIFTI_ECODE_FREESURFER      14  /* http://surfer.nmr.mgh.harvard.edu  */
#define NIFTI_ECODE_PYPICKLE        16  /* embedded Python objects http://niftilib.sourceforge.net/pynifti                     */
/* LONI MiND codes: http://www.loni.ucla.edu/twiki/bin/view/Main/MiND */
#define NIFTI_ECODE_MIND_IDENT      18  /* Vishal Patel: vishal.patel@ucla.edu*/
#define NIFTI_ECODE_B_VALUE         20
#define NIFTI_ECODE_SPHERICAL_DIRECTION 22
#define NIFTI_ECODE_DT_COMPONENT    24
#define NIFTI_ECODE_SHC_DEGREEORDER 26  /* end LONI MiND codes                */
#define NIFTI_ECODE_VOXBO           28  /* Dan Kimberg: www.voxbo.org         */
#define NIFTI_ECODE_CARET           30  /* John Harwell: john@brainvis.wustl.edu http://brainvis.wustl.edu/wiki/index.php/Caret:Documentation:CaretExtension */
#define NIFTI_MAX_ECODE             30  /******* maximum extension code *******/
#define LNI_MAX_NIA_EXT_LEN 100000  /* consider a longer extension invalid */

#pragma mark Extension Class
class Extension {
	private:
		int m_code;          //!< Extension code, one of the NIFTI_ECODE_ values
		vector<char> m_data; //!< Raw data, with no byte swapping (length is esize-8)
	
	public:
		static const string &CodeName(const int code);
		
		Extension(int code, vector<char> data);
		Extension(int size, int code, char *data);
		const int rawSize() const;
		const int size() const;
		const int padding() const;
		const int code() const;
		const string &codeName() const;
		void setCode(int code);
		
		const vector<char> &data() const;
		void setData(const vector<char> &data);
};

} // End namespace Nifti

#endif // NIFTI_EXTENSION
