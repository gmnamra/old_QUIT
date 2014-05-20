#ifndef QUIT_H
#define QUIT_H

#include <string>

//******************************************************************************
#pragma mark Versioning etc.
//******************************************************************************
#include "version"

const std::string credit_shared {
"Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const std::string credit_me {
"Written by tobias.wood@kcl.ac.uk. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

//******************************************************************************
#pragma Rest of the library
//******************************************************************************
#include "QUIT/MultiArray.h"
#include "QUIT/Volume.h"
#include "QUIT/VolumeFunctions.h"
#include "QUIT/ThreadPool.h"
#include "QUIT/Util.h"

#endif // QUIT_H
