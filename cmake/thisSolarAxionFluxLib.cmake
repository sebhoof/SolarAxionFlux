# Write thisSolarAxionFluxLib.[c]sh to INSTALL directory

install(CODE
"
file( WRITE ${PROJECT_SOURCE_DIR}/bin/thisSolarAxionFluxLib.sh

\"\#!/bin/bash

\# Check active shell by checking for existence of _VERSION variable
if [[ -n \\\"\\\${BASH_VERSION}\\\" ]]; then
    THIS_DIR=\\\$(cd \\\$(dirname \\\${BASH_ARGV[0]}); pwd)
elif [[ -n \\\"\\\${ZSH_VERSION}\\\" ]]; then
    THIS_DIR=\\\$(cd \\\$(dirname \\\$0); pwd)
else
    echo \\\"Invalid shell! Either source with bash or zsh!\\\"
    return 1
fi

if [ \\\$SOLAXFLUXLIB_PATH ] ; then
echo Switching to SolarAxionFlux library installed in \\\${THIS_DIR}
_PATH=`echo \\\$PATH | sed -e \\\"s\#\\\${SOLAXFLUXLIB_PATH}/bin:\#\#g\\\"`
_LD_LIBRARY_PATH=`echo \\\$LD_LIBRARY_PATH | sed -e \\\"s\#\\\${SOLAXFLUXLIB_PATH}/lib:\#\#g\\\"`
else
_PATH=\\\$PATH
_LD_LIBRARY_PATH=\\\$LD_LIBRARY_PATH
fi

export SOLAXFLUXLIB_PATH=\\\${THIS_DIR}
export SOLAXFLUXLIB_INCLUDE_PATH=\\\${THIS_DIR}/include
export SOLAXFLUXLIB_LIBRARY_PATH=\\\${THIS_DIR}/lib
export SOLAXFLUXLIB_DATA_PATH=\\\${THIS_DIR}/data

export PATH=\\\${SOLAXFLUXLIB_PATH}/bin:\\\$_PATH
export LD_LIBRARY_PATH=\\\${SOLAXFLUXLIB_LIBRARY_PATH}:\\\$_LD_LIBRARY_PATH
export LIBRARY_PATH=\\\$LIBRARY_PATH:\\\${SOLAXFLUXLIB_LIBRARY_PATH}

\"
)
"
)



install( CODE
"
file( WRITE ${PROJECT_SOURCE_DIR}/bin/thisSolarAxionFluxLib.csh

\"\#!/bin/csh

setenv SOLAXFLUXLIB_PATH ${PROJECT_SOURCE_DIR}
setenv SOLAXFLUXLIB_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/include
setenv SOLAXFLUXLIB_LIBRARY_PATH ${PROJECT_SOURCE_DIR}/lib
setenv SOLAXFLUXLIB_DATA_PATH ${PROJECT_SOURCE_DIR}/data

setenv PATH \\\${SOLAXFLUXLIB_PATH}/bin:\\\$PATH
setenv LD_LIBRARY_PATH \\\${SOLAXFLUXLIB_LIBRARY_PATH}:\\\$LD_LIBRARY_PATH
if ( \\\$?LIBRARY_PATH == 0 ) setenv LIBRARY_PATH
setenv LIBRARY_PATH \\\${LIBRARY_PATH}:\\\${SOLAXFLUXLIB_LIBRARY_PATH}

\"
)
"
)
