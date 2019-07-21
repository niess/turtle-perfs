#!/bin/bash
basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

bindir=${basedir}/bin
[[ "${PATH}" =~ "${bindir}" ]] || export PATH=${bindir}:$PATH

bindir=${basedir}/scripts
[[ "${PATH}" =~ "${bindir}" ]] || export PATH=${bindir}:$PATH

libdir=${basedir}/lib
[[ "${LD_LIBRARY_PATH}" =~ "${libdir}" ]] || export LD_LIBRARY_PATH=${libdir}:$LD_LIBRARY_PATH
