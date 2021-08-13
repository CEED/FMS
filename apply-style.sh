#!/bin/bash

# Copyright (c) 2021, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

function main()
{
  cd $(dirname "$0")
  fms_astyle_file="fms.astylerc"
  if [[ ! -r "${fms_astyle_file}" ]]; then
    echo "FMS's astyle format file not found: '${fms_astyle_file}'. Stop."
    exit 21
  fi

  find_astyle

  local old_IFS="${IFS}"
  IFS=$'\n'
  format_files=($(git ls-files "*.[ch]" "*.[ch]pp"))
  if [[ "$?" -ne 0 ]]; then
    echo "Error getting list of C/C++ source files from Git. Stop."
    exit 22
  fi
  IFS="${old_IFS}"

  if ${astyle_bin} --options="${fms_astyle_file}" "${format_files[@]}" | \
     grep "Formatted"; then
    printf "\nPlease make sure the changes are committed.\n\n"
    return 1
  else
    printf "All source files are properly formatted.\n"
  fi
  return 0
} # end of function 'main'

function find_astyle()
{
  astyle_req_version="Artistic Style Version 3.1"
  astyle_bin_list=("${ASTYLE_BIN:-astyle}" astyle-3.1)
  for astyle_bin in "${astyle_bin_list[@]}"; do
    if ! command -v "${astyle_bin}" > /dev/null 2>&1; then
      continue
    fi
    astyle_version="$("${astyle_bin}" --version)"
    if [[ "${astyle_version}" != "${astyle_req_version}" ]]; then
      continue
    fi
    return 0
  done
  echo "Required astyle version not found: '${astyle_req_version}'."
  printf "Astyle commands tried:"
  printf " '%s'" "${astyle_bin_list[@]}"
  printf ".\n"
  exit 23
} # end of function 'find_astyle'


# Invoke the 'main' function
main "$@"
