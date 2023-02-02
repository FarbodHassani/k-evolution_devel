#!/bin/bash

# Find all files with the "_cdm" in the name in subdirectories
for file in $(find . -name "*_cdm*" -type f)
do
  # Execute the gadget_hdr command on the current file
  ./gadget_hdr "${file}" "r" #gadget_hdr is the executable from gadget_hdr.cpp
done
