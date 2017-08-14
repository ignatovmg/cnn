#!/bin/bash

PDB=${1}
pdbprep.pl ${PDB}
pdbnmd.pl ${PDB} '?'

