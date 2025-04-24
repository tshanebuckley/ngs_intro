#!/bin/bash
# run this script to reset the training repo without the solutions completed.

USER="nextflow-io"
REPO="training"
FOLDER="hello-nextflow"

cd hello_world
rm -rf hello_nextflow
mkdir hello_nextflow
rm -rf training
git clone --depth 1 --filter=blob:none --sparse https://github.com/${USER}/${REPO}.git
cd "$REPO"
git sparse-checkout set "$FOLDER"

cp -a hello-nextflow/. ../hello_nextflow
cd ..
rm -rf training
