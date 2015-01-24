#!/bin/bash

cp ~/GitHub/autoBlock/autoBlock_utils.R ~/GitHub/automated-blocking-examples/
cp ~/GitHub/autoBlock/modelfiles/*.RData ~/GitHub/automated-blocking-examples/modelfiles/
cp ~/GitHub/autoBlock/runscripts/* ~/GitHub/automated-blocking-examples/runscripts/

sed -i "" "s/GitHub\/autoBlock/GitHub\/automated-blocking-examples/" ~/GitHub/automated-blocking-examples/runscripts/*
