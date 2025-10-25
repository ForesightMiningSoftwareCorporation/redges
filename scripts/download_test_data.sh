#!/bin/bash
# Download test data for mesh repair testing

set -e

# Create directories
mkdir -p tests/fixtures/stanford
mkdir -p tests/fixtures/synthetic
mkdir -p tests/fixtures/defective

# Download Stanford Bunny (simplified version)
echo "Downloading Stanford Bunny..."
curl -L -o tests/fixtures/stanford/bunny.obj \
  https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/master/data/stanford-bunny.obj

# Download other common test models
echo "Downloading Armadillo..."
curl -L -o tests/fixtures/stanford/armadillo.obj \
  https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/master/data/armadillo.obj || echo "Armadillo not available"

echo "Downloading Cow..."
curl -L -o tests/fixtures/stanford/cow.obj \
  https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/master/data/cow.obj || echo "Cow not available"

echo "Download complete!"
