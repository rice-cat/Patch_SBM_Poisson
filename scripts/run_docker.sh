#!/bin/bash
# Use deal.II docker image
docker run --rm -t -v $(pwd):/builds/app dealii/dealii:v9.5.0-focal /bin/bash -c "cd /builds/app && rm -rf build && mkdir -p build && cd build && cmake .. && make"
