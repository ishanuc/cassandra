#!/bin/bash

parallel --gnu  -j 8 ./geteps_.sh ::: `ls -d ../_D__*[0-9]`

