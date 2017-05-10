#!/usr/bin/env bash

ffmpeg -framerate 100 -pattern_type glob -i 'build/*.png' -r 60 -c:v libx264 -pix_fmt yuv420p output.mp4
