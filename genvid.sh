#!/usr/bin/env bash

ffmpeg -framerate 60 -pattern_type glob -i 'renders/*.png' -r 60 -c:v libx264 -pix_fmt yuv420p output.mp4
