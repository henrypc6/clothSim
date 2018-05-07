#!/bin/sh

# Cut videos
# ffmpeg -i out.mp4 -ss 00:00:01 -t 00:00:20 -strict -2 0.9.mp4

# Add text to Videos
# ffmpeg -i 0.9.mp4 -vf drawtext="text='COR = 0.9': fontcolor=white: fontsize=24: box=1: boxcolor=black@0.5: \
# boxborderw=5: x=(w-text_w)/2: y=100" -codec:a copy cor_0.9.mp4 -y

# Concatenate videos
ffmpeg -i "concat:cor_0.3.mp4|cor_0.6.mp4|cor_0.9.mp4" -c copy output.mp4