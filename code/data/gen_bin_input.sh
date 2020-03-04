#!/usr/bin/env bash

SIZE=$1
STEPS=$2
THETA=$3
futhark dataset --seed=42 --binary \
        -g ${STEPS}i32 \
        -g 50f32 \
        -g 0.5f32 \
        -g ${THETA}f32 \
        --f32-bounds=0:1000 -g [$SIZE]f32 \
        -g [$SIZE]f32 \
        -g [$SIZE]f32 \
        -g [$SIZE]f32 \
        --f32-bounds=30:50 -g [$SIZE]f32 \
        --f32-bounds="(-10)":10 -g [$SIZE]f32 \
        -g [$SIZE]f32 \
        -g [$SIZE]f32
