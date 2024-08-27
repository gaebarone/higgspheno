#!/bin/bash

# Run draw_eff_table.C with the specified arguments
root -b <<EOF
 EOF
.x draw_eff_table.C("wpwmhjj", "zzhjj", "wpwmjj", "zzjj",  "hwpwmjj", "hzzjj")
.q
EOF

# Run draw_normalized_dists.C with the specified arguments
root -b <<EOF
 EOF
.x draw_normalized_dists.C("wpwmhjj", "zzhjj", "wpwmjj", "zzjj", "hwpwmjj", "hzzjj")
.q
EOF

# Run draw_stacks.C for "wpwmhjj" with scaling factor 1
root -b <<EOF
 EOF
.x draw_stacks.C("wpwmhjj", 1)
.q
EOF

# Run draw_stacks.C for "wpwmjj" with scaling factor 1
root -b <<EOF
 EOF
.x draw_stacks.C("wpwmjj", 1)
.q
EOF

# Run draw_stacks.C for "zzhjj" with scaling factor 1
root -b <<EOF
 EOF
.x draw_stacks.C("zzhjj", 1)
.q
EOF

# Run draw_stacks.C for "zzjj" with scaling factor 1
root -b <<EOF
 EOF
.x draw_stacks.C("zzjj", 1)
.q
EOF


root -b <<EOF
 EOF
.x draw_stacks.C("hwpwmjj", 1)
.q
EOF


root -b <<EOF
 EOF
.x draw_stacks.C("hzzjj", 1)
.q
EOF
