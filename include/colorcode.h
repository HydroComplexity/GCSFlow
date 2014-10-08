#pragma once

#define RS "\e[m"       // reset color

// L1 and L2 are used in Printf command for printing result
#define L1 "\x1b[31m"   // red
#define L2 "\x1b[33m"   // yellow

// Change the default color code by replacing L1 and L2 definition by 
//  the codes below if necessary.
#define L3 "\x1b[32m"   // green
#define L4 "\x1b[34m"   // blue
#define L5 "\x1b[35m"   // magenta
#define L6 "\x1b[36m"   // cyan
