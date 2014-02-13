#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_
#include <string>
#include <iostream>
#include <vector>
#include "../NucConversion.h"

#include "../defs.h"

//
// Objects for representing alignments.
//

#include "../datastructures/alignment/Path.h"
#include "../datastructures/alignment/AlignmentBlock.h"
#include "../datastructures/alignment/AlignmentMap.h"
#include "../datastructures/alignment/Alignment.h"


//
// Utilities for manipulating alignments.
//

#include "alignment/AlignmentUtils.h"
#include "alignment/AlignmentPrinter.h"

//
// Methods for computing alignments.
//


#include "alignment/ScoreMatrices.h"
#include "alignment/SWAlign.h"
#include "alignment/KBandAlign.h"
#include "alignment/SDPAlign.h"

#endif
