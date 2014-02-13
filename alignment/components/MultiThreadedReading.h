#ifndef COMPONENTS_ALIGNMENT_MULTI_THREADED_READING_H_
#define COMPONENTS_ALIGNMENT_MULTI_THREADED_READING_H_

#include "files/ReaderAgglomerate.h"
#include "SMRTSequence.h"
#include "Semaphores.h"
#include "datastructures/alignment/AlignmentContext.h"

bool GetNextReadThroughSemaphore(ReaderAgglomerate &reader, MappingParameters &params, SMRTSequence &read, AlignmentContext &context) {

  //
  // Grab the value of the semaphore for debugging purposes.
  //
  // uncomment when needed
  // int semvalue;
  // if (params.nProc > 1) {
  //  sem_getvalue(&semaphores.reader, &semvalue);
  // }

  //
  // Wait on a semaphore
  if (params.nProc > 1) {
    sem_wait(&semaphores.reader);
  }

  bool returnValue = true;
  //
  // CCS Reads are read differently from other reads.  Do static casting here
  // of this.
  //
  if (reader.GetFileType() == HDFCCS) {
    if (reader.GetNext((CCSSequence&)read) == 0) {
      returnValue = false;
    }
  }
  else {
    if (reader.GetNext((SMRTSequence&)read) == 0) { 
      returnValue = false;
    }
  }

  //
  // Set the read group id before releasing the semaphore, since other
  // threads may change the reader object to a new read group before
  // sending this alignment out to printing. 
  context.readGroupId = reader.readGroupId;
  
  if (params.nProc > 1) {
    sem_post(&semaphores.reader);
  }


  return returnValue;
}

#endif
