#ifndef SAM_FLAG_H_
#define SAM_FLAG_H_

bool IsMultipleSegment(unsigned int flag) {
  return flag & 0x1;
}

bool AllSegmentsMapped(unsigned int flag) {
  return flag & 0x2;
}

bool SegmentUnmapped(unsigned int flag) {
  return flag & 0x4;
}

bool NextUnmapped(unsigned int flag) {
  return flag & 0x8;
}

bool IsReverseComplement(unsigned int flag) {
  return flag & 0x10;
}

bool IsNextReverseComplement(unsigned int flag) {
  return flag & 0x20;
}

bool IsSegemntFirst(unsigned int flag) {
  return flag & 0x40;
}

bool IsSegmentLast(unsigned int flag) { 
  return  flag & 0x80;
}

bool IsSecondaryAlignment(unsigned int flag) {
  return flag & 0x100;
}

bool IsNotPassedQuality(unsigned int flag) {
  return flag & 0x200;
}

bool IsDuplicate(unsigned int flag) {
  return flag & 0x400;
}

#endif
