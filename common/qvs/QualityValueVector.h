#ifndef QVS_QUALITY_VALUE_VECTOR_H_
#define QVS_QUALITY_VALUE_VECTOR_H_

#include "QualityValue.h"
#include "utils/ProtectedNew.h"

template<typename T_QV>
class QualityValueVector {
 public:
  T_QV   *data;
  QVScale qvScale;

  T_QV &operator[](unsigned int pos) {
    return data[pos];
  }

  QualityValueVector<T_QV>() {
    data = NULL;
    // Default to phred.
    qvScale = PHRED;
  }

  QualityProbability ToProbability(unsigned int pos) {
    return QualityValueToProbability(data[pos], qvScale);
  }
  T_QV ToPhred(unsigned int pos) {
    if (qvScale == PHRED) {
      return data[pos];
    }
    else {
      return PacBioQVToPhred(data[pos]);
    }
  }
  void Copy(const QualityValueVector<T_QV> &rhs, const DNALength length) {
    Free();
    if (rhs.Empty()) { 
      return;
    }
    Allocate(length);
    memcpy(data, rhs.data, length * sizeof(T_QV));
  }
    
  void Free() {
    if (data != NULL) {
      delete[] data;
      data = NULL;
    }
  }

  void Allocate(unsigned int length) {
    data = ProtectedNew<T_QV>(length);
  }

  bool Empty() const {
    return data == NULL;
  }
  
  void ShallowCopy(const QualityValueVector<T_QV> &ref, int pos = 0) {
    data = &ref.data[pos];
    qvScale = ref.qvScale;
  }

};

#endif
