#ifndef UTILS_PHRED_UTILS_H_
#define UTILS_PHRED_UTILS_H_

double InversePhred(double phred) {
	float num = 1;
	float denom = (1+pow(10,phred/10));
	return num/denom;
}

double Phred(double pvalue) {
  double v = log10(pvalue);
  return -10*v;
}

#endif
