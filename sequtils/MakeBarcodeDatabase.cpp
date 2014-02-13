#include "FASTASequence.h"
#include "FASTAReader.h"
#include "utils.h"

#include <string>
int main(int argc, char* argv[]) {
  string barcodeFileName, insertFileName, outputFileName;
  if (argc != 4) {
    cout << "usage: makeBarcodeDatabase insert.fasta barcodes.fasta output.fasta" << endl;
    exit(1);
  }
  insertFileName = argv[1];
  barcodeFileName = argv[2];
  outputFileName  = argv[3];

  FASTAReader barcodeReader, insertReader;
  barcodeReader.Initialize(barcodeFileName);
  insertReader.Initialize(insertFileName);
  
  ofstream barcodedOut;
  CrucialOpen(outputFileName, barcodedOut, std::ios::out);

  vector<FASTASequence> forwardBarcodes, reverseBarcodes;
  FASTASequence barcodeSequence, reverseBarcodeSequence;
  while(barcodeReader.GetNext(barcodeSequence)) {
    forwardBarcodes.push_back(barcodeSequence);
    barcodeSequence.MakeRC(reverseBarcodeSequence);
    reverseBarcodes.push_back(reverseBarcodeSequence);
  }
  
  FASTASequence insert;
  insertReader.GetNext(insert);
  
  int i;
  for (i = 0; i < forwardBarcodes.size(); i++) {
    FASTASequence barcodedInsert;
    barcodedInsert.Resize(forwardBarcodes[i].length * 2 + insert.length);
    stringstream titleStrm;
    titleStrm << insert.title << "|ff|" << forwardBarcodes[i].title;
    barcodedInsert.CopyTitle(titleStrm.str());
    memcpy(&barcodedInsert.seq[0], &forwardBarcodes[i].seq[0], forwardBarcodes[i].length);
    memcpy(&barcodedInsert.seq[forwardBarcodes[i].length], insert.seq, insert.length);
    memcpy(&barcodedInsert.seq[forwardBarcodes[i].length + insert.length], forwardBarcodes[i].seq, forwardBarcodes[i].length);
    barcodedInsert.PrintSeq(barcodedOut);

    titleStrm.str("");
    titleStrm << insert.title << "|fr|" << forwardBarcodes[i].title;
    barcodedInsert.CopyTitle(titleStrm.str());
    memcpy(&barcodedInsert.seq[0], &forwardBarcodes[i].seq[0], forwardBarcodes[i].length);
    memcpy(&barcodedInsert.seq[forwardBarcodes[i].length], insert.seq, insert.length);
    memcpy(&barcodedInsert.seq[forwardBarcodes[i].length + insert.length], reverseBarcodes[i].seq, reverseBarcodes[i].length);
    barcodedInsert.PrintSeq(barcodedOut);


    titleStrm.str("");
    titleStrm << insert.title << "|rf|" << forwardBarcodes[i].title;
    barcodedInsert.CopyTitle(titleStrm.str());
    memcpy(&barcodedInsert.seq[0], &reverseBarcodes[i].seq[0], reverseBarcodes[i].length);
    memcpy(&barcodedInsert.seq[reverseBarcodes[i].length], insert.seq, insert.length);
    memcpy(&barcodedInsert.seq[reverseBarcodes[i].length + insert.length], forwardBarcodes[i].seq, forwardBarcodes[i].length);
    barcodedInsert.PrintSeq(barcodedOut);


    titleStrm.str("");
    titleStrm << insert.title << "|rr|" << forwardBarcodes[i].title;
    barcodedInsert.CopyTitle(titleStrm.str());
    memcpy(&barcodedInsert.seq[0], &reverseBarcodes[i].seq[0], reverseBarcodes[i].length);
    memcpy(&barcodedInsert.seq[reverseBarcodes[i].length], insert.seq, insert.length);
    memcpy(&barcodedInsert.seq[reverseBarcodes[i].length + insert.length], reverseBarcodes[i].seq, reverseBarcodes[i].length);
    barcodedInsert.PrintSeq(barcodedOut);
  }
}
  
    

  

