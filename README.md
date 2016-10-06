##Installation##

###1. Requirements###

  To build BLASR, you must have hdf 1.8.0 or above installed and
  configured with c++ support (you should have the library
  libhdf5_cpp.a).  If you are intalling the entire PacBio secondary
  analysis software suite, appropriate hdf libraries are already
  distributed and no configuration is necessary.  Otherwise, it is
  necessary to point two environment variables, HDF5INCLUDEDIR and
  HDF5LIBDIR to the locations of the HDF5 libraries.

  For example:

+    > export HDF5INCLUDEDIR=/usr/include/hdf
+    > export HDF5LIBDIR=/usr/lib/hdf

###2. Build the source tree###

+    > make


###3. The executable will be in alignment/bin/blasr###

+    > cd alignment/bin/blasr


##Running BLASR##

Typing blasr -h or blasr -help on the command line will give you a
list of options.  At the least, provide a fasta, fastq, or bas.h5 file,
and a genome.

*Some typical use cases:*

+    Align reads from reads.bas.h5 to ecoli_K12 genome, and output in SAM format.

		 blasr reads.bas.h5  ecoli_K12.fasta -sam

+    Same as above, but with soft clipping

		 blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft

+    Create sam output that may be used to resolve structural variation using local assembly

		 blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping subread -bestn 2 


+    Use multiple threads

		 blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft -out alignments.sam -nproc 16

+    Include a larger minimal match, for faster but less sensitive alignments

		 blasr reads.bas.h5  ecoli_K12.fasta -sam -clipping soft -minMatch 15

+    Produce alignments in a pairwise human readable format

		 blasr reads.bas.h5  ecoli_K12.fasta -m 0

+    Use a precomputed suffix array for faster startup

		 sawriter hg19.fasta.sa hg19.fasta #First precompute the suffix array
		 blasr reads.bas.h5 hg19.fasta -sa hg19.fasta.sa

+    Map assembled contigs (multiple megabases) to a reference

		 blasr human.ctg.fasta  hg19.fasta -alignContigs -sam -out alignments.sam

+    Use a precomputed BWT-FM index for smaller runtime memory footprint, but slower alignments.

		 sa2bwt hg19.fasta hg19.fasta.sa hg19.fasta.bwt
		 blasr reads.bas.h5 hg19.fasta -bwt hg19.fasta.bwt

## Output formats ##

The most universally compatible output is the SAM format, specified by ''-sam''. Other formats specified by the ''-m'' option conform to different applications, and as such the meanings of columns are not consistent between formats.   Alignments reported on the reverse strand may be converted to the forward strand using forward_start = length - reverse_end, reverse_start = length  - forward_start.  All output except for SAM is half-open zero based. 

+		 -m 0

A human readable version 

+		 -m 1

<table>
<tr> <td> 1.  </td> <td> query name </td> </tr>
<tr> <td> 2.  </td> <td> ref contig name </td> </tr>
<tr> <td> 3.  </td> <td> query strand </td> </tr>
<tr> <td> 4.  </td> <td> ref strand </td> </tr>
<tr> <td> 5.  </td> <td> align score </td> </tr>
<tr> <td> 6.  </td> <td> alignment percent identity </td> </tr>
<tr> <td> 7.  </td> <td> ref align start  </td> </tr>
<tr> <td> 8.  </td> <td> ref align end  </td> </tr>
<tr> <td> 9.  </td> <td> ref length </td> </tr>
<tr> <td> 10. </td> <td> query align start </td> </tr>
<tr> <td> 11. </td> <td> query align end </td> </tr>
<tr> <td> 12. </td> <td> query length </td> </tr>
<tr> <td> 13. </td> <td> alignment space usage </td> </tr>
</table>

Reverse strand alignments are reported starting at the 3' end of the reverse strand.

-m 2

XML based output. Reverse strand alignments are reported starting at the 3' end of the reverse strand.

+    -m 3
 
VULGAR alignment format from EXONERATE (deprecated)

+    -m 4

<table>
<tr> <td> 1.  </td> <td> query name </td> </tr>
<tr> <td> 2.  </td> <td> ref contig name </td> </tr>
<tr> <td> 3.  </td> <td> align score </td> </tr>
<tr> <td> 4.  </td> <td> alignment percent identity </td> </tr>
<tr> <td> 5.  </td> <td> query strand </td> </tr>
<tr> <td> 6. </td> <td> query align start </td> </tr>
<tr> <td> 7. </td> <td> query align end </td> </tr>
<tr> <td> 8. </td> <td> query length </td> </tr>
<tr> <td> 9.  </td> <td> ref strand </td> </tr>
<tr> <td> 10.  </td> <td> ref align start  </td> </tr>
<tr> <td> 11.  </td> <td> ref align end  </td> </tr>
<tr> <td> 12.  </td> <td> ref length </td> </tr>
<tr> <td> 13. </td> <td> alignment space usage </td> </tr>
<tr> <td> 14. </td> <td> mapping quality </td> </tr>
</table>
Reverse strand alignments are reported starting at the 3' end of the reverse strand.

+		 -m 5

This alignment format contains the full representation of the pairwise alignment of the two sequences in a verbose (easily parsed) stick format. 


<table>
<tr> <td> 1.  </td> <td> query name </td> </tr>
<tr> <td> 2.  </td> <td> query length </td> </tr>
<tr> <td> 3.  </td> <td> query align start </td> </tr>
<tr> <td> 4.  </td> <td> query align end </td> </tr>
<tr> <td> 5.  </td> <td> query strand </td> </tr>
<tr> <td> 6.  </td> <td> ref name </td> </tr>
<tr> <td> 7. </td> <td> ref length </td> </tr>
<tr> <td> 8. </td> <td> ref align start </td> </tr>
<tr> <td> 9. </td> <td> ref align end </td> </tr>
<tr> <td> 10.  </td> <td> ref strand </td> </tr>
<tr> <td> 11.  </td> <td> score  </td> </tr>
<tr> <td> 12.  </td> <td> nMatch  </td> </tr>
<tr> <td> 13.  </td> <td> nMismatch </td> </tr>
<tr> <td> 14. </td> <td> nIns </td> </tr>
<tr> <td> 15. </td> <td> nDel </td> </tr>
<tr> <td> 16. </td> <td> mapping quality </td> </tr>
<tr> <td> 17. </td> <td> ref align string </td> </tr>
<tr> <td> 18. </td> <td> query align string </td> </tr>
<tr> <td> 19. </td> <td> stick string </td> </tr>
<tr> <td> 20. </td> <td> ref align string </td> </tr>
</table>

For reverse strand alignments, the coordinates are reported starting at the 5' end of the forward strand. 

[![Build Status](https://travis-ci.org/zeeev/blasr.svg?branch=master)](https://travis-ci.org/zeeev/blasr.svg?branch=master)

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/036412483bfb92d2f18c1e62a34586b4 "githalytics.com")](http://githalytics.com/PacificBiosciences/blasr)
