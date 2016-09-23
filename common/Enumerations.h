#ifndef ENUMERATIONS_H_
#define ENUMERATIONS_H_

typedef enum T_FileType {  Fasta,
													 Fastq,
													 HDFPulse,
													 Fourbit,
													 HDFBase,
													 HDFCCS,
													 HDFCCSONLY,
													 SAM_READ,
													 None
} FileType;
typedef enum T_Strand {Forward, Reverse} Strand;

typedef enum T_PlatformType {
	Springfield,
	Astro,
  NoPlatformType,
} PlatformType;

typedef enum T_RegionType {
    GlobalAccuracy,
    HQRegion,
    Adapter,
    Insert, 
    ArtifactRegion
} RegionType;

typedef enum T_PulseMetricType {
	QualityValueMetric,ClassifierQVMetric,StartTimeMetric,PulseWidthMetric,WidthInFramesMetric,
	pkmidMetric,IPDMetric,LightMetric, PreBaseFramesMetric
} PulseMetricType;


typedef enum T_AlignMode { NoAlignMode, Fullread, Subread, CCSDeNovo, CCSFullPass, CCSAllPass, ZmwSubreads } AlignMode;

#endif
