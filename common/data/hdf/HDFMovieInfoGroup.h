#ifndef HDF_MOVIE_INFO_GROUP_H_
#define HDF_MOVIE_INFO_GROUP_H_

#include "HDFGroup.h"
#include "HDFArray.h"

class HDFMovieInfoGroup {
 public:
	
	HDFGroup movieInfoGroup;
	HDFArray<UInt> idArray;
	HDFArray<string>  nameArray;
	
	HDFArray<string> whenStartedArray;
	HDFArray<float> frameRateArray;
	HDFArray<string> bindingKitArray;
	HDFArray<string> sequencingKitArray;
	HDFArray<string> softwareVersionArray;

	~HDFMovieInfoGroup() {
		movieInfoGroup.Close();
	}

  bool Create(HDFGroup &parentGroup) {
    parentGroup.AddGroup("MovieInfo");
		if (movieInfoGroup.Initialize(parentGroup.group, "MovieInfo") == 0) { return 0; }
    idArray.Create(movieInfoGroup, "ID");
    nameArray.Create(movieInfoGroup, "Name");
		bindingKitArray.Create(movieInfoGroup, "BindingKit");
		sequencingKitArray.Create(movieInfoGroup, "SequencingKit");
		softwareVersionArray.Create(movieInfoGroup, "SoftwareVersion");
    return true;
  }

	int Initialize(HDFGroup &parentGroup) {
		if (movieInfoGroup.Initialize(parentGroup.group, "MovieInfo") == 0) { return 0; }
		if (idArray.Initialize(movieInfoGroup, "ID") == 0) { return 0; }
		if (nameArray.Initialize(movieInfoGroup, "Name") == 0) { return 0; }
		if (bindingKitArray.Initialize(movieInfoGroup, "BindingKit") == 0) {return 0;}
		if (sequencingKitArray.Initialize(movieInfoGroup, "SequencingKit") == 0) {return 0;}
		if (softwareVersionArray.Initialize(movieInfoGroup, "SoftwareVersion") == 0) { return 0;}
		return 1;
	}
	
	void Read(MovieInfo &movieInfo) {
		int nId = idArray.arrayLength;
		movieInfo.id.resize(nId);
		idArray.Read(0, nId, &movieInfo.id[0]);
		
		int nName = nameArray.arrayLength;
		movieInfo.name.resize(nName);
		movieInfo.sequencingKit.resize(nName);
		movieInfo.bindingKit.resize(nName);
		movieInfo.softwareVersion.resize(nName);
		int i;
		for (i = 0; i < nName; i++ ){
			nameArray.Read(i,i+1,&movieInfo.name[i]);
			sequencingKitArray.Read(i,i+1,&movieInfo.sequencingKit[i]);
			bindingKitArray.Read(i,i+1,&movieInfo.bindingKit[i]);
			softwareVersionArray.Read(i,i+1,&movieInfo.softwareVersion[i]);
		}
	}

  int AddMovie(string &movieName) {
    nameArray.Write(&movieName, 1);
    unsigned int id = nameArray.size();
    idArray.Write(&id, 1);
    return id;
  }


  int AddMovie(string &movieName, string &sequencingKit, string &bindingKit, string &softwareVersion) {
    nameArray.Write(&movieName, 1);
		sequencingKitArray.Write(&sequencingKit, 1);
		bindingKitArray.Write(&bindingKit, 1);
		softwareVersionArray.Write(&softwareVersion, 1);
    unsigned int id = nameArray.size();
    idArray.Write(&id, 1);
    return id;
  }

  void StoreFrameRate(int movieIndex, float frameRate) {
    if (movieIndex < 0) {
      cout << "ERROR. Invalid movie index " << movieIndex << endl;
      exit(1);
    }

    if (!frameRateArray.IsInitialized()) {
      if (!movieInfoGroup.ContainsObject("FrameRate")) {
        frameRateArray.Create(movieInfoGroup, "FrameRate");
      }
      else {
        frameRateArray.Initialize(movieInfoGroup, "FrameRate");
      }
    }
    frameRateArray.WriteToPos(&frameRate, 1, movieIndex);
  }
  
};
			
			
			
				
#endif
