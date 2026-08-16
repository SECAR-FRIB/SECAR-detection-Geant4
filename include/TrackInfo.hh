#ifndef TrackInfo_h
#define TrackInfo_h 1

#include "G4VUserTrackInformation.hh"

class TrackInfo : public G4VUserTrackInformation {
public:
  TrackInfo() = default;
  ~TrackInfo() override = default;

  bool recordedTarget = false;
  bool recordedStrip  = false;

};

#endif