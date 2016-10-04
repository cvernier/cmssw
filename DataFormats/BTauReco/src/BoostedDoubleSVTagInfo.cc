#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"

#include <cstring>

// ------------------------------------------------------------------------------------------
StoredBoostedDoubleSVTagInfo::StoredBoostedDoubleSVTagInfo() 
{
  memset(this, 0, sizeof(StoredBoostedDoubleSVTagInfo));
}

// ------------------------------------------------------------------------------------------
StoredBoostedDoubleSVTagInfo::~StoredBoostedDoubleSVTagInfo() 
{
}

// ------------------------------------------------------------------------------------------
BoostedDoubleSVTagInfo::BoostedDoubleSVTagInfo() 
{
  memset(this, 0, sizeof(BoostedDoubleSVTagInfo));
}

// ------------------------------------------------------------------------------------------
BoostedDoubleSVTagInfo::~BoostedDoubleSVTagInfo() 
{
}
