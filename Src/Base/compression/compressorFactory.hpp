/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/
#pragma once

#ifdef COMPRESSION

#include "compressorInterface.hpp"

namespace compress
{

class CompressorFactory
{
  public:
	static CompressorInterface * createCompressor(std::string compressorName)
	{
	  #ifdef CBENCH_HAS_BLOSC
		if (compressorName == "BLOSC")
		  return new BLOSCCompressor();
	  #endif

	  #ifdef CBENCH_HAS_SZ
		if (compressorName == "SZ")
		  return new SZCompressor();
	  #endif

		  return NULL;
	}
};

} // compress namespace

#endif // #ifdef COMPRESSION