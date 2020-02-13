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

#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>


namespace compress
{

class CompressorInterface
{
  public:
    std::unordered_map<std::string, std::string> compressorParameters;

  protected:
    std::string compressorName;		// internal compressor name
    std::stringstream log;			// logfile stream
    size_t cbytes;					// compressed stream size in bytes

  public:
    virtual void init() = 0;
    virtual int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n) = 0;
    virtual int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n) = 0;
    virtual void close() = 0;

    std::string getCompressorInfo();
    std::string getCompressorName(){ return compressorName; }
    std::string getLog() { return log.str(); }
    size_t getCompressedSize(){ return cbytes; }
    std::string getParamsInfo();
	void clearLog() { log.str(""); }
};


inline std::string CompressorInterface::getCompressorInfo()
{
    std::stringstream dataInfo;
    dataInfo << "\nCompressor: " << compressorName << std::endl;

    return dataInfo.str();
}


inline std::string CompressorInterface::getParamsInfo()
{
    std::string paramString = "";
    for (auto it=compressorParameters.begin(); it!=compressorParameters.end(); it++)
    {
        if (paramString != "")
            paramString += "_";
        paramString += (*it).first + (*it).second;
    }

    return paramString;  
}

} // compress namespace

#endif // #ifdef COMPRESSION