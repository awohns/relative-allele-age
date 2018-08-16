//
//  LoadSim.hpp
//  rvage
//
//  Created by Patrick Albers on 13/02/2017.
//  Copyright © 2017 pkalbers. All rights reserved.
//

#ifndef LoadSim_hpp
#define LoadSim_hpp

#include <stdio.h>
#include <exception>

#include "Reader.hpp"

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"

#include "IBD_SIM.hpp"


class LoadSim
{
public:
	
	// construct
	LoadSim(const std::string &);
	
	// forward to next line
	bool next();
	
	// parse current line
	bool parse(size_t &, IBD::SIM::Truth &);
	
private:
	
	Reader stream; // stream of GEN file
	
	bool good; // flag status, false = exit
};


#endif /* LoadSim_hpp */
