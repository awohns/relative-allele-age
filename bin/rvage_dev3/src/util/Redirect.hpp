//
//  RedirectLogErr.hpp
//  ship
//
//  Created by Patrick Albers on 17/02/2016.
//  Copyright © 2016 Patrick K. Albers. All rights reserved.
//

#ifndef Redirect_hpp
#define Redirect_hpp

#include <fstream>
#include <iostream>


//
// Redirect stream to output file stream
//
class Redirect
{
public:
	
	// constructs
	Redirect(std::ostream &, const std::string &);
	Redirect(std::ostream &, std::streambuf *);
	
	// destruct
	~Redirect();
	
	
private:
	
	std::ostream   & source;
	std::streambuf * buffer;
	std::ofstream    output;
};


#endif /* Redirect_hpp */

