#include <string>
#include <vector>
#include "DataFormats/Common/interface/Wrapper.h"

namespace 
{      
   struct dictionary 
     {	
	std::vector<std::vector<bool > > dummy2;
	edm::Wrapper< std::vector<std::vector<bool> > > dummy3;
     };   
}
