#ifndef GFFPROC_H
#define GFFPROC_H
#include "utils.h"
#include <string>
#include "log4cpp/Category.hh"
#include <utility>
#include "email.h"

class GFFPROC {
    
    ////// overload it for string?!? - always throw in constructor - but never in destructor?!?
    static const char summary_template[];
    std::string efingercmd;
    std::string gffdoccmd;

    log4cpp::Category & log4;
    DB_PARAMS* dbp;

public:

    GFFPROC (std::string _ss, std::string _s, log4cpp::Category & _l) 
      : gffdoccmd(_ss), efingercmd(_s), log4(_l) {}

    // int validate (char*, int, DB_PARAMS&, std::string&); // tmp extra message for admin?!?
    int validate (char*, int, DB_PARAMS&, std::string&, std::string&);

    bool gffcheck (const char*, std::string&);

    int gffload (const char*, std::string&);
    // int gffload (const char*);
    
    std::pair<int,int> eclean (const char*, int);

};

int int_from_mysql_query (DB_PARAMS *, const char * query); // sheer lazyness to get simple acc

#endif

