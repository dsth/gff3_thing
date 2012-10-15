#include <iostream>
#include "utils.h"
#include <openssl/md5.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "boost/filesystem.hpp"

using boost::filesystem::exists;
// using boost::filesystem3::exists;

#define SPC "<< \" \" <<"
#define DSN_PATTERN "The dsn string must be of the form: "

///// typedef type ~alias... e.g. in match_result : typedef sub_match<BidirectionalIterator> value_type;

using std::string;

typedef boost::smatch::value_type subm;

std::ostream& operator<<(std::ostream & os, const subm & c) {
//return os << "pink shit biyach!?!?!?!\n\n";
    std::string auto_tmp(c);
    return os << auto_tmp << " ";
}

DB_PARAMS::DB_PARAMS(std::string dsn) {
    // overloaded non-inline costructor will take dsn and generate object..
    // DBI:mysql:database=dsth_Gmor_Cap_Testing;host=mysql-eg-devel-3.ebi.ac.uk;port=4208,ensrw,scr1b3d3
    boost::regex reg_dsn(DSN_REGEX);
    boost::smatch match_obj;
    if (boost::regex_match(dsn,match_obj,reg_dsn)){

        //// for some reason part assigning with temporaries just scambles stuff?!? - iterating?!?
        //// - try iterating through?!? - either pre-assign in batch to strings directly or use stringstream?!?

        std::stringstream pinky(std::stringstream::in|std::stringstream::out);
        pinky << match_obj[1] << match_obj[2] << match_obj[3] << match_obj[4] << match_obj[5];
        std::string db, h, u, p;
        int P;
        pinky >> db >> h >> P >> u >> p;
        // really shouldn't do this?!? - i.e. cast away const?!? - but it's just a temporary so whatever - can't be const as would
        // need to be handled with initialiser list - if want it constant still just make private and give public getter...
        
        ///y while we can store as std::string and take advantage of implicit casting of char* to string then when we either need 
        ///y to make getters or exatract c string whenever using mysql etc.. - or just 
        
        /*
        dbname = (char*)db.c_str();
        host = (char*)h.c_str();
        this->user = (char*)u.c_str();
        //user = (char*)u.c_str();
        pass = (char*)p.c_str();
        port = P;
        */

        dbname_ = db;
        host_ = h;
        this->user_ = u;
        //user = (char*)u.c_str();
        pass_ = p;
        port_ = P;

        /*
        boost::smatch::value_type sit;
        boost::smatch::iterator bit;
        std::cout << "m1:dbname="<<match_obj[1]<<std::endl;
        // string stream then assign en masse and then store!?!
        dbname = (char*)std::string(match_obj[1]).c_str();
        host = (char*)std::string(match_obj[2]).c_str(); 
        this->port = atoi(std::string(match_obj[3]).c_str()); // should use boost for this?!?
        user = (char*)std::string(match_obj[4]).c_str(); 
        pass = (char*)std::string(match_obj[5]).c_str(); 
        std::cout << "m1:dbname="<<match_obj[1]<<std::endl;

        std::cout << "match[1]=" << typeid(dbname).name() << std::endl;
        std::cout << "dbname="<<dbname<<std::endl;
        */

    } else throw std::runtime_error(DSN_PATTERN DSN_REGEX);
}

// opts are either from cli or sqlite?!?
CAP_META::CAP_META(char* capdb) {

/*
    boost::regex reg_dsn(DSN_REGEX);
    boost::smatch match_obj;
        pinky << match_obj[1] << match_obj[2] << match_obj[3] << match_obj[4] << match_obj[5];
        std::string db, h, u, p;
        int P;
        pinky >> db >> h >> P >> u >> p;
        dbname_ = db;
        host_ = h;
        this->user_ = u;
        //user = (char*)u.c_str();
        pass_ = p;
        port_ = P;
    } else throw std::runtime_error("Problem getting config data!?!");
        */

        std::cout << "let's grab those values"<<std::endl; 

        std::cout << "PINKERTON"<<std::endl; 
}

void filepath(string& filename) {
    static unsigned char md5_result[MD5_DIGEST_LENGTH];
    MD5((unsigned char*)filename.c_str(), filename.size(), md5_result); 
    char ca[3];
    string buf = std::move(filename);
    filename.clear();
    for(int i=0; i <2; i++) { // for(i=0; i <MD5_DIGEST_LENGTH; i++) {
        sprintf(ca,"%02X/",md5_result[i]);
        filename += ca;
    }                
    filename += buf;
    return;
}    
                                                             
// void filedir(string& filename, const string& dir) {
void filedir(string& filename, string dir) { // dir taken by value so can do whatever with the copy...

    static unsigned char md5_result[MD5_DIGEST_LENGTH];
    MD5((unsigned char*)filename.c_str(), filename.size(), md5_result); 
    char ca[3];
    string buf = std::move(dir);
    buf += "/";
    dir.clear();
    for(int i=0; i <2; i++) { // for(i=0; i <MD5_DIGEST_LENGTH; i++) {
        sprintf(ca,"%02X/",md5_result[i]);
        dir += ca;
    }                
    buf += dir;
    // int status;
    // status = mkdir("/home/cnd/mod1", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // boost::filesystem::create_directories("/tmp/a/b/c");
    if (!exists(buf)) boost::filesystem::create_directories(buf);

    buf += filename;
    filename = std::move(buf);

    return;
}    

