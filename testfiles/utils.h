#ifndef UTILS_H
#define UTILS_H

#define SCHEMA_20SEP2012 "CREATE TABLE `cap_files` (\
  `sbm_id` INTEGER PRIMARY KEY AUTOINCREMENT,\
  `submitter_name` varchar(40) NOT NULL,\
  `user` varchar(15) NOT NULL,\
  `user_email` varchar(40) NOT NULL,\
  `species` varchar(40) NOT NULL,\
  `file_name` varchar(255) NOT NULL,\
  `file_type` varchar(40) NOT NULL,\
  `file_md5` varchar(40) NOT NULL,\
  `file_desc` varchar(255) NOT NULL,\
  `file_size` int(11) NOT NULL,\
  `rts` timestamp DATETIME DEFAULT null,\
  `ts` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,\
  `sbm_status` varchar(40) DEFAULT NULL\
);\
CREATE TABLE `cap_meta` (       `key` varchar(40) PRIMARY KEY NOT NULL,       `value` varchar(255) NOT NULL );\
CREATE TABLE `cap_species` (\
  `species` varchar(40) PRIMARY KEY NOT NULL,\
  `dsn` varchar(255) NOT NULL,\
  `prefix` varchar(10) NULL,\
  `max_id` int(11) NULL\
);\
CREATE TABLE `cap_status` (\
  `status_id` INTEGER PRIMARY KEY AUTOINCREMENT,\
  `ts` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,\
  `status` int NOT NULL\
);"

/// v. bad practice... but... for extra-redirection e.g. get preprocessors own values need to stringfy...
#define RETRY_NAP 180
#define LAZY_SLEEP(x) "In retry loop, sleeping for " STRINGIFY(x) " seconds."
#define COOKIE "cookie.txt"
#define STRING_SIZE 100
#define LONG_STRING_SIZE 500
#define REALLY_LONG_STRING_SIZE 2500
#define MED_STRING_SIZE 750 // should put check on file name length?!?

#define COMMAND "loop"

/// always use cout with ANSI so that they never go into the log files!?!
#define ANSI_COFF "\033[0m"

#define ANSI_FUNKYB "\033[1;35m"
#define ANSI_FUNKY "\033[0;35m"
#define ANSI_YELLOWB "\033[1;31m"
#define ANSI_YELLOW "\033[0;31m"
#define ANSI_BLUEB "\033[1;36m" 
#define ANSI_BLUE "\033[0;36m" 
#define ANSI_WHITEB "\033[1;37m" 
#define ANSI_WHITE "\033[0;37m" 

//////////// wtf are these doing?!?

// ploncker!?! you're using strings without including it!?!
#include <string>
#include <exception>
#include <ostream>
#include "boost/regex.hpp" 
#include <sstream>
#include <signal.h>

#define DSN_REGEX "DBI:mysql:database=(\\w+?);host=([\\w\\.-]+?);port=(\\d+?),(\\w+?),(\\w+?)"

#define CMD_GFFDOC "GffDoc.pl "
#define CMD_EFINGER "efingerprint "

#define CMD_GFFDOC_OPTS " -type contig=ignore -type match=ignore -type match_part=ignore \
-type pseudogenic_tRNA=mRNA:pseudogenic_tRNA -type ncRNA=mRNA:ncRNA \
-type tRNA=mRNA:tRNA -type miRNA=mRNA:miRNA -type rRNA=mRNA:rRNA \
-mRNA_biotype_dominant \
-non_protein_coding_types pseudogenic_tRNA \
-coordsystem toplevel -non_coding_cds -non_protein_coding -leafnonunique "

// struct DB_PARAMS {
class DB_PARAMS {

    /* 
    * if you want them non-rewrittable but still not const (i.e. can only set with initialiser list
    * simply make them private and provide readonly accessor!?!
    const char * host;
    const char * user; 
    const char * pass; 
    const int port;
    // std::string dbname;
    const char * dbname;

    char * host;
    char * user; 
    char * pass; 
    char * dbname;

    */

    std::string host_;
    std::string user_; 
    std::string pass_; 
    std::string dbname_;
    unsigned int port_;
    DB_PARAMS (); // prevent default consctruction

public:

    //DB_PARAMS (const char * _h, const char * _u, const char * _pw, int _p, const char * _db )
    //DB_PARAMS (char * _h, char * _u, char * _pw, int _p, char * _db )
    DB_PARAMS (std::string _h, std::string _u, std::string _pw, int _p, std::string _db)
      : host_(_h), user_(_u), pass_(_pw), port_(_p), dbname_(_db) {}

    DB_PARAMS (std::string); // just pass by value and copy... 

    ///y just inline these
    const char * dbname() { return dbname_.c_str(); } // by value so can't be screwed with
    const char * host() { return this->host_.c_str(); }
    const char * user() { return user_.c_str(); }
    const char * pass() { return pass_.c_str(); }
    unsigned int port() { return port_; }

};

class CAP_META {

/*
if(opt::converter == NULL)  opt::converter  =   "parsexlsNfasta2gff3.pl";
if(opt::uploader == NULL)   opt::uploader   =   "gffCapProcessing_MultiThreaded.pl";
//// don't use this - put in the appropriate version of perl in #!... should really use env?!? if(opt::perl == NULL)     opt::perl        =   "/homes/dsth/dev/localperl/bin/perl";
if(opt::capdb == NULL)  opt::capdb  =     SQLITE_DB_NAME;
static int startemail = 0;
static int realemail = 0;
static int emaillist = 0;
static int das = 0;
static int iterations = 0;
static int max_retries = 10;
static int max_server_retries = 10;
static int max_json_retries = 10;
static int max_db_retries = 10;
static int port = 4208;
static char *host;
static char *user;
static char *capdb;
static char *pass;
static char *dbname;
static char *converter;
static char *execdir;
static char *gffcheck;
static char *tempdir;
static char *uploader;
static char *perl;
static char *bin;
std::vector<std::string> emailvector;
}

namespace retry {
static int retries = 0;
static int server_retries = 0;
static int json_retries = 0;
static int db_retries = 0;
}

if(opt::to_id == NULL)      opt::to_id      =   "dsth@ebi.ac.uk";
if(opt::from_addr == NULL)  opt::from_addr  =   EMAIL_CAP;
if(opt::xemails == NULL)  opt::xemails  =     EMAIL_LISTSERVE;
if(opt::smtp == NULL)       opt::smtp       =   "smtp.ebi.ac.uk";
static int restart = 900;
static int sleep = 300;
static int xcurl = 0;
static char *loglevel;
static char *sbmdir;
static char *smtp;
static char *url;
static char *from_addr;
static char *xemails;
static char *to_id;
if(opt::loglevel == NULL)    opt::loglevel    =   "info";
if(opt::url == NULL)        opt::url        =   "http://vectorbase-cap.ensemblgenomes.org/?q=uploaded";
if(opt::sbmdir == NULL)     opt::sbmdir     =   "submissions";
if(opt::bin == NULL)        opt::bin        =   "/net/isilon3/production/panda/ensemblgenomes/development/dsth/NewCap/backend/monitor"; if(opt::execdir == NULL)    opt::execdir    =   "/net/isilon3/production/panda/ensemblgenomes/development/dsth/NewCap/backend/monitor";
*/

    CAP_META (); // avoid silly typos...
    std::string local_bindir_;
    std::string local_sbmdir_;
    std::string email_to_; 
    std::string email_from_; 
    std::string email_xtras_; 
    std::string capmon_loglevel_; 
    std::string remote_jsonurl_; 
    int capmon_restart_;
    int capmon_sleep_;
    bool capmon_xcurl_;

public:

    CAP_META (char*);

    ///y just inline these - simply return the address of the c_str component?!?
    const char* local_bindir() { return local_bindir_.c_str(); } // by value so can't be screwed with
    const char* local_sbmdir() { return this->local_sbmdir_.c_str(); }
    const char* email_to() { return email_to_.c_str(); }
    const char* email_from() { return email_from_.c_str(); }
    const char* email_xtras() { return email_xtras_.c_str(); }
    const char* capmon_loglevel() { return capmon_loglevel_.c_str(); }
    const char* remote_jsonurl() { return remote_jsonurl_.c_str(); }
    unsigned int capmon_restart() { return capmon_restart_; }
    unsigned int capmon_sleep() { return capmon_sleep_; }
    unsigned int capmon_xcurl() { return capmon_xcurl_; }

};

struct COORDS {

    int start;
    int end;
    int strand; // not actually using this atm.?!?

    COORDS(int _a,int _b,int _c) : start(_a), end(_b), strand(_c) {} // inline this

private:
    
    COORDS(); // no need for public default ctor?!?


};

void filepath(std::string&);
void filedir(std::string&, std::string);
// void filedir(std::string&, const std::string&);
// void terminate (int param)  __attribute__((noreturn)); 
void terminate (int sig, siginfo_t *info, void *ptr)  __attribute__((noreturn));
void sig_action_function(int sig, siginfo_t *info, void *ptr);
// size_t my_curl_write(void * ptr, size_t size, size_t nmemb, std::stringstream & stream) {

#endif
