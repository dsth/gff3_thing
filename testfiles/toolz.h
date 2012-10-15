#ifndef TOOLZ_H
#define TOOLZ_H
#include <curl/curl.h>
#include <stdio.h>
#include <typeinfo>
#include <string>
#include <sstream>
#include <algorithm>
#include <sqlite3.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <exception>
#include <type_traits>
#include <sys/stat.h> /* struct stat, fchmod (), stat (), S_ISREG, S_ISDIR */
#include "utils.h"

#define STRING_LENGTH 250
#define SHORT_STRING 100
#define MAGIC_NUMBER 2500 

/// perhaps start using boost::filesystem?!?
/// mysql?!?
/// sqlite insert?!?

/// dsn_by_species?!?
/// insert_speces?!? 
/// ...

/// insert_meta_key/value
/// update_meta_key/value


typedef std::vector<std::string> ROW;
typedef std::vector<std::vector<std::string>> TABLE;

inline size_t my_curl_write(void * ptr, size_t size, size_t nmemb, std::stringstream & stream) {
    std::string blah((char*)ptr);
    stream << blah;//(char*)ptr; // if you deref this you'll get just the one char!?!
    //// must return length or truncates... size vs. length?!?
    return blah.size();
} 

inline void set_cap_status (sqlite3 * db, int status) {
    char buf[STRING_LENGTH];
    sprintf(buf,SQL_CAP_STATUS,status);
    char * zErrMsg = NULL;
    if(sqlite3_exec(db, buf, NULL, NULL, &zErrMsg)!=SQLITE_OK)
      throw std::runtime_error(zErrMsg);
    sqlite3_free(zErrMsg);
}

// to stop continued connections can (1) new/delete the resource - i.e. instead of using auto var
// can use connect/disconnect with static accessor as follows - forces singleton!?!?!?
// give accessor for handling object - can use it to invoke overloaed member function query?!?
// USE template member function!?!? - i.e. can simply assign for appropriate return type?!?

/* ITS A HEADER!?!?!
using std::runtime_error;
using std::true_type;
using std::cout;
using std::endl;
using std::true_type;
using std::false_type;
*/

namespace toolz {
    
// should use boost/filesystem.hpp but want something entirely header driven atm?!?
inline bool file_exists(const char* dbname) {
    struct stat sts;
    ///// cannot return directly as presumably resource isn't properly free'd?!? - i.e. will get false on next check?!?
    bool x = (stat(dbname, &sts) == -1 && errno == ENOENT);
    return !x;
}

////// probably best way would be to have friend function we call on object to grab pointer to meta - thus can populate it.
////// then we have a grab meta member function and/or cache it - either way need to do if(meta.find("blah")==meta.end()) to make sure
////// it exists as [] will generate it?!?
///r qw / qe qr
/////// since the meta stuff is being taken out of checksqlite_forloop_populatemetamap also take out the check capdb status out - i.e.
/////// that should be a member function of sqlite thingy here?!?
/////// then rename it check_perl_conf?!? - and activate the checks on the major perl module requirements!??!

// inline static int internal_meta_callback(void *NotUsed, int argc, char **argv, char **azColName){
inline static int internal_meta_callback(void *blah, int argc, char **argv, char **azColName){
    
    ///r since callback is not templated - i.e. accepts structure to populate using C esque void* need to cast back type character!?!?
    // std::map<std::string,std::string>* _meta = (std::map<std::string,std::string>*)blah;
    std::map<std::string,std::string>* _meta = static_cast<std::map<std::string,std::string>*>(blah);
    _meta->insert(std::pair<std::string,std::string>(argv[0],argv[1])); // _meta[argv[0]]=argv[1];
    // std::cout << "key=" << argv[0] << "\tvalue="<<argv[1]<<std::endl; 
    return 0;
}

/* class ADAPTOR {
    virtual ~ADAPTOR = 0;
}; */

//// now make the meta_value routine check _meta first. then try the db - then throw?!? - i.e. cached values...

///r sqlite and mysql adaptor should be made into explicit instantiations of an un-defined but declared primary template?!?

//// then can put in the perl_conf stuff?!?
class MYSQL_ADAPTOR {

  private:

    static MYSQL_ADAPTOR* me_ptr;
    MYSQL* conn;
    DB_PARAMS& dp;

    MYSQL_ADAPTOR(DB_PARAMS& _dp) : dp(_dp) { 
        // cout << "\n\nopening connection\n\n";
        conn = mysql_init(NULL);
        if(mysql_real_connect(conn,dp.host(),dp.user(),dp.pass(),dp.dbname(),dp.port(),NULL,0) == NULL)
          throw std::runtime_error("could not connect to species-specific e! cap database");
    }

    ~MYSQL_ADAPTOR() {
        // cout << "\n\nclosing connection\n\n";
        mysql_close(conn);
    }

    MYSQL_ADAPTOR(const MYSQL_ADAPTOR&) = delete;

  public:


    static MYSQL_ADAPTOR* get_instance() {
        if(!me_ptr) throw std::runtime_error("pointer to sqlite object is not initialised");
        return me_ptr;
    }
    
    static void connect(DB_PARAMS& dp) {
    // std::cout << "\ngonna new it up\n";
        if(!me_ptr) me_ptr = new MYSQL_ADAPTOR(dp);
        // return me_ptr;
    }

    static void disconnect() {
    // std::cout << "\ngonna delete it\n";
        if(me_ptr) delete me_ptr;
        me_ptr = NULL; // return me_ptr;
    }

    //// really must run this stuff through valgrind?!?

    void func_helper2 (std::string q, std::string& x) { 

        MYSQL_RES *result;
        MYSQL_ROW row;
        MYSQL_FIELD *field;

        // cout <<"HEY -\n gonna \n   grab stuff\n";

        if (mysql_query(conn,q.c_str()))
          throw std::runtime_error("unable to query mysql instance");

        if(!(result = mysql_store_result(conn)))
          throw std::runtime_error("result set pointer is null");
        else {
            row = mysql_fetch_row(result); // if just selecting then this would indicate no entry - but not v. safe better to count 
            if (row == 0) throw std::runtime_error("result is null");
            x=row[0];
        }
        mysql_free_result(result);
    }

    template <typename T> T generic_query(std::string y) {
        T x;
        func_helper2(y, x); 
        return std::move(x);
    }

};

/// actually start to use RAII as you really ought to pretty much always do?!?
class SQLITE_ADAPTOR { // possibly shouldn't have quite all the member functions inlined - i.e. within class definition

  private:

    //// the callback is inlined - should all be accesible?!?
    std::map<std::string,std::string> _meta;

    static SQLITE_ADAPTOR* me_ptr;
    sqlite3* database;
    char* db_name;
    sqlite3_stmt* statement;

    SQLITE_ADAPTOR(char* dbname) : db_name(dbname) { 
    // SQLITE_ADAPTOR(const char* dbname) { 
    
        if(!file_exists(dbname)) throw std::runtime_error("the db doesn't exist " + std::string(db_name));
        // if(!file_exists(dbname)) throw runtime_error("the db doesn't exist " + std::string(dbname));

// std::cout << "\nctor! - just opened connection\n";
        sqlite3_open(dbname, &database); 

    }

    ~SQLITE_ADAPTOR() { sqlite3_close(database); /* std::cout << "\ndtor! - just closed connection\n"; */ }
    SQLITE_ADAPTOR(const SQLITE_ADAPTOR&) = delete;

  public:



    // static SQLITE_ADAPTOR* connect(char* dbname) {
    //     if(!me_ptr) me_ptr = new SQLITE_ADAPTOR(dbname);
    //     return me_ptr;
    // }
    //
    static SQLITE_ADAPTOR* get_instance() {
        if(!me_ptr) throw std::runtime_error("pointer to sqlite object is not initialised");
        return me_ptr;
    }
    
    static void connect(char* dbname) {
    // std::cout << "\ngonna new it up\n";
        if(!me_ptr) me_ptr = new SQLITE_ADAPTOR(dbname);
        // return me_ptr;
    }

    static void disconnect() {
    // std::cout << "\ngonna delete it\n";
        if(me_ptr) delete me_ptr;
        me_ptr = NULL; // return me_ptr;
    }

    /*
    template <typename T> T helper(T,true_type) {
        std::cout << "\nintegral form!\n";
    }

    ///// need helper function indirection?!? - will throw it away?!?
    template<typename T> T pullme() { 
        helper(typename std::is_integral<T>::type());
    }
    // defaults to string...

    // template<typename T> T pullme(char*); // not gonna actually define the primary template!?!?
    */

    void pull_files(const char*, int&) {
        std::cout << "this is the single int form\n";
    }

    void pull_files(const char*, std::string&) {
        std::cout << "this is the single string form\n";
    }

//==============================================================================
void func_helper2 (std::string q, double& x); // { // full specialisation for if-true of main_func
    // std::cout << " - double="<<x <<" and " << q << "\n";
//}

void func_helper2 (std::string q, std::string& x) { // full specialisation for if-true of main_func
    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
        std::cout << "using "<<q<<"\n";
        if(sqlite3_step(statement)== SQLITE_ROW) x=reinterpret_cast<const char*>(sqlite3_column_text(statement, 0));
        else throw std::runtime_error("no rows!?!");
    } else {
        sqlite3_finalize(statement);
        throw std::runtime_error("error in response");
    }
    sqlite3_finalize(statement);
}
void func_helper2 (std::string q, int& x) { // full specialisation for if-true of main_func
    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {

        // std::cout << "using "<<q<<"\n";
        if(sqlite3_step(statement)== SQLITE_ROW) {
            x=atoi(reinterpret_cast<const char*>(sqlite3_column_text(statement, 0)));
        } else throw std::runtime_error("no rows!?!");
            
    } else {
        sqlite3_finalize(statement);
        throw std::runtime_error("error in response");
    }
    sqlite3_finalize(statement);
}

void func_helper2 (std::string q, ROW& x) { 

    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
        int cols = sqlite3_column_count(statement);
        if(sqlite3_step(statement)== SQLITE_ROW) 
          for(int col = 0; col < cols; col++) 
            x.push_back(reinterpret_cast<const char*>(sqlite3_column_text(statement, col)));
        else throw std::runtime_error("no rows!?!");
        sqlite3_finalize(statement);
    } else {
        sqlite3_finalize(statement);
        throw std::runtime_error("error in response");
    }
    sqlite3_finalize(statement);
    // return std::move(x);
}

void func_helper2 (std::string q, TABLE& x) { // full specialisation for if-true of main_func
// void func_helper2 (std::string q, std::vector<std::vector<std::string>>& x) { // full specialisation for if-true of main_func
    if(sqlite3_prepare_v2(database, q.c_str(), -1, &statement, 0) == SQLITE_OK) {
            int cols = sqlite3_column_count(statement);
            int result = 0;
            int row =0;
            while(true) {
                result = sqlite3_step(statement);
                if(result == SQLITE_ROW) {
                    row++;
                    std::vector<std::string> s;
                    for(int col = 0; col < cols; col++) 
                      s.push_back(reinterpret_cast<const char*>(sqlite3_column_text(statement, col)));
                    x.push_back(std::move(s));
                } else break;  
            }
        }
    sqlite3_finalize(statement);
        /* for(int i=0;i<4;i++) {

    std::vector<std::string> s;

      for (int y=0;y<3;y++) 
        s.push_back(std::to_string((long long)i));
        x.push_back(std::move(s));
    }
    */
}

template <typename T> T generic_query(std::string y) {
    T x;
    // to make this work we NEED an actual object and NOT a typename - hence the final paren' 'std::is_integral<T>::type()'!?!?
    // std::cout << "THIS IS IS : " << y << " ";
    func_helper2(y, x); 
    return std::move(x);
}
//==============================================================================









template <typename T> T func_helper (T t, std::true_type) { // full specialisation for if-true of main_func
    std::cout << t << " : true" << std::endl;
}
 
template <typename T> T func_helper (T t, std::false_type) { // if-false
    std::cout << t << " : false" << std::endl;
}
 
template <typename T> T main_func(std::string y) { // the if part?!?
    T x;
    std::cout << "here we are = "<< y<<"\n";
    x = func_helper(x, typename std::is_integral<T>::type()); // need to instantiate a temporary object hence paren-paren
    return x;
}








    // void pull_files(std::vector<std::string>&);
    // just force inline by definition within class!?!? inline void SQLITE_ADAPTOR::pull_files(std::vector<std::string>& v) {
    void pull_files(const char*, std::vector<std::string>&) {


    if(sqlite3_prepare_v2(database, "select * from cap_files", -1, &statement, 0) == SQLITE_OK) {
            int cols = sqlite3_column_count(statement);
            int result = 0;
            int row =0;
            while(true) {
                result = sqlite3_step(statement);
                
                if(result == SQLITE_ROW) {
                    row++;
    //                std::cout << "ROW:" << row << std::endl;

                    for(int col = 0; col < cols; col++) {
                        std::cout << "  [" << col<<"] " << sqlite3_column_text(statement, col);
                        // can thus cast to whatever!?!
                        std::string s = (char*)sqlite3_column_text(statement, col);
                        //do something with it
                    }


                    std::cout << "\n";

                } else break;  

break;

            }
            
            sqlite3_finalize(statement);
        }
    }

    unsigned int check_capdb_status () {

        // int cap_status = query_integer(db,const_cast<char*>(SQL_CHECK_DB));
        if (!me_ptr) throw std::runtime_error("you have to connect first!");
        int cap_status = me_ptr->generic_query<int>(SQL_CHECK_DB);

        if (cap_status==0) {
            std::cout << "Cap_status table appears to be empty\nexiting" << std::endl;
            exit(0);
        } else if (cap_status>0) {

           //  cap_status = query_integer(db,const_cast<char*>(SQL_QUERY_STATUS));
            cap_status = me_ptr->generic_query<int>(SQL_QUERY_STATUS);

        } else throw Sqlite3Error("hmm");

        switch (cap_status) { //EnumParser<current_status> parser;
        case CAP_STATUS_READY:
            std::cout << "Capdb has status 'ready'\n";
            std::cout << "Setting status to 'looping'\n";
///// get rid of this!?!?!?
            set_cap_status(database,CAP_STATUS_LOOPING);
            break;
        case CAP_STATUS_LOOPING:
            std::cout << "The database wasn't not shut down shut down properly - run the ready command" << std::endl;
            exit(0);
            break;
        case CAP_STATUS_LOCKED:
            std::cout << "Please reset database state with ready command" << std::endl;
            exit(0);
            break;
        default:
            std::cout << "Unknown status" << std::endl;
            exit(1);
        }

    }

    void update(const char* q) {

        char* zErrMsg = NULL;
        if(sqlite3_exec(database, q, NULL, NULL, &zErrMsg)!=SQLITE_OK) // throw std::runtime_error(zErrMsg);
          throw std::runtime_error(THROW_DEFAULT("unable to update sqlite3 db"));
        sqlite3_free(zErrMsg);

    }

    void meta_insert(const char* k,const char* v) { // it's really not appropriate to handle problems with normal flow - i.e. bool in this context - just throw?!?
        char query[SHORT_STRING];
        sprintf(query,"insert into cap_meta values ('%s','%s');",k,v);
        char* zErrMsg = NULL;
        if(sqlite3_exec(database, query, NULL, NULL, &zErrMsg)!=SQLITE_OK) // throw std::runtime_error(zErrMsg);
          throw Sqlite3Error(THROW_DEFAULT("unable to insert meta key/value pair into capdb"));
        sqlite3_free(zErrMsg);
    }


/*  prepare/step for insert vs. exec?!?
 *  if(sqlite3_prepare_v2(sqlitedb,sqlstring.c_str(),-1,&stmt,0)!=SQLITE_OK) 
    if(sqlite3_step(stmt)!=SQLITE_DONE)
*/

    void generic_insert(const std::string& q) { // it's really not appropriate to handle problems with normal flow - i.e. bool in this context - just throw?!?
        char* zErrMsg = NULL;
        if(sqlite3_exec(database, q.c_str(), NULL, NULL, &zErrMsg)!=SQLITE_OK) // throw std::runtime_error(zErrMsg);
          throw Sqlite3Error(THROW_DEFAULT("unable to insert into capdb"));
        sqlite3_free(zErrMsg);
    }

    std::string meta_value(const char* key) {

        // clearly if _meta is populated no need to check db
        if(_meta.find(key)!=_meta.end()) // safe to use [] without generating entry
          return _meta[key];

        char q[SHORT_STRING];
        sprintf(q,"select value from cap_meta where key = '%s'",key);
        // std::cout << q << std::endl;

        std::string local;
        if(sqlite3_prepare_v2(database, q, -1, &statement, 0) == SQLITE_OK) {

            if(sqlite3_step(statement)== SQLITE_ROW) {
                local=reinterpret_cast<const char*>(sqlite3_column_text(statement, 0));
            } else {
                // sqlite3_finalize(statement); // ?!?
                throw std::runtime_error("Cannot find entry for " + std::string(key)); // just use sprintf?!?
            }
                
        } else {
            sqlite3_finalize(statement);
            throw std::runtime_error("error in response");
        }

        // if(local=="local.perl5lib") std::for_each(perl5lib.begin(),perl5lib.end(),[](char& n){if(n==',')n=':';});

        sqlite3_finalize(statement);
        return std::move(local);

    }

    void meta_populate() {

        char *zErrMsg = 0; 
        int rc = sqlite3_exec(database, "select * from cap_meta", internal_meta_callback, &_meta, &zErrMsg);
        // int rc = sqlite3_exec(database, "select * from cap_meta", internal_meta_callback, 0, &zErrMsg);
        
        /* std::function<int(void*,int,char**,char**)> make_offseter(void *NotUsed, int argc, char **argv, char **azColName)
        {
        return        [_meta&](void *NotUsed, int argc, char **argv, char **azColName)->int{
    // _meta[argv[0]]=argv[1];
    std::cout << "key=" << argv[0] << "\tvalue="<<argv[1]<<std::endl; 
    return 0;}
        } */
        // for (auto i = _meta.begin(); i!=_meta.end();i++)
          //  cout << "YAY key = "<< i->first << " and value = "<< i->second << "\n";
    // std::cout << "key=" << argv[0] << "\tvalue="<<argv[1]<<std::endl; 
        
        if(rc!=SQLITE_OK)
        fprintf(stderr, "SQL error: %s\n", zErrMsg);

    }
/*
 *



 */ 



};

class CURL_ADAPTOR {

    CURL* curl;
    bool xcurl;
    std::string cookie_file; // just let underlying string be statically allocated?!?

    static CURL_ADAPTOR* me_ptr;

    CURL_ADAPTOR(bool xc,std::string ck) : xcurl(xc), cookie_file(ck) { 
        if(!xcurl) {
            // duh?!? you've defined a stack version so the proper handle never gets initialise?!?!?!? CURL * curl;  
            curl_global_init(CURL_GLOBAL_ALL);
            curl = curl_easy_init();
        }
    }

    ~CURL_ADAPTOR() { 
        curl_easy_cleanup(curl); 
        if(me_ptr) delete me_ptr;
        me_ptr = NULL; 
    }

    CURL_ADAPTOR(const CURL_ADAPTOR&) = delete;

  public:

    static CURL_ADAPTOR* get_instance(bool xterm = false, std::string ck = COOKIE) {
        if(!me_ptr) me_ptr = new CURL_ADAPTOR(xterm, ck);
        return me_ptr;
    }

    void pull(std::string uri, std::stringstream& strstrm) {
        if (!xcurl) {
            curl_easy_setopt(this->curl,    CURLOPT_URL,                uri.c_str()); 
            curl_easy_setopt(curl,          CURLOPT_COOKIEFILE,         cookie_file.c_str()); // COOKIE);
            curl_easy_setopt(curl,          CURLOPT_COOKIEJAR,          cookie_file.c_str()); // COOKIE);
            curl_easy_setopt(curl,          CURLOPT_NOPROGRESS,         1);         
            curl_easy_setopt(curl,          CURLOPT_WRITEFUNCTION,      my_curl_write); 
            curl_easy_setopt(curl,          CURLOPT_WRITEDATA,          &strstrm);      
            
            curl_easy_perform(curl);
        } else {
            char tmp[MAGIC_NUMBER];
            std::string cmd("curl -s -c " + cookie_file + " -b " + cookie_file);
            cmd += uri;
            FILE * f2 = popen(cmd.c_str(), "r");
            while(( fgets( tmp, MAGIC_NUMBER, f2 )) != NULL ) strstrm << tmp;
            pclose(f2);
        }
    }

};






/// multiple inclusion issue - have no compiled stuff here?!? - i.e. make inline?!?!
inline int query_integer (char * query) {

    int db_return(0);
    sqlite3 * database;

    db_return = sqlite3_open("new_sqlite.db", &database);
    sqlite3_stmt * statement;

    if(sqlite3_prepare_v2(database, query, -1, &statement, 0) == SQLITE_OK) {

        if(sqlite3_step(statement)== SQLITE_ROW) {
            db_return = atoi((char*)sqlite3_column_text(statement, 0));
        } else throw std::runtime_error("no rows!?!");
            
    } else throw std::runtime_error("error in response");

    sqlite3_finalize(statement);

   // sqlite3_close(database);
    return db_return;
}

}






// static so much define it!?!?
// goes in main?!? toolz::SQLITE_ADAPTOR* toolz::SQLITE_ADAPTOR::me_ptr = 0; // seems std::nullptr not yet implemented?!?

#endif
