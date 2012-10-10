#include <fstream>
#include <iostream>
#include <sstream>
#include "boost/regex.hpp"
#include <set>
#include <algorithm>
#include "exceptions.h"
// #include "boost/lambda/lambda.hpp"

#include <iomanip> // manipulators
#include "utils.h"
#include "feature.h"

#ifdef CAPMON_EXT
#include <mysql/mysql.h> 
#endif

#define MAX_ID_LENGTH 50


///r whether or not use this for multiple situations by allowing to run simply externally from cap will make doing/extending tests simpler?!?
//
//
///r BUT don't factor out things like mysql - just have all relevant bits - i.e. little bit at end simply conditionally included?!?

///y make sure it handles parent=x,y,z;!?!

using std::cout;
using std::endl;
using std::runtime_error;

typedef unsigned int uint;
typedef unsigned char uchar; 
typedef boost::regex regex;
typedef boost::smatch smatch;
typedef feature_min feature;

// must put back check for non-unique names?!?
struct name_holder {
    std::set<std::string> scaffolds;
    std::set<std::string> gene_ids;
    std::set<std::string> parents;
    std::set<std::string> gene_transcript_ids;
};

// for header?!?
class print_row {

    long long bf;
    long long bm;
    bool blah;
    std::string issue;
    std::string item;

  public:

	print_row(std::string i, std::string a, long long _bf, long long _bm) : item(i), issue(a), bf(_bf), bm(_bm) {}
	print_row(std::string i, std::string a, bool b) : item(i), issue(a), blah(b) {}

	friend std::ostream& operator<< (std::ostream& os, const print_row& bobj) {

        bool problem = (bobj.bf&bobj.bm); // just to allow for testing of multiple bits
        // bool problem = (bitflag&bitmask)==bitmask;

        /// really lazy - clearly can't further overload << so just do it this horrible way?!?
        /// should really use polymorphism for this?!?

        std::string td, s; 
        if (bobj.bm==0) {
            td = bobj.blah ? "#FF927D" : "#e5f1ff";
            s =  bobj.blah ? "TRUE":"";
        } else { 
            td = problem ? "#FF927D" : "#e5f1ff"; // std::string td(problem?"#FF927D":"#e5f1ff");
            s =  problem ? "___________________TRUE____________" : "FALSE"; // std::string s(problem?"TRUE":"");
            // s =  problem ? "TRUE" : ""; // std::string s(problem?"TRUE":"");
        }

        td="<td BGCOLOR=\"" + td + "\" style=\"padding:10px;\" >";

		return os << " <tr> " << td 
          << bobj.item << " </td> " << td << bobj.issue 
          << " </td>" << td 
          << s << "</td></tr>\n"; //r break-up table with '\n' to avoid that odd html email formating issue?!?



	}

};

//////// really ought to abolish the use of set and make it a map so that can use same structure for both set difference and checking the gff
//////// do we do the extent during parsing or do we store everything in one structure that we can then use to add more tests too later?!?

/////// start storing strand?!? either way all we actually need for incomplete test is start/end of each exon
/////// could start checking for undefined exons - i.e. just make sure there's always an exon for every cds
/////// and check it's extension...?!?

std::ostream & nl(std::ostream& os) { return os << " \\\n"; } // nullary maninulator - should probably use an effector with arg...

#define PERMITTED_TRANSCRIPT_TYPES_REGEX    "mRNA|tRNA|pseudogenic_tRNA|rRNA|miRNA|ncRNA|pseudogene"
#define IGNORE_TYPES_REGEX                  "contig|supercontig|match|match_part"
#define CORE_PERMITTED                      "gene|exon|CDS|"
#define ALL_PERMITTED                       CORE_PERMITTED PERMITTED_TRANSCRIPT_TYPES_REGEX
#define MAX_BAD_LINES 5

/*  
#define REGEX_EMBL          "FT\\s.*"
#define REGEX_GFFFASTA      "#{1,2}FASTA.*"
#define REGEX_FASTA         ">\\w+.*"

#define REGEX_COM           "\\s*#.*"
#define REGEX_LINE_ENDS     ".*[^;]\\s*"
#define REGEX_BLANK         "\\s*"
#define REGEX_CR            "\\r"

#define REGEX_FEAT          "([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"
#define REGEX_NEGSTART      "([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t-\\d+\\t\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"
#define REGEX_NEGEND        "([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t-\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)"

#define REGEX_SCFNAME       "(\\w+?\\d+?):\\d+?[:-]\\d+?"
#define REGEX_ID            "ID=([^;]+)"
#define REGEX_PARENT        "Parent=([^;]+)"
*/

#define SPRINTF_STRING      "select count(1) from seq_region where name = '%s'"
#define SQL_STRING          "select * from seq_region where name = ''"

/// update to use c++0x?!?

/// clearly instead return int with 0 as validation 1 okay?!?! - thus just print of as desired!??
/// so these below are a sort of val1a set, then fragmentation/exon-overlap require features so more of val1b?!?

// should prolly use const int?!?
#define APOLLO_SCF_NAMES                            (1<<0) // (0x1<<0) 0x1
#define NAMES_HAVE_SPACES                           (1<<1)
#define LINES_WO_9COLS                              (1<<2)
#define LINE_ENDINGS                                (1<<3)
#define ILEGAL_FEAT_TYPES                           (1<<4)

#define EXCESS_GENE_CONSISTENCY_PROB                (1<<5)
#define EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB      (1<<6)
#define EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB  (1<<7)
#define EXCESS_CDS_EXON_CONSISTENCY_PROB            (1<<8)
#define NON_PERMITTED_BIOTYPES                      (1<<9)
#define ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE       (1<<10)
#define PARENT_WITHOUT_ID_NOT_CDS_EXON              (1<<11)

#define UNKNOWN_SCAFFOLD                            (1<<12)

#define NEGATIVE_COORDINATES                        (1<<13)
#define NON_PRINTING_X0D                            (1<<14)
#define BLANK_LINES                                 (1<<15)

#define NON_UNIQUE_ID                               (1<<16)

// no_CDS_without_pseudogene = false; // only processing protein coding via cap
#define PSEUDOGENES_PRESENT                         (1<<17)
#define CDS_PRESENT                                 (1<<18)

#define EMBL_FORMAT                                 (1<<19)
#define GFF_FASTA_HEADER                            (1<<20)
#define FASTA_HEADER                                (1<<21)

#define PARTIAL_MODEL                               (1<<22)

#define NO_FEATLINES                                (1<<23)
#define NO_GENES                                    (1<<24)
#define NO_EXON_CDS                                 (1<<25)
#define NO_TRANSCRIPTS                              (1<<26)

//// perhaps inline the individual checks?!?

static long long PROBLEM = EXCESS_GENE_CONSISTENCY_PROB | EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB 
    | EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB | EXCESS_CDS_EXON_CONSISTENCY_PROB | ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE 
    | PARENT_WITHOUT_ID_NOT_CDS_EXON | UNKNOWN_SCAFFOLD|NON_PERMITTED_BIOTYPES | NEGATIVE_COORDINATES 
    | NON_UNIQUE_ID | EMBL_FORMAT | GFF_FASTA_HEADER | FASTA_HEADER | PARTIAL_MODEL
    | NO_FEATLINES | NO_GENES | NO_EXON_CDS | NO_TRANSCRIPTS;

// PROBLEM = ~PROBLEM;
static long long NO_CDS_AND_NO_PSEUDOGENES = CDS_PRESENT | PSEUDOGENES_PRESENT;

namespace gff {

    // make these static, nest in separate namespace?!??!?
    /// no need for the regex's to be in headers?!?
    static regex reg_embl       ("FT\\s.*");
    static regex reg_gfffasta   ("#{1,2}FASTA.*");
    static regex reg_sequence   ("[actgnACTGN]+");
    static regex reg_fasta      (">\\w+.*");
    static regex reg_comment        ("\\s*#.*");
    static regex reg_line_ends  (".*[^;]\\s*");
    static regex reg_blank      ("\\s*");
    static regex reg_carriageret         ("\\r");

    static regex reg_feat       ("([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t\\d+\\t[^\\t]+\\t[+-]+\\t[^\\t]+\\t([^\\t]+)");
    static regex reg_negstart   ("([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t-\\d+\\t\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)");
    static regex reg_negend     ("([^\\t]+)\\t[^\\t]+\\t([^\\t]+)\\t\\d+\\t-\\d+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t([^\\t]+)");

    static regex reg_scfname    ("(\\w+?\\d+?):\\d+?[:-]\\d+?");
    static regex reg_id         ("ID=([^;]+)");
    static regex reg_parent     ("Parent=([^;]+)");

    static regex reg_ignore_types(IGNORE_TYPES_REGEX);
    static regex reg_allow_types(ALL_PERMITTED);
    static regex reg_transcript_biotypes(PERMITTED_TRANSCRIPT_TYPES_REGEX);
}

//// don't even need the gene names except to group transcripts if desired?!?

struct feature_ext : public feature_min {
    // explicit feature_ext(const feature_min& fm) : feature_min(fm) {}
    // char parent[MAX_ID_LENGTH]; // if anything stack smashes it's gonna be here!?!?
    feature_ext() = delete;
    feature_ext(const feature_min& fm) = delete;
    ~feature_ext () {}
    feature_ext(std::string p, uint a, uint b, unsigned char s) : _parent(p), feature_min(a,b), _strand(s) {} 
    std::string parent() const { return _parent; }   // don't ever return a non-const ref on a private member - i.e. you've just abolished privateness?!?
    unsigned char strand() const { return _strand; }
    // just forward args and default initialise?!?
    // feature_min(uint a, uint b) : _gstart(a), _gend(b), tstart(0), tend(0) {} 
    //
private :
    std::string _parent;
    unsigned char _strand;
};

struct gff_holder {
    //y overlap, fragmentation and names checks at..
    std::map<std::string,feature_ext> mrna_coords;
    // std::map<std::string,feature_min> mrna_coords;
    std::multimap<std::string,feature_min> exonbymrna_coords;
    //y check cds have exons and/or check all mRNAs have CDS and exons
    std::multimap<std::string,feature_min> cdsbymrna_coords;
};

/// let them do silly things like have different stranded-ness between exons, genes etc., - i.e. only use transcript strand - check consistency?!?
/// but add back non-unique id checks for gene/transcript?!?

long long gff_basic_validation_1a_gff_parse (const char* filename, std::stringstream& strstrm, name_holder& nh, gff_holder& gh) {
// long long gff_basic_validation_1a_gff_parse (const char* filename, std::stringstream& strstrm, name_holder& nh, gff_holder& gh, DB_PARAMS* dbp = 0) {

    cout << "\n[1a] parsing\n";

    // perhaps make this static?! - nah?!?
    long long bitflag = NO_FEATLINES|NO_GENES|NO_EXON_CDS|NO_TRANSCRIPTS;

    smatch match_obj;
    smatch match_obj2;

    std::ifstream in(filename);
    if(in == 0) throw runtime_error("problem opening gff file");

///////// needs planning?!?
//////// what do we want to test for really?!?
/////// names not adding up?!? - clearly a set is simplest
////// fragmentation - need start/end of mRNA via id and exons via parent
///// overlapping exons - same as above - i.e. simply the exons via parent
//// however, id for mrna_coords clearly makes transcript_ids superfluous - i.e. just extract the names from there?!?
/// also if start using extended structures for mRNA and/or storing cds then you've stored the same info multiple times e.g. transcript_parents in map of mrna...

/// some checks can be done at once e.g. cds->exon mapping - if no cds to map clearly there's an issue?!?

    // reduncant : std::set<std::string> transcript_parents;
    // redundant : std::set<std::string> transcript_ids;
    // get length of file: // is.seekg (0, ios::end); // length = is.tellg(); // is.seekg (0, ios::beg);

///y parse features
        
    std::string s;
    unsigned char error_counter = 0;
    long long ln = 0;
    while(getline(in, s)) { // Discards newline char

/// 1 : 'ordered'-ish series of checks of actual line format/headers?!?

//fstart

        // force ignoring of cr - but why if removing below and restarting?!? 
        if ((bitflag&NON_PRINTING_X0D)==NON_PRINTING_X0D) s.erase(s.end()-1, s.end()); //std::replace( s.begin(), s.end(), '\r', '');

        // horrible line endings - i.e. without final ';' can get a little messy
        if(regex_search(s,gff::reg_carriageret)) {
            if ((bitflag&NON_PRINTING_X0D)==NON_PRINTING_X0D) throw runtime_error ("this ain't good!"); // wtf?!?
            char cmdtmp[MED_STRING_SIZE];
            sprintf(cmdtmp,"perl -i -pe 's/\\r//' %s",filename); // wtf?!?
            if (system(cmdtmp) != 0) throw runtime_error ("problem removing carriage returns!");
            bitflag |= NON_PRINTING_X0D; // non_printing_x0d = true;
            in.seekg (0, std::ios::beg); //cout << "restarting : " << in.tellg() << endl;
        }

        ++ln;

        // inline bool _boring_checks(std::string& s, 
        // #define _boring_checks(std::string& s, 

        // blank lines and comments
        if(regex_match(s,gff::reg_blank)) {
          bitflag |= BLANK_LINES;
          continue;
        } 

        if(regex_match(s,gff::reg_line_ends)) bitflag |= LINE_ENDINGS; // line_endings without ';' - grrrr?!?;

        // silly headers 
        if(regex_match(s,gff::reg_gfffasta)) {
            bitflag |= GFF_FASTA_HEADER;
            // gff_fasta_header = true;
            strstrm << "<p>All lines must be valid gff feature lines or comments. Consequently, we do not accept submissions containing fasta sequence!</p>\n";
            continue;
        } else if(regex_match(s,gff::reg_fasta)) {
          bitflag |= FASTA_HEADER;
          continue;
        } else if(regex_match(s,gff::reg_embl)) {
          bitflag|=EMBL_FORMAT;
          continue;
        }
       
        // odd apollo thing of putting in negative coords - i.e. can't be bothered making feature matching too restrictive so catch it here?!?
        if (regex_match(s,gff::reg_negstart) || regex_match(s,gff::reg_negend)) {
            bitflag |= NEGATIVE_COORDINATES; 
            continue;
        }

        // do last or mis-diagnose the annoyance above?!?
        if(regex_match(s,gff::reg_comment)) continue;

        // if (regex_match(s,match_obj,gff::reg_feat)) { // don't do this just do a continue if not?!? - no need to nest?!?
        if (!regex_match(s,match_obj,gff::reg_feat)) { 
            bitflag |= LINES_WO_9COLS;

            // be more helpful?!? already know the line isn't properly formed - i.e. with 8 tabs etc.?!?
            if(error_counter<MAX_BAD_LINES) { 
              strstrm << "<p>line " << ln << " does not have gff 9-col, tab-delimited format :</p>\n<p>";

                int tabs = count(s.begin(),s.end(),'\t');
                int spaces = count(s.begin(),s.end(),' ');

                if(tabs==8 && spaces < 7) // if(tabs==8) {
                strstrm << "- start & end must be numeric and strand must conform to '+' or '-' " << ln << " : ";
                else if (tabs > 5 && tabs < 11) // put in check on number of spaces?!?
                strstrm << "- please check all 9 columns are present " << ln << " : ";
                else if (spaces > 7 && spaces < 9 && tabs < 3) // tabs check?!?
                strstrm << "- is this line delimited with spaces and not tabs? " << ln << " : ";
                else if (regex_match(s,gff::reg_sequence)) 
                strstrm << "- this line looks like sequence. " << ln << " : ";
                else 
                strstrm << "- i'm really not sure what this is? " << ln << " : ";

                strstrm << s << "</p>\n";

            } else if (error_counter==MAX_BAD_LINES) 
              strstrm<<"<p>Non 9-column, tab-delimited format errors truncated.</p>\n";

            error_counter++;
            continue;
        }

//fend

/// 2 : we have a feature

//fstart
        //y really not worth hassle of parsing col9 without regex?!?

        bitflag &= ~NO_FEATLINES; // featlines++;

        std::stringstream featurestrm(s); // lazy but we know it adheres to correct format so can't be bothered matching?!?
        std::string scfname, score, strand, source, type, ignore, annot;
        int start = 0, end = 0; // wtf?!? - wasn't initialised?!?!
        featurestrm >> scfname >> source >> type >> start >> end >> score >> strand >> ignore >> annot;

        // std::string annot; if (match_obj[3].matched) annot = match_obj[3]; ///////// what the bollocks is this?!? wtf - why match too?!?

        if (regex_match(type,gff::reg_ignore_types)) continue; // ignore contigs etc.

        if (!regex_match(type,gff::reg_allow_types)) {
            strstrm << "<p>Ilegal feature type='" << type << "'</p>\n";
            bitflag |= NON_PERMITTED_BIOTYPES;
        }
        
        // horrible apollo ids with coords inserted?!?
        if (!scfname.empty() && regex_match(scfname,match_obj,gff::reg_scfname)) {
            nh.scaffolds.insert(match_obj[1]);
            bitflag|=APOLLO_SCF_NAMES; 
        } else if (!scfname.empty()) nh.scaffolds.insert(scfname); // completely unecessary conditional

//fend 

/// 3 : we have a 2o or 3o feature (have id and parent) : transcript or exon/CDS features

//fstart 

        if (regex_search(annot,match_obj,gff::reg_id)&&regex_search(annot,match_obj2,gff::reg_parent)){ 

            std::string id = match_obj[1];
            std::string parent = match_obj2[1];

            // we actually ifnogre this...
            if(id.find(' ')!=std::string::npos || parent.find(' ')!=std::string::npos) bitflag|=NAMES_HAVE_SPACES;

            if(type == "CDS") bitflag|=CDS_PRESENT;

            if(type == "exon") gh.exonbymrna_coords.insert(std::pair<std::string,feature_min>(parent,feature_min(start,end))); 

            if (type == "exon" || type == "CDS") {
                bitflag &= ~NO_EXON_CDS;
                nh.parents.insert(parent);
            } else if (regex_match(type,gff::reg_transcript_biotypes)) {
                bitflag &= ~NO_TRANSCRIPTS;
                // handled by mrna_coords below : transcript_ids.insert(id);
                // handled by mrna_coords below : transcript_parents.insert(parent);
                gh.mrna_coords.insert(std::pair<std::string,feature_ext>(id,feature_ext(parent,start,end,strand=="+"?1:0))); 
                nh.gene_transcript_ids.insert(id);
                
                /*
                // if(transcript_ids.count(id)) non_unique_id = true;
                if(transcript_ids.count(id)) {
                    non_unique_id = true;
                    strstrm << "> id " << id << " is not unique"<< endl;
                }
                /
                // book-keeping for gene/mRNA name overlap
                // if(gene_transcript_ids.count(id)) non_unique_transcript_id = true;
                if(gene_transcript_ids.count(id)) { 
                    non_unique_id = true;
                    strstrm << "> id " << id << " is not unique"<< endl;
                }
                */

            } else {
                bitflag |= NON_PERMITTED_BIOTYPES;
                strstrm << "<p>Not accepting biotype " << type << " at transcript level</p>\n";
            }

//fend

/// 4 : we have 1o features - only have id

//fstart 

        } else if (regex_search(annot,match_obj,gff::reg_id)) {   

            std::string id = match_obj[1];

            if(id.find(' ')!=std::string::npos) bitflag|=NAMES_HAVE_SPACES;

            if(type == "pseudogene") bitflag|=PSEUDOGENES_PRESENT;

            if (type == "gene" || type == "pseudogene") {

                bitflag &= ~NO_GENES;
                nh.gene_ids.insert(id);
                nh.gene_transcript_ids.insert(id);

                /*
                // if(gene_ids.count(id)) non_unique_gene_id = true;
                if(gene_ids.count(id)) {
                    non_unique_id = true;
                    strstrm << "> id " << id << " is not unique"<< endl;
                }

                // book-keeping for gene/mRNA name overlap
                // if(gene_transcript_ids.count(id)) non_unique_transcript_id = true;
                if(gene_transcript_ids.count(id)) {
                    non_unique_id = true;
                    strstrm << "> id " << id << " is not unique"<< endl;
                }
                */

            } else {
                cout << "REALLY NASTY : line=" << __LINE__ << " type=" << type << " file=" << __FILE__ << endl;
                bitflag |= ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE;
                strstrm << "<p>The following is not a gene or pseudogene and thus must have a Parent tag: "<<s<<"</p>\n";
            }

//fend

/// 5 : we have 3o feautres - only parent

//fstart 

        } else if (regex_search(annot,match_obj2,gff::reg_parent)) {

            std::string parent = match_obj2[1];

            if(parent.find(' ')!=std::string::npos) bitflag|=NAMES_HAVE_SPACES;

            if(type == "CDS") bitflag|=CDS_PRESENT;

            if(type == "exon") gh.exonbymrna_coords.insert(std::pair<std::string,feature_min>(parent,feature_min(start,end))); 

            if (type == "exon" || type == "CDS") {
                bitflag &= ~NO_EXON_CDS;
                nh.parents.insert(parent);
            } else {
                cout << "REALLY NASTY : line=" << __LINE__ << " type=" << type << " file=" << __FILE__ << endl;
                bitflag |= PARENT_WITHOUT_ID_NOT_CDS_EXON;
                strstrm << "<p>The following is not a CDS or exon and thus must have an ID tag: "<< s << "</p>\n";
            }

        } else {}

//fend

    } // end of file 

    return bitflag;

}

long long gff_basic_validation_1b_gff_name_checks (long long bitflag, std::stringstream& strstrm, name_holder& nh, gff_holder& gh) {

///y simple checks for feature id/parent names...

cout << "\n[1b] doing id/parent name checks!??!\n";

    // should truncate messages here?!? - i.e. just have an externally scoped counter passed by reference?!?
    // error_counter = 0;

//fstart

    // using namespace boost::lambda;
    // [](){ cout << "poo";}();
    // put this in an inline function?!?

    std::set<std::string> transcript_parents;
    for_each(gh.mrna_coords.begin(),gh.mrna_coords.end(),[&transcript_parents](std::pair<std::string,feature_ext> e){ transcript_parents.insert(e.second.parent()); });

    std::vector<std::string> top_level_problems;
    set_difference(nh.gene_ids.begin(),nh.gene_ids.end(),transcript_parents.begin(),transcript_parents.end(), std::back_inserter(top_level_problems));
    if (top_level_problems.size() > 0) {
        strstrm <<"<p>Top level set difference from genes_ids with no transcript_parents (missing transcript/excess gene features) = " <<top_level_problems.size()<<"</p>\n";
        bitflag |= EXCESS_GENE_CONSISTENCY_PROB;
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * Excess 'gene' ID=" << s << "</p>\n"; });
      // strstrm << constant("<p> * Excess 'gene' ID=") << _1 << constant("</p>\n")

    top_level_problems.clear();

    set_difference(transcript_parents.begin(),transcript_parents.end(),nh.gene_ids.begin(),nh.gene_ids.end(),std::back_inserter(top_level_problems));
    if (top_level_problems.size() > 0) {
        strstrm << "<p>Top level set difference from transcript_parents with no gene_ids (excess transcript/missing gene features) = " <<top_level_problems.size() << "</p>\n";
        // "with no gene_ids (excess transcript/missing gene features) = " <<top_level_problems.size() << "</p>\n";
        bitflag |= EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB;
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * No 'gene' corresponding to transcript Parent=" << s << "</p>\n"; });
      // strstrm << constant("<p> * No 'gene' corresponding to transcript Parent=") << _1 << constant("</p>\n")

    top_level_problems.clear();

    std::set<std::string> transcript_ids;
    for_each(gh.mrna_coords.begin(),gh.mrna_coords.end(),[&transcript_ids](std::pair<std::string,feature_ext> e){ transcript_ids.insert(e.first); });

    set_difference(transcript_ids.begin(),transcript_ids.end(),nh.parents.begin(),nh.parents.end(),std::back_inserter(top_level_problems));
        // std::inserter(top_level_problems, top_level_problems.end()) - was using set
    if (top_level_problems.size() > 0) {
        strstrm << "<p>Top level set difference from transcript_ids with no cds/exon_parents  (missing cds-exon/excess transcript features) = "
        <<top_level_problems.size() << "</p>\n";
        bitflag |= EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB;
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * Excess 'transcript' ID=" << s << "</p>\n"; } );
      // strstrm << constant("<p> * Excess 'transcript' ID=") << _1 << constant("</p>\n")

    top_level_problems.clear();

    set_difference(nh.parents.begin(),nh.parents.end(),transcript_ids.begin(),transcript_ids.end(),std::back_inserter(top_level_problems));
        // std::inserter(top_level_problems, top_level_problems.end()) - was using set
    if (top_level_problems.size() > 0) {
        strstrm << "<p>Top level set difference from cds/exon parents with no transcript_id (excess cds-exon/missing transcript features) = "<<top_level_problems.size() << "</p>\n";
        bitflag |= EXCESS_CDS_EXON_CONSISTENCY_PROB; 
    }
    for_each(top_level_problems.begin(),top_level_problems.end(),[&strstrm](std::string s) { strstrm << "<p> * No 'transcript' corresponding to cds/exon Parent=" << s << "</p>\n"; });


//fend

    return bitflag;

}

long long gff_basic_validation_1d_individual_model_checks (long long bitflag, std::stringstream& strstrm, gff_holder& gh) {

    cout << "\n[1b] individual model checks!??!\n";

///y fragmneted model check

//fstart

    //// e! api does NOT use mrna/gene coords and so allows genes to be truncated without problem!?! THUS must check this?!?
    
    ///// this can replace most of the set_diff stuff?!? - i.e. if we don't find the elements!?!
    for (auto it = gh.mrna_coords.begin(); it!=gh.mrna_coords.end(); ++it) {

        int count = gh.exonbymrna_coords.count(it->first);

        switch (count) {

            case 0: std::cout << "move the book-keeping stuff to here?!?"<<std::endl; 
            break;

            case 1: { 

                // following two are equivalent - but not according to standard!?!?
                // std::_Rb_tree_iterator<std::pair<const std::string, feature_min> > exon_it = exonbymrna_coords.find(it->first);
                std::multimap<std::string,feature_min>::iterator exon_it = gh.exonbymrna_coords.find(it->first);
                
                if ((it->second).gstart()!=(exon_it->second).gstart()||(it->second).gend()!=(exon_it->second).gend()) {
                    std::cout << "[single exon] FRAGMENTED MODEL "<<it->first << std::endl; 
                    strstrm << "<p>Transcript "<< it->first << " is fragmented (exon features do not fully account for transcript extension; "
                      <<"transcript: "<<(it->second).gstart()<<"-"<<(it->second).gend()<<", exons: "<<(exon_it->second).gstart()<<"-"<<(exon_it->second).gend()<<").</p>\n";
                    bitflag |= PARTIAL_MODEL;
                }

            break; }

            default: /* declaring vars in switch!?! */ { //std::cout << "check the extent" <<std::endl;

                auto mmit = gh.exonbymrna_coords.equal_range(it->first);

                int min_start = mmit.first->second.gstart();
                int max_end = mmit.first->second.gend();

                int i=0;
                for (auto exon_it = ++mmit.first ; exon_it != mmit.second ; ++exon_it) { 
                    if((exon_it->second).gstart()<min_start) min_start=(exon_it->second).gstart();
                    if((exon_it->second).gend()>max_end) max_end=(exon_it->second).gend();
                }

                if ((it->second).gstart()!=min_start||(it->second).gend()!=max_end) {
                    std::cout << "[multi exon] FRAGMENTED MODEL " << it->first << std::endl; 
                    strstrm << "<p>Transcript "<< it->first << " is fragmented (exon features do not fully account for transcript extension; "
                      <<"transcript: "<<(it->second).gstart()<<"-"<<(it->second).gend()<<", exons: "<<min_start<<"-"<<max_end<<").</p>\n";
                    bitflag |= PARTIAL_MODEL;
                }
            break; }
        };

    }

//fend


///y overlapping exon check?!?

//fstart
//fend

    return bitflag;

}

long long gff_basic_validation_1c_scfnames (long long bitflag, std::stringstream& strstrm, name_holder& nh, DB_PARAMS* dbp) {

///y scaffold names - needs to be converted to use mysql adaptor?!? - do once integrating into capmon?!? - changes must be only local additions

#ifdef CAPMON_EXT

cout << "\n[1c] scaffold name checks\n";

//fstart 

///  the mysql part - i.e. e! dependent checks are purely about scaffolds being known - i.e. non-coding cds is via loader?!?
/// either way these are fairly simple to go from either flat-file of db but need to be separate to stop dependencies
/// just required std::set<std::string>::iterator of names to check for?!?

    try {

        MYSQL *conn;
        MYSQL_RES *result;
        MYSQL_ROW row;
        MYSQL_FIELD *field;
        conn = mysql_init(NULL);

        if(mysql_real_connect(conn,dbp->host(),dbp->user(),dbp->pass(),dbp->dbname(),dbp->port(),NULL,0) == NULL)
          throw MySqlError("couldn't connect to database. exiting.");

        for (std::set<std::string>::iterator sit = nh.scaffolds.begin() ; sit != nh.scaffolds.end() ; sit++) { // prolly ought to be a functor...

            // options: (1) bind parameters/statement prepare - but the actual binding part is pretty darned ugly!?! - http://dev.mysql.com/doc/refman/5.0/en/mysql-stmt-execute.html
            // http://dev.mysql.com/doc/refman/5.0/en/mysql-stmt-execute.html
            char qbuf[STRING_SIZE]; // (2)  sprintf...
            sprintf(qbuf, SPRINTF_STRING, sit->c_str());

            if (mysql_query(conn,qbuf))
              throw MySqlError("unable to access seq_region table - is this really a e!/cap db?!?");

            if(!(result = mysql_store_result(conn))) { //y get the result set
                throw runtime_error ("unable to query database for seq_region!");
            } else {
                row = mysql_fetch_row(result);
                /// if just selecting then this would indicate no entry - but that's not v. safe
                /// better to do a count which guarantees a result and allows other f'ups to be distinguished
                if (row == 0) {
                    throw runtime_error ("null pointer returned on db query!");
                } else {
                    if (atoi(row[0]) == 0) {
                        strstrm << "<p>There is no scaffold named " << sit->c_str() << " in cap db</p>\n";
                        bitflag |= UNKNOWN_SCAFFOLD;
                    }
                }
            }
            mysql_free_result(result);
        }
        ///y freeing here will give invalid ponters!?! - duh
        mysql_close(conn);

    } catch (std::runtime_error e) {
        std::string prob("Exception thrown : ");
        throw; // propagate up...
    }

//fend

#else

cout << "\n[1c] MUST compile with CAPMON_EXT for scaffold name checks aagainst mysql/e! db instance\n";

#endif

    return bitflag;

}

void gff_basic_validation_1x_file_cleanup (long long bitflag,const char* filename,std::stringstream& strstrm) {

///y local copy cleanup?!?

    cout << "\n[1x] file cleanup\n";

//fstart

///// why?!? - at least put post-processing clean-up elsewhere and make it conditional on capmon?!?
    //y must clean this first as otherwise putting in line-endings means parent name will include \r etc..
    if ((bitflag&LINE_ENDINGS)==LINE_ENDINGS) {
        char cmdtmp[MED_STRING_SIZE];
        sprintf(cmdtmp,"perl -i -pe 's/(.*[^;]\\s*)\\n/${1};\\n/' %s",filename);
        // strstrm << "> fixing line endings in-place" << endl;
        // strstrm << cmdtmp << endl;
        //if (system(cmdbuf) == -1) {
        if (system(cmdtmp) != 0) 
          throw runtime_error ("problem fixing line endings!");
        // perror("Error executing input");
    }

    if ((bitflag&APOLLO_SCF_NAMES)==APOLLO_SCF_NAMES) { // if (bitflag&APOLLO_SCF_NAMES) {  // if (apollo_scf_names) {

        // regex gff::reg_scfname("(\\w+?\\d+?):\\d+?[:-]\\d+?");
        char cmdtmp[MED_STRING_SIZE];
        // sprintf(cmdtmp,"perl -i -pe 's/^(\\w+\\d+):\\d+-\\d+\\t/${1}\\t/' %s
        sprintf(cmdtmp,"perl -i -pe 's/^(\\w+\\d+):1-\\d+\\t/${1}\\t/' %s",filename);
        strstrm << "<p>Fixing apollo scaffold name endings in-place</p>\n";
        if (system(cmdtmp) != 0) 
          throw runtime_error("problem fixing apollo scaffold name endings!");
    }

    //r make sure we fix this last - i.e. after ';' etc..
    if (bitflag&NAMES_HAVE_SPACES==NAMES_HAVE_SPACES) { // if (names_have_spaces) {

        char cmdtmp[MED_STRING_SIZE];
        sprintf(cmdtmp,"perl -i -pe 'if(/^(.*ID=)([ \\S]+?)(;.*)$/) { $c=$1;$d=$2;$e=$3;$d=~s/ /_/g; $_= qq{$c$d$e\\n} }' %s",filename);
        strstrm << "<p>Removing spaces from IDs in-place</p>\n";
        if (system(cmdtmp) != 0) 
          throw runtime_error("problem removing spaces from IDs!");
        sprintf(cmdtmp,"perl -i -pe 'if(/^(.*Parent=)([ \\S]+?)(;.*)$/) { $c=$1;$d=$2;$e=$3;$d=~s/ /_/g; $_= qq{$c$d$e\\n} }' %s",filename);
        strstrm << "<p>Removing spaces from Parents in-place</p>\n";
        if (system(cmdtmp) != 0) 
          throw runtime_error("problem removing spaces from Parents!");
    }
///////////////////////////////////////////////////////////////////////////////////////////////////

//fend


}

long long generic_validation_and_store (const char* filename, std::string& report, DB_PARAMS* dbp = 0) {
    gff_holder gd; // if scaffold names with fasta perhaps put that in else directive with mysql stuff?!?
    // could call gff_basic_validation_1a_gff_parse then perhaps even do fasta checks?!?
    // should the mysql stuff go into it's own function - i.e. that has empty body without inclusion?!?
}

// have summary separate so have single routine that returns bitflag and data structure that capmon calls - use it for test all flags against appropriate files?!?

bool capmon_gff_validation (const char* filename, std::string& report, DB_PARAMS* dbp = 0) {

    gff_holder gh; // just let it go out of scope?!?

    std::stringstream strstrm(std::stringstream::out);

    long long bitflag = 0;

{ // just avoiding persistance of name_holder for no reason?!?

    name_holder nh;

    bitflag = gff_basic_validation_1a_gff_parse (filename, strstrm, nh, gh);
    // bitflag = gff_basic_validation_1a_gff_parse (filename, strstrm, nh, gh, dbp);
    // long long bitflag = gff_basic_validation_1a_gff_parse (filename, strstrm, dummy, dbp);

    bitflag = gff_basic_validation_1b_gff_name_checks (bitflag, strstrm, nh, gh);

//     if(dbp) // running directly without CAPMON_EXT will just moan at you - i.e. don't want unecessary linking dependencies?!?
      bitflag = gff_basic_validation_1c_scfnames (bitflag, strstrm, nh, dbp);

    //y could put in fasta checking?!?

}

    bitflag = gff_basic_validation_1d_individual_model_checks(bitflag, strstrm, gh);

    gff_basic_validation_1x_file_cleanup (bitflag,filename,strstrm);






    /////// summary and capmon decision?!?
       
    std::stringstream strstrm2(std::stringstream::out);

    // print_row using simple | thus if any of the flags are true will report
    strstrm2 << "<p></p><table style=\"width:500px;\">\n"

      << print_row("A01","Contains EMBL qualifiers",bitflag,EMBL_FORMAT)

      << print_row("A02","Contains fasta entries",bitflag,GFF_FASTA_HEADER|FASTA_HEADER)
      // << print_row("A02","Contains fasta entries",gff_fasta_header||fasta_header)

      << print_row("A03","Does not contain gene features",bitflag,NO_GENES) 
      << print_row("A04","Does not contain transcript features",bitflag,NO_TRANSCRIPTS)
      << print_row("A05","Does not contain exon/CDS features",bitflag,NO_EXON_CDS)
      // << print_row("A03","Does not contain gene features",genes==0) 
      // << print_row("A04","Does not contain transcript features",transcripts==0)
      // << print_row("A05","Does not contain exon/CDS features",exon_cds==0)

// could check number of tabs and/or spaces in line?!? i.e. if there are 7-15 spaces state they should check delimiter
// if no spaces or tabs - are they all actgnACTGN - looks like you have sequence data?!?
// if 6-15 tabs but not proper format - check all columns are present
// else tell them to check the individual column values?!?

    // put in message about what each column must be?!?
    // about making sure there are all columns
    // and about tabs not spaces?!?
      << print_row("A06","There are non-comment, non-9-col format lines present",bitflag,LINES_WO_9COLS)
      // << print_row("A06","There are non-comment, non-9-col format lines present (make sure 'all' columns are present and file is delimited with tab not space)",bitflag,LINES_WO_9COLS)
      << print_row("A07","Gene id excess vs. transcript parent",bitflag,EXCESS_GENE_CONSISTENCY_PROB)
      << print_row("A08","Transcript parent excess vs. gene id",bitflag,EXCESS_TRANS_REL_GENE_CONSISTENCY_PROB)
      << print_row("A09","Transcript id excess vs. exon/CDS parent",bitflag,EXCESS_TRANS_REL_CDS_EXON_CONSISTENCY_PROB)
      << print_row("A10","Exon/CDS parent excess vs. transcript id",bitflag,EXCESS_CDS_EXON_CONSISTENCY_PROB)



      << print_row("A12","Non-permitted biotypes",bitflag,NON_PERMITTED_BIOTYPES)
      << print_row("A13","Negative coordinates - known Apollo bug",bitflag,NEGATIVE_COORDINATES)

      // really crude but also needs a specific check on the combination
      // << print_row("A14","Missing CDS - may be bug in previous Artemis versions",bitflag,CDS_PRESENT|PSEUDOGENES_PRESENT)
      // << print_row("A14","Missing CDS - may be bug in previous Artemis versions",!CDS_present&&!pseudogenes_present)

      << print_row("A15","Non-gene features with ID but no Parent ",bitflag,ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE)
      << print_row("A16","Non-CDS/exon features with Parent but no ID",bitflag,PARENT_WITHOUT_ID_NOT_CDS_EXON)
      << print_row("A17","Broken (partial) transcript(s)",bitflag,PARTIAL_MODEL)

// requires specific combination but either way is horribly crude - should test all the genes individually?!?
// probably better done in gff_checks_1b?!? - i.e. don't have that sort of thing on by default?!?
// checking of names - this clearly might as well get done first?!?
// exons w/o cds?!?
// overlapping exons?!?
// fragmentation?!?

      << "</table><p></p>";

//// seem to be merging reports?!?
     // strstrm2 << strstrm.rdbuf();
    if (!strstrm.str().empty()) {
        std::string blah = "\n\n<h3>DETAILED REPORT:</h3>\n" + strstrm.str();
        strstrm2 << blah;
        // strstrm2 << "</pre>";
    }
    // cout << strstrm.str() << endl;

    report = strstrm2.str();

/// this actually goes later - i.e. in caller?!?
/// to do this can either turn off ALL bits apart from the ones of interest and test?!? e.g. bm & ~(mask1|mask2...)!=0
/// alternatively - and probably more readibly just a compound mask directly and test e.g. if (bm & multiplemask)==multiplemask)
/// i.e. exactly same as testing for a single bit?!? - BUT THAT IS FOR TESTING THEM ALL AT ONCE!?!?

//       |(!cds_PRESENT&&!PSEUDOGENES_PRESENT)
    
    ///r will return the bitflag at this point - then interogate it in caller?!?

    bool no_cds_and_no_pseudogenes = (bitflag & NO_CDS_AND_NO_PSEUDOGENES)==NO_CDS_AND_NO_PSEUDOGENES; // i.e. testing the exact combination?!?
    // i.e. if 'any' of the bit masks are active stop
    bool bitwise_stop = bitflag & PROBLEM; // PROBLEM = ~PROBLEM; //  bool bitwise_proceed = bitflag & PROBLEM;

    // if (genes==0||transcripts==0||exon_cds==0||no_cds_and_no_pseudogenes||bitwise_stop) {
    if (no_cds_and_no_pseudogenes||bitwise_stop) return false;
    else return true;
}

int main () {

/// so we get this going in isolation and then it can be used more generally - and allow new tests to be added much more easily?!?!

// run the file formating only checks
// vector<features> x?!?
    // long long bitflag = gff_basic_validation_1a_gff_parse(x);
    //
    // gff_checks_1b - exon coverage?!?
    // overlapping exons?!?...
    //
    // i.e. have a validation level - low stringency only does 1a
    // higher stringency does the actual model checks?!?
    // 1a only bothers populating a structure of features if validation level is higher?!?
    //
    // these names are likely to be keys in map datastructure : bool unknown_scaffold = false;
//       << print_row("A11","Unknown scaffold names" ,unknown_scaffold) 





    cout << "gonna check file = ";
    std::string summary;

    cout << " " << std::boolalpha << capmon_gff_validation("exonprob.gff",summary)<<"\n\n";

    cout << " " << std::boolalpha << capmon_gff_validation("Yb.gff3",summary)<<"\n\n";

    // cout << " " << std::boolalpha << gff_basic_validation_1a_gff_parse("../gff/Yb_superscaffold_v1-0.gff3",summary)<<"\n\n";

    cout << summary;

    cout << "\n\ndone\n\n";
}

    /*
    while(getline(in, s)) { // Discards newline char
        // problem is that the file is buffered and we now need to wipe the damned thing every single line internally - better close n re-open?!?
        // cout << "at position: " << in.tellg() << endl;
        // really should use full match on .*\r\n - i.e. \r\n$...
        if(regex_search(s,gff::reg_cr))
          non_printing_x0d = true;
    }
    if (non_printing_x0d) {
        char cmdtmp[STRING_SIZE];
        //y should specifically use \\r\\n$
        sprintf(cmdtmp,"perl -i -pe 's/\\r//' %s",filename);
        strstrm << "> removing carriage returns in-place" << endl;
        if (system(cmdtmp) != 0) 
            throw runtime_error ("problem removing carriage returns!");
    }
    in.close();
    in.open();
    */
    /*
    ///y for set_difference simply need sorted ranges (associative containers are always sorted!?!) and an output 
    ///y iterator to put the difference in - i.e. one that supports * = ++ - deref, assigment and incrementation pointer operations etc..
    ///b hence can shove back into a set, but since it must be unique as we're using sets there's little point
    std::vector<std::string> hey;
    // std::vector<std::string>::iterator it = set_difference(
    set_difference(
        gene_ids.begin(),gene_ids.end(),
        transcript_parents.begin(),transcript_parents.end(), 
        ///r back_inserter() - can be used if the recipient container supports the push_back() member function
        //back_inserter(container_supporting_push_back) 
        //front_inserter(container_supporting_push_front) 
        std::back_inserter(hey));
    */