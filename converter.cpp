#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include "feature.h"
using std::cout;
using std::endl;
using std::runtime_error;

/// genomic coords always absolute/transcript strand dependent?!?

class Object {
    struct ObjectConcept {  
        virtual ~ObjectConcept() {} 
        virtual bool has_concept_of_some_property() const = 0;
        virtual std::string name() const = 0;
    };

    template<typename T> struct ObjectModel : ObjectConcept {
        ObjectModel( const T& t ) : object( t ) {}
        virtual ~ObjectModel() {} // not needed it's already virtual in base?!?
        virtual bool has_concept_of_some_property() const { return object.has_some_property_blah(); } 
        virtual std::string name() const { return typeid(object).name(); }
      private:
        T object; // the actual stored type
    };
    std::shared_ptr<ObjectConcept> object; 
  public:
    template<typename T> Object(const T& obj) : object( new ObjectModel<T>( obj ) ) {}
    std::string name() const { return object->name(); }
    bool has_concept_of_some_property() const { return object->has_concept_of_some_property(); }
};

// inheritance here is not for polymorphism it just forces common minimal features with same types?!?
struct feature_converter : public feature_min {

//// have a conversion ctor from feature_min!?!?!

    explicit feature_converter(const feature_min& fm) : feature_min(fm) {}

    feature_converter() = delete;
    ~feature_converter () {}
    feature_converter(uint a, uint b) : feature_min(a,b), tstart(0), tend(0) {} // just forward args and default initialise?!?
    // feature_min(uint a, uint b) : _gstart(a), _gend(b), tstart(0), tend(0) {} 
    uint tstart;
    uint tend;
};

// should prolly keep transcript start/stop - i.e. to make sure the transcript has exons spanning full extension?!?
class transcript {

    transcript()=delete;
  private:

    uchar _strand;
    uint _gpstart;
    uint _gpend;
    uint _peptide_transcript_offset;
    uint _peptide_transcript_length;
    std::vector<feature_converter> exons;
    /// seems horrible to store cds too but alternative it to iterate through all of them for all transcritps for no reason?!?
    std::vector<feature_converter> cds;

  public:

    transcript(uchar s) : _strand(s), _gpstart(0), _gpend(0), _peptide_transcript_offset(0), _peptide_transcript_length(0) { if (_strand!=1&&_strand!=0) 
      throw runtime_error("must be 1 or 0"); };

    void addexon(const feature_converter& f) { exons.push_back(std::move(f)); }
    void addcds(const feature_converter& f) { cds.push_back(std::move(f)); }

    void printexons() {
        for_each(exons.begin(),exons.end(), [](feature_converter& f) { cout << "transcript: "<< f.gstart() <<" "<<f.gend() << " " << f.tstart << " " << f.tend<< "\n"; });
    }

    uint convertt2g(uint i) {

        if (!exons.at(0).tstart) {
            _sortexons();
            _fill_in_coords();

            if(cds.size()>0) _fill_in_peptide();
        } 

        // cout <<"will convert " << i<<"\n";
        for(auto it=exons.begin(); it!=exons.end(); it++) {

            if(it->tstart<=i&&it->tend>=i) {
                // cout << "we have our exon " << it->tstart << "-"<<it->tend<<"\n";
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
                if (_strand) return it->gstart()+i-it->tstart;
                else return it->gend()-i+it->tstart;
            }


        }
        
        throw runtime_error("argument must fall within range of transcript");
        
    }

    uint convertg2t(uint i) {

        if (!exons.at(0).tstart) {
            _sortexons();
            _fill_in_coords();

            if(cds.size()>0) _fill_in_peptide();
            // do we fill in cds details here?!? - treat them in the same way or just give exons extra parameter?!?
            // and/or phase to allow them to deal with peptide location?!?
            //
            // just find peptide absolute genomic positions?!?
            // do we do any atual qc on the models?!?
        } 


        // cout <<"will convert " << i<<"\n";
        /// how to bother with introns?!? - could be horribly lazy and use two loops?!?
        for(auto it=exons.begin(); it!=exons.end(); it++) {

            if(it->gstart()<=i&&it->gend()>=i) {
                cout << "we have our exon " << it->tstart << "-"<<it->tend<<"\n";
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
                if (_strand) return i-it->gstart()+it->tstart;
                else return it->gend()-i+it->tstart;
            }


        }

        auto it2=exons.begin();
        std::vector<feature_converter>::iterator it1;
        int intron = 0;
        while (1) {

            // auto it1=it2++;
            it1=it2++;
            intron++;
            if (it2==exons.end()) break;
             cout << "using "<<it2->gend()+1 <<" - "<<it1->gstart()-1 << "\n";
        
            if((_strand&&it1->gend()+1<=i&&it2->gstart()-1>=i)||(!_strand&&it2->gend()-1<=i&&it1->gstart()+1>=i)) {
                cout << "it's in intron " << intron << "\n";// : "; //  << 
                // it1->gend()+1 << "-" << it2->gstart()-1<<"\n";
                // _strand?(it1->gend()+1):(it2->gend()-1) << "-" << _strand?(it2->gstart()-1):(it1->gstart()+1) << "\n";
                return 0;
                // return 0;
                // cout <<"and "<<it->gstart()<<" " <<it->tstart<<"\n";
            }
        }
        
        throw runtime_error("argument must fall within range of transcript : " + std::to_string(i));

    }

  private:
    // template<typename T> inline void doit(T& it, T& eit) {

    void _sortexons() {
        // could use std::function?!?
        // cout << "we have strand="<<_strand<<"\n";
        if (_strand) sort(exons.begin(),exons.end(),[](feature_converter a, feature_converter b) { return (a.gstart()<b.gstart()); });
        else sort(exons.begin(),exons.end(),[](feature_converter a, feature_converter b) { return (a.gstart()>b.gstart()); });
    }

    void _sortcds() {
        // if (_strand) 
            sort(cds.begin(),cds.end(),[](feature_converter a, feature_converter b) { return (a.gstart()<b.gstart()); });
        // else sort(cds.begin(),cds.end(),[](feature_converter a, feature_converter b) { return (a.gstart()>b.gstart()); });
    }

    void _fill_in_coords() {
        int i=0,last_te=0;
        for (auto it=exons.begin() ; it!=exons.end(); it++,i++) { // for ( ; it!=eit; it++,i++) {
            // cout <<"i="<<i<< " " << it->gstart() << "\n";
            if(i==0) {
                it->tstart=1; // already 0 initialised tstart?!?!
                // use ternaries e.g. it->tend = strand ? it->gend-it->gstart+1 : ...?!?
                it->tend=it->gend()-it->gstart()+1;
            } else {
                it->tstart=last_te+1;
                it->tend=it->gend()-it->gstart()+it->tstart; 

            }
            last_te=it->tend;
        }
    }

    void _fill_in_peptide() {

        // if we're not gonna do model qc no actual reason to order the cds?!?
        _sortcds();
/// this checking should be an option?!?
        /* do we really care about this?!? - just use first and last?!?
        for (auto cit=cds.begin() ; cit!=cds.end(); cit++) { 
            
            cout << "check cds " << cit->gstart() << "-" << cit->gend() << "\n";
            bool found=false;
        
            for (auto eit=exons.begin() ; eit!=exons.end(); eit++) { 
                if (cit->gstart()>=eit->gstart()&&cit->gend()<=eit->gend()) {
                    cout << " heya!!!\n";
                    found=true;
                }
            }
            if (!found) throw std::string("unable to find corresponding exon - this gff model doesn't make sense!?!");
        }
        */

        feature_converter ff = cds.front(), fb = cds.back();

        ////// all we really want is peptide transcript offset and peptide length!?!?
        ///// all else follows from that?!?
        //// offset is trivial - length is then just full transcript length minus other offset!?!?

        // _pstart = cds.front().gstart();
        // _pend = cds.back().gend();
        // sorted 5'->3' 
        _gpstart = ff.gstart();
        _gpend = fb.gend();

        /// can clearly use convertg2t - but it there really any point - i.e. if we have corresponding first/last exons?!?
        cout << "peptide 5' start is "<< _gpstart << "\n";
        cout << "peptide 3' end is "<< _gpend << "\n";

        printexons();

        for (auto eit=exons.begin() ; eit!=exons.end(); eit++) { 
            ///// do 5'
            if (ff.gstart()>=eit->gstart()&&ff.gend()<=eit->gend()) {
                cout << "checkin "<<eit->gstart()<<"-"<<eit->gend()<<"\n";
              // always the 5' start - will be transcript start if strand is positive?!?
              if (_strand==1) _peptide_transcript_offset = eit->gend()-_gpstart+1;
              // if (_strand) _peptide_transcript_length=_gpstart-eit->gstart()+1;
              if (_strand==0) _peptide_transcript_length= eit->tend+_gpstart-eit->gend()-1;
                
            }
            ///// do 3'
            if (fb.gstart()>=eit->gstart()&&fb.gend()<=eit->gend()) {
              // if (_strand) _peptide_transcript_offset = _gpstart-eit->gstart()+1;
              if (_strand) _peptide_transcript_length= eit->tend+_gpend-eit->gend()-1;
              if (!_strand) _peptide_transcript_offset = eit->gend()-_gpend+1;
              // if (_strand) _peptide_transcript_length=eit->tend-(eit->gend()-_gpend+1);
            }
              //_peptide_transcript_length = _strand==1? eit->gend()-_gpend+1:_gpstart-eit->gstart()+1;
        }
        cout << "peptide transcript offset is "<< _peptide_transcript_offset << "\n";
        cout << "peptide transcript length is "<< _peptide_transcript_length << "\n";
        if (_peptide_transcript_offset==0) throw runtime_error("unable to find corresponding exon for most 5' cds - this gff model doesn't make sense!?!");
        if (_peptide_transcript_length==0) throw runtime_error("unable to find corresponding exon for most 3' cds - this gff model doesn't make sense!?!");

    }

};


class iter {
     // this is not the way to do it - i.e. wanna i.e. too many dereferences?!?
// have it store forward or reverse iterator?!?
};

inline void doit(feature_converter& f,int &i, int& last_te) {
}

template<typename T> inline void doit(T& it, T& eit) {
    int i=0,last_te=0;
    for ( ; it!=eit; it++,i++) {
    // for ( ; it!=eit; it++,i++) {
         cout <<"i="<<i<< " " << it->gstart() << "\n";

         // if (i==0) 
         if(i==0) {
             // already 0 initialised tstart?!?!
             it->tstart=1;
             // use ternaries e.g. it->tend = strand ? it->gend-it->gstart+1 : ...?!?
             it->tend=it->gend()-it->gstart()+1;
         } else {
             it->tstart=last_te+1;
             it->tend=it->gend()-it->gstart()+it->tstart; 

         }
        
        last_te=it->tend;
    }
}

int main () {

// cout<<"here="<<sizeof(bool)<<"\nhere="<<sizeof(signed char)<<"\n";


cout <<"\n\n\n";
    transcript trans(1);
    trans.addexon(feature_converter(14500,15000));
    trans.addexon(feature_converter(12000,14000));
    trans.addexon(feature_converter(10000,11000));
    trans.addcds(feature_converter(14500,14600));
    trans.addcds(feature_converter(12000,14000));
    trans.addcds(feature_converter(10500,11000));

    int pos =1254;
    trans.printexons();
    cout << "G2T : is " << pos << " = " <<trans.convertt2g(pos)<<"\n";;
    trans.printexons();
    pos = trans.convertt2g(pos);
    cout << "AND REVERSE " << pos << " = " <<trans.convertg2t(pos)<<"\n";;
    trans.printexons();
    cout << "and intron " << pos << " = " <<trans.convertg2t(11500)<<"\n";;
    // for_each(feats.begin(),feats.end(), [](feature_converter& f) { cout << f.gstart <<" "<<f.gend << "\n"; });

    transcript trans1(0);
    trans1.addexon(feature_converter(14500,15000));
    trans1.addexon(feature_converter(12000,14000));
    trans1.addexon(feature_converter(10000,11000));
    trans1.addcds(feature_converter(14500,14600));
    trans1.addcds(feature_converter(12000,14000));
    trans1.addcds(feature_converter(10500,11000));

    pos =1254;
    trans1.printexons();
    cout << "G2T : is " << pos << " = " <<trans1.convertt2g(pos)<<"\n";;
    trans1.printexons();
    pos = trans1.convertt2g(pos);
    cout << "AND REVERSE " << pos << " = " <<trans1.convertg2t(pos)<<"\n";;
    cout << "and intron " << pos << " = " <<trans1.convertg2t(11500)<<"\n";;
    trans1.printexons();

    // for_each(feats.begin(),feats.end(), [](feature_converter& f) { 

   //  if (1) while ((first!=last)&&(first!=--last)) swap (*first++,*last);
   // this really shouldn't be necessary - need to?!?
//    if(strand) reverse(feats.begin(),feats.end());
    // clearly this is strand dependent - reverse order and use of start/end

    // std::vector<feature_converters>::iterator sit, eit;
    return 0;

// don't want to reverse the elements?!?
// std::vector<feature_converter>::iterator it, eit;
// deque and vector are random access?!?

//////////////////// perhaps use type erasure to allow iterators to be used generaically?!?

}


