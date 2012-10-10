#include <iostream>
#include <bitset>
using namespace std;

#define APOLLO_SCF_NAMES                            (1<<0) // (0x1<<0) 0x1
#define NAMES_HAVE_SPACES                           (1<<1)
#define LINES_WO_9COLS                              (1<<2)
#define LINE_ENDINGS                                (1<<3)
#define ILEGAL_FEAT_TYPES                           (1<<4)

#define EXCESS_GENE_CONSISTENCY_PROB                (1<<5)
#define EXCESS_mRNA_REL_GENE_CONSISTENCY_PROB       (1<<6)
#define EXCESS_mRNA_REL_CDS_EXON_CONSISTENCY_PROB   (1<<7)
#define EXCESS_CDS_EXON_CONSISTENCY_PROB            (1<<8)
#define NON_PERMITTED_BIOTYPES                      (1<<9)
#define ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE       (1<<10)
#define PARENT_WITHOUT_ID_NOT_CDS_EXON              (1<<11)

#define UNKNOWN_SCAFFOLD                            (1<<12)

#define NEGATIVE_COORDINATES                        (1<<13)
#define NON_PRINTING_X0D       (1<<14)
#define BLANK_LINES          (1<<15)

#define NON_UNIQUE_ID        (1<<16)

// no_CDS_without_pseudogene = false; // only processing protein coding via cap
#define PSEUDOGENES_PRESENT (1<<17)
#define CDS_PRESENT (1<<18)

#define EMBL_FORMAT (1<<19)
#define GFF_FASTA_HEADER (1<<20)
#define FASTA_HEADER (1<<21)

#define PARTIAL_MODEL (1<<22)

static long long PROBLEM = EXCESS_GENE_CONSISTENCY_PROB | EXCESS_mRNA_REL_GENE_CONSISTENCY_PROB 
    | EXCESS_mRNA_REL_CDS_EXON_CONSISTENCY_PROB | EXCESS_CDS_EXON_CONSISTENCY_PROB | ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE 
    | PARENT_WITHOUT_ID_NOT_CDS_EXON | UNKNOWN_SCAFFOLD|NON_PERMITTED_BIOTYPES | NEGATIVE_COORDINATES 
    | NON_UNIQUE_ID | EMBL_FORMAT | GFF_FASTA_HEADER | FASTA_HEADER | PARTIAL_MODEL;

// PROBLEM ~= PROBLEM;
// PROBLEM = ~PROBLEM;

inline bool blah() { return true; }

int main() {


for (int i=100 ; i>0 ; i--) {

    cout << "hey ";

    if(blah()) continue;

    cout << "hoo \n";
}

return 0;


    long long bm = 0;


    bm = CDS_PRESENT | BLANK_LINES | PSEUDOGENES_PRESENT;
    bool stop = bm & PROBLEM;
    cout << "and with mask = "<< hex << bitset<32>(bm) << " = stop = " << stop << "\n";


    long long no_cds_and_no_pseudogenes = CDS_PRESENT | PSEUDOGENES_PRESENT;

    cout << "and with mask = "<< hex << bitset<32>(no_cds_and_no_pseudogenes) << " = stop = " << stop << "\n";
    return 0;

    for (long long i = 0x1, c = 0 ; i <= 0x800000; i <<=1, ++c ) {
        int x = bm&i;
        cout << "i="<< hex << i << dec << " - " << c << " - " << bitset<16>(i) << " and bm &=i " << x<< "\n";
    }

    bm |= EXCESS_GENE_CONSISTENCY_PROB;

    for (long long i = 0x1, c = 0 ; i <= 0x800000; i <<=1, ++c ) {
        int x = bm&i;
        cout << "i="<< hex << i << dec << " - " << c << " - " << bitset<16>(i) << " and bm &=i " << x<< "\n";
    }

    bm = 0x80;

//     if (ID_WITHOUT_PARENT_NOT_GENE_PSEUDOGENE

        cout << "\n\nprob="<< hex << PROBLEM << dec << " - " << bitset<32>(PROBLEM) << " and bm &=i\n";

 //        PROBLEM = ~PROBLEM;
        cout << "\n\nprob="<< hex << PROBLEM << dec << " - " << bitset<32>(PROBLEM) << " and bm &=i\n\n";

    for (long long i = 0x1, c = 0 ; i <= 0x8000000; i <<=1, ++c ) {
        bool stop = i & PROBLEM;
        cout << "i="<< hex << i << dec << " - " << c << " - " << bitset<32>(i) << " and bm &=i " << stop << "\n";
    }


         PROBLEM = ~PROBLEM;
    for (long long i = 0x1, c = 0 ; i <= 0x8000000; i <<=1, ++c ) {
        bool proceed = i & PROBLEM;
        cout << "i="<< hex << i << dec << " - " << c << " - " << bitset<32>(i) << " and bm &=i " << proceed << "\n";
    }
}
