#ifndef FEATURE_H
#define FEATURE_H

typedef unsigned int uint;
typedef unsigned char uchar;
// typedef unsigned char byte;

class feature_min { 

    // make them all private?!?!?
    uint _gstart;
    uint _gend;

  public:

 //   feature_min() = delete; i.e. for conversion ctor - perhaps should use operator_feature_converter?!?
    feature_min(uint a, uint b) : _gstart(a), _gend(b) {} // , tstart(0), tend(0) {} 
    virtual ~feature_min() {} // undefined not to have virtual dtor in base!?!

    uint gstart() const { return _gstart; } // doesn't need setter!?!?
    uint gend() const { return _gend; } 

};

#endif
