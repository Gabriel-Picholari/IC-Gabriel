#ifndef MyQuark_h
#define MyQuark_h

class MyQuark;

class MyQuark: public TObject
{
 public: MyQuark() :
    fPt(0), fEta(0), fPhi(0), fMass(0), fPx(0), fPy(0), fPz(0), fE(0), fCharge(0), fPdg(0), lqpT(0), lqEta(0), lqPhi(0) {;}

 public:

    Double32_t fPt;
    Double32_t fEta;
    Double32_t fPhi;
    Double32_t fMass;
    Double32_t fPx;
    Double32_t fPy;
    Double32_t fPz;
    Double32_t fE;
    Short_t fPdg;
    Double32_t fCharge;

    // Last quark variables

    Double32_t lqpT;
    Double32_t lqEta;
    Double32_t lqPhi;

  ClassDef(MyQuark, 1)
};

#endif