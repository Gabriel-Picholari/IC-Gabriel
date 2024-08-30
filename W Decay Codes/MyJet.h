#ifndef MyJet_h
#define MyJet_h

class MyJet;

class MyJet: public TObject
{
 public: MyJet() :
  TObject(), fPt(0), fEta(0), fPhi(0), fMass(0), fPx(0), fPy(0), fPz(0), fE(0) {;}

 public:

    Float_t    fPt;
    Float_t    fEta;
    Float_t    fPhi;
    Float_t    fMass;
    Float_t    fPx;
    Float_t    fPy;
    Float_t    fPz;
    Float_t    fE;
    Float_t    fnConst;
    Float_t    pT_LeadConst;

  ClassDef(MyJet, 1)
};

#endif
