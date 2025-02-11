#ifndef MyJet_h
#define MyJet_h

class MyJet;

class MyJet: public TObject
{
 public: MyJet() :
  TObject(), jetPt(0), jetEta(0), jetPhi(0), jetMass(0), jetPx(0), jetPy(0), jetPz(0), jetE(0), nConstituent(0), pTLeadConstituent(0){;}

 public:
    Double32_t    jetPt;
    Double32_t    jetEta;
    Double32_t    jetPhi;
    Double32_t    jetMass;
    Double32_t    jetPx;
    Double32_t    jetPy;
    Double32_t    jetPz;
    Double32_t    jetE;
    Int_t         nConstituent;
    Double32_t    pTLeadConstituent;

  ClassDef(MyJet, 1)
};

#endif