#ifndef MyPart_h
#define MyPart_h

class MyPart;

class MyPart: public TObject
{
 public: MyPart() :
  TObject(), fPt(0), fEta(0), fPhi(0), fMass(0), fPx(0), fPy(0), fPz(0), fE(0), fCharge(0), fFirstMother(0), fPdg(0){;}

 public:
  Double32_t    fPt;               //[0,0,16] pt
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Double32_t    fMass;             //[0,0,16] mass
  Double32_t    fPx;               //[0,0,16] px
  Double32_t    fPy;               //[0,0,16] py
  Double32_t    fPz;               //[0,0,16] pz
  Double32_t    fE;                //[0,0,16] energy
  Double32_t    fCharge;
  Int_t         fFirstMother;
  Short_t       fPdg;              //pdg code

  ClassDef(MyPart, 1)
};

#endif