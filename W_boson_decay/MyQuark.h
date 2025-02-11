#ifndef MyQuark_h
#define MyQuark_h

class MyQuark;

class MyQuark: public TObject
{
 public: MyQuark() :
    qPdg(0), qpT(0), qEta(0), qPhi(0) {;}

 public:

    Int_t         qPdg;
    Float_t       qpT;
    Float_t       qEta;
    Float_t       qPhi;

  ClassDef(MyQuark, 1)
};

#endif
