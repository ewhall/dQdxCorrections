// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dEdxConversionDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Cal_Test.h"
#include "Calibration_File.h"
#include "ConvertdEdx.h"
#include "DEdx_Test.h"
#include "MakeDEdxTracks.h"
#include "Projection2D.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLCal_Test_Dictionary();
   static void larlitecLcLCal_Test_TClassManip(TClass*);
   static void *new_larlitecLcLCal_Test(void *p = 0);
   static void *newArray_larlitecLcLCal_Test(Long_t size, void *p);
   static void delete_larlitecLcLCal_Test(void *p);
   static void deleteArray_larlitecLcLCal_Test(void *p);
   static void destruct_larlitecLcLCal_Test(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::Cal_Test*)
   {
      ::larlite::Cal_Test *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::Cal_Test));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::Cal_Test", "Cal_Test.h", 26,
                  typeid(::larlite::Cal_Test), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecLcLCal_Test_Dictionary, isa_proxy, 0,
                  sizeof(::larlite::Cal_Test) );
      instance.SetNew(&new_larlitecLcLCal_Test);
      instance.SetNewArray(&newArray_larlitecLcLCal_Test);
      instance.SetDelete(&delete_larlitecLcLCal_Test);
      instance.SetDeleteArray(&deleteArray_larlitecLcLCal_Test);
      instance.SetDestructor(&destruct_larlitecLcLCal_Test);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::Cal_Test*)
   {
      return GenerateInitInstanceLocal((::larlite::Cal_Test*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::larlite::Cal_Test*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLCal_Test_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::Cal_Test*)0x0)->GetClass();
      larlitecLcLCal_Test_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLCal_Test_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLCalibration_File_Dictionary();
   static void larlitecLcLCalibration_File_TClassManip(TClass*);
   static void *new_larlitecLcLCalibration_File(void *p = 0);
   static void *newArray_larlitecLcLCalibration_File(Long_t size, void *p);
   static void delete_larlitecLcLCalibration_File(void *p);
   static void deleteArray_larlitecLcLCalibration_File(void *p);
   static void destruct_larlitecLcLCalibration_File(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::Calibration_File*)
   {
      ::larlite::Calibration_File *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::Calibration_File));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::Calibration_File", "Calibration_File.h", 43,
                  typeid(::larlite::Calibration_File), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecLcLCalibration_File_Dictionary, isa_proxy, 0,
                  sizeof(::larlite::Calibration_File) );
      instance.SetNew(&new_larlitecLcLCalibration_File);
      instance.SetNewArray(&newArray_larlitecLcLCalibration_File);
      instance.SetDelete(&delete_larlitecLcLCalibration_File);
      instance.SetDeleteArray(&deleteArray_larlitecLcLCalibration_File);
      instance.SetDestructor(&destruct_larlitecLcLCalibration_File);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::Calibration_File*)
   {
      return GenerateInitInstanceLocal((::larlite::Calibration_File*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::larlite::Calibration_File*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLCalibration_File_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::Calibration_File*)0x0)->GetClass();
      larlitecLcLCalibration_File_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLCalibration_File_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLConvertdEdx_Dictionary();
   static void larlitecLcLConvertdEdx_TClassManip(TClass*);
   static void *new_larlitecLcLConvertdEdx(void *p = 0);
   static void *newArray_larlitecLcLConvertdEdx(Long_t size, void *p);
   static void delete_larlitecLcLConvertdEdx(void *p);
   static void deleteArray_larlitecLcLConvertdEdx(void *p);
   static void destruct_larlitecLcLConvertdEdx(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::ConvertdEdx*)
   {
      ::larlite::ConvertdEdx *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::ConvertdEdx));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::ConvertdEdx", "ConvertdEdx.h", 42,
                  typeid(::larlite::ConvertdEdx), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecLcLConvertdEdx_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::ConvertdEdx) );
      instance.SetNew(&new_larlitecLcLConvertdEdx);
      instance.SetNewArray(&newArray_larlitecLcLConvertdEdx);
      instance.SetDelete(&delete_larlitecLcLConvertdEdx);
      instance.SetDeleteArray(&deleteArray_larlitecLcLConvertdEdx);
      instance.SetDestructor(&destruct_larlitecLcLConvertdEdx);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::ConvertdEdx*)
   {
      return GenerateInitInstanceLocal((::larlite::ConvertdEdx*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::larlite::ConvertdEdx*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLConvertdEdx_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::ConvertdEdx*)0x0)->GetClass();
      larlitecLcLConvertdEdx_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLConvertdEdx_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLDEdx_Test_Dictionary();
   static void larlitecLcLDEdx_Test_TClassManip(TClass*);
   static void *new_larlitecLcLDEdx_Test(void *p = 0);
   static void *newArray_larlitecLcLDEdx_Test(Long_t size, void *p);
   static void delete_larlitecLcLDEdx_Test(void *p);
   static void deleteArray_larlitecLcLDEdx_Test(void *p);
   static void destruct_larlitecLcLDEdx_Test(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::DEdx_Test*)
   {
      ::larlite::DEdx_Test *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::DEdx_Test));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::DEdx_Test", "DEdx_Test.h", 28,
                  typeid(::larlite::DEdx_Test), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecLcLDEdx_Test_Dictionary, isa_proxy, 0,
                  sizeof(::larlite::DEdx_Test) );
      instance.SetNew(&new_larlitecLcLDEdx_Test);
      instance.SetNewArray(&newArray_larlitecLcLDEdx_Test);
      instance.SetDelete(&delete_larlitecLcLDEdx_Test);
      instance.SetDeleteArray(&deleteArray_larlitecLcLDEdx_Test);
      instance.SetDestructor(&destruct_larlitecLcLDEdx_Test);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::DEdx_Test*)
   {
      return GenerateInitInstanceLocal((::larlite::DEdx_Test*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::larlite::DEdx_Test*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLDEdx_Test_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::DEdx_Test*)0x0)->GetClass();
      larlitecLcLDEdx_Test_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLDEdx_Test_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLMakeDEdxTracks_Dictionary();
   static void larlitecLcLMakeDEdxTracks_TClassManip(TClass*);
   static void *new_larlitecLcLMakeDEdxTracks(void *p = 0);
   static void *newArray_larlitecLcLMakeDEdxTracks(Long_t size, void *p);
   static void delete_larlitecLcLMakeDEdxTracks(void *p);
   static void deleteArray_larlitecLcLMakeDEdxTracks(void *p);
   static void destruct_larlitecLcLMakeDEdxTracks(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::MakeDEdxTracks*)
   {
      ::larlite::MakeDEdxTracks *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::MakeDEdxTracks));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::MakeDEdxTracks", "MakeDEdxTracks.h", 26,
                  typeid(::larlite::MakeDEdxTracks), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecLcLMakeDEdxTracks_Dictionary, isa_proxy, 0,
                  sizeof(::larlite::MakeDEdxTracks) );
      instance.SetNew(&new_larlitecLcLMakeDEdxTracks);
      instance.SetNewArray(&newArray_larlitecLcLMakeDEdxTracks);
      instance.SetDelete(&delete_larlitecLcLMakeDEdxTracks);
      instance.SetDeleteArray(&deleteArray_larlitecLcLMakeDEdxTracks);
      instance.SetDestructor(&destruct_larlitecLcLMakeDEdxTracks);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::MakeDEdxTracks*)
   {
      return GenerateInitInstanceLocal((::larlite::MakeDEdxTracks*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::larlite::MakeDEdxTracks*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLMakeDEdxTracks_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::MakeDEdxTracks*)0x0)->GetClass();
      larlitecLcLMakeDEdxTracks_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLMakeDEdxTracks_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLCal_Test(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::Cal_Test : new ::larlite::Cal_Test;
   }
   static void *newArray_larlitecLcLCal_Test(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::Cal_Test[nElements] : new ::larlite::Cal_Test[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLCal_Test(void *p) {
      delete ((::larlite::Cal_Test*)p);
   }
   static void deleteArray_larlitecLcLCal_Test(void *p) {
      delete [] ((::larlite::Cal_Test*)p);
   }
   static void destruct_larlitecLcLCal_Test(void *p) {
      typedef ::larlite::Cal_Test current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::Cal_Test

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLCalibration_File(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::Calibration_File : new ::larlite::Calibration_File;
   }
   static void *newArray_larlitecLcLCalibration_File(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::Calibration_File[nElements] : new ::larlite::Calibration_File[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLCalibration_File(void *p) {
      delete ((::larlite::Calibration_File*)p);
   }
   static void deleteArray_larlitecLcLCalibration_File(void *p) {
      delete [] ((::larlite::Calibration_File*)p);
   }
   static void destruct_larlitecLcLCalibration_File(void *p) {
      typedef ::larlite::Calibration_File current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::Calibration_File

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLConvertdEdx(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::ConvertdEdx : new ::larlite::ConvertdEdx;
   }
   static void *newArray_larlitecLcLConvertdEdx(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::ConvertdEdx[nElements] : new ::larlite::ConvertdEdx[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLConvertdEdx(void *p) {
      delete ((::larlite::ConvertdEdx*)p);
   }
   static void deleteArray_larlitecLcLConvertdEdx(void *p) {
      delete [] ((::larlite::ConvertdEdx*)p);
   }
   static void destruct_larlitecLcLConvertdEdx(void *p) {
      typedef ::larlite::ConvertdEdx current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::ConvertdEdx

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLDEdx_Test(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::DEdx_Test : new ::larlite::DEdx_Test;
   }
   static void *newArray_larlitecLcLDEdx_Test(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::DEdx_Test[nElements] : new ::larlite::DEdx_Test[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLDEdx_Test(void *p) {
      delete ((::larlite::DEdx_Test*)p);
   }
   static void deleteArray_larlitecLcLDEdx_Test(void *p) {
      delete [] ((::larlite::DEdx_Test*)p);
   }
   static void destruct_larlitecLcLDEdx_Test(void *p) {
      typedef ::larlite::DEdx_Test current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::DEdx_Test

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLMakeDEdxTracks(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::MakeDEdxTracks : new ::larlite::MakeDEdxTracks;
   }
   static void *newArray_larlitecLcLMakeDEdxTracks(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::MakeDEdxTracks[nElements] : new ::larlite::MakeDEdxTracks[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLMakeDEdxTracks(void *p) {
      delete ((::larlite::MakeDEdxTracks*)p);
   }
   static void deleteArray_larlitecLcLMakeDEdxTracks(void *p) {
      delete [] ((::larlite::MakeDEdxTracks*)p);
   }
   static void destruct_larlitecLcLMakeDEdxTracks(void *p) {
      typedef ::larlite::MakeDEdxTracks current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::MakeDEdxTracks

namespace ROOT {
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary();
   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p);
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<double> >*)
   {
      vector<vector<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<double> >", -2, "vector", 447,
                  typeid(vector<vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<double> >) );
      instance.SetNew(&new_vectorlEvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<double> >*)0x0)->GetClass();
      vectorlEvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> > : new vector<vector<double> >;
   }
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> >[nElements] : new vector<vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete ((vector<vector<double> >*)p);
   }
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete [] ((vector<vector<double> >*)p);
   }
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p) {
      typedef vector<vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<double> >

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 447,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_libdEdxConversion_Impl() {
    static const char* headers[] = {
"Cal_Test.h",
"Calibration_File.h",
"ConvertdEdx.h",
"DEdx_Test.h",
"MakeDEdxTracks.h",
"Projection2D.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/Cellar/root/6.18.04/include/root",
"/Users/mwhall/dllee_unified/larlite/core",
"/Users/mwhall/dllee_unified/LArCV/build/include",
"/usr/local/opt/opencv@3/include",
"/usr/local/Cellar/python/3.7.5/Frameworks/Python.framework/Versions/3.7/include/python3.7m",
"/usr/local/Cellar/python/3.7.5/Frameworks/Python.framework/Versions/3.7/include/python3.7m",
"/usr/local/lib/python3.7/site-packages/numpy/core/include",
"/Users/mwhall/dllee_unified/larlite/core",
"/Users/mwhall/dllee_unified/larlite/UserDev",
"/Users/mwhall/dllee_unified/LArCV/app/UBWireTool",
"/usr/local/Cellar/root/6.18.04/include/root",
"/Users/mwhall/dllee_unified/larlite/UserDev/RecoTool/dEdxConversion/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libdEdxConversion dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$Cal_Test.h")))  Cal_Test;}
namespace larlite{class __attribute__((annotate("$clingAutoload$Calibration_File.h")))  Calibration_File;}
namespace larlite{class __attribute__((annotate("$clingAutoload$ConvertdEdx.h")))  ConvertdEdx;}
namespace larlite{class __attribute__((annotate("$clingAutoload$DEdx_Test.h")))  DEdx_Test;}
namespace larlite{class __attribute__((annotate("$clingAutoload$MakeDEdxTracks.h")))  MakeDEdxTracks;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libdEdxConversion dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Cal_Test.h"
#include "Calibration_File.h"
#include "ConvertdEdx.h"
#include "DEdx_Test.h"
#include "MakeDEdxTracks.h"
#include "Projection2D.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::Cal_Test", payloadCode, "@",
"larlite::Calibration_File", payloadCode, "@",
"larlite::ConvertdEdx", payloadCode, "@",
"larlite::DEdx_Test", payloadCode, "@",
"larlite::MakeDEdxTracks", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libdEdxConversion",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libdEdxConversion_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libdEdxConversion_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libdEdxConversion() {
  TriggerDictionaryInitialization_libdEdxConversion_Impl();
}
