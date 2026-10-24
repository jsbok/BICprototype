// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__rootIO
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include/RootInterface.h"
#include "/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include/koBICInterface.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *tuplelEfloatcOfloatcOfloatgR_Dictionary();
   static void tuplelEfloatcOfloatcOfloatgR_TClassManip(TClass*);
   static void *new_tuplelEfloatcOfloatcOfloatgR(void *p = nullptr);
   static void *newArray_tuplelEfloatcOfloatcOfloatgR(Long_t size, void *p);
   static void delete_tuplelEfloatcOfloatcOfloatgR(void *p);
   static void deleteArray_tuplelEfloatcOfloatcOfloatgR(void *p);
   static void destruct_tuplelEfloatcOfloatcOfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::tuple<float,float,float>*)
   {
      ::tuple<float,float,float> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::tuple<float,float,float>));
      static ::ROOT::TGenericClassInfo 
         instance("tuple<float,float,float>", "tuple", 609,
                  typeid(::tuple<float,float,float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &tuplelEfloatcOfloatcOfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(::tuple<float,float,float>) );
      instance.SetNew(&new_tuplelEfloatcOfloatcOfloatgR);
      instance.SetNewArray(&newArray_tuplelEfloatcOfloatcOfloatgR);
      instance.SetDelete(&delete_tuplelEfloatcOfloatcOfloatgR);
      instance.SetDeleteArray(&deleteArray_tuplelEfloatcOfloatcOfloatgR);
      instance.SetDestructor(&destruct_tuplelEfloatcOfloatcOfloatgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("tuple<float,float,float>","koBICInterface::threeVector"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("tuple<float,float,float>","std::tuple<float, float, float>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::tuple<float,float,float>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *tuplelEfloatcOfloatcOfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::tuple<float,float,float>*>(nullptr))->GetClass();
      tuplelEfloatcOfloatcOfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void tuplelEfloatcOfloatcOfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *koBICInterfacecLcLkoBICSiPMData_Dictionary();
   static void koBICInterfacecLcLkoBICSiPMData_TClassManip(TClass*);
   static void *new_koBICInterfacecLcLkoBICSiPMData(void *p = nullptr);
   static void *newArray_koBICInterfacecLcLkoBICSiPMData(Long_t size, void *p);
   static void delete_koBICInterfacecLcLkoBICSiPMData(void *p);
   static void deleteArray_koBICInterfacecLcLkoBICSiPMData(void *p);
   static void destruct_koBICInterfacecLcLkoBICSiPMData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::koBICInterface::koBICSiPMData*)
   {
      ::koBICInterface::koBICSiPMData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::koBICInterface::koBICSiPMData));
      static ::ROOT::TGenericClassInfo 
         instance("koBICInterface::koBICSiPMData", "koBICInterface.h", 28,
                  typeid(::koBICInterface::koBICSiPMData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &koBICInterfacecLcLkoBICSiPMData_Dictionary, isa_proxy, 4,
                  sizeof(::koBICInterface::koBICSiPMData) );
      instance.SetNew(&new_koBICInterfacecLcLkoBICSiPMData);
      instance.SetNewArray(&newArray_koBICInterfacecLcLkoBICSiPMData);
      instance.SetDelete(&delete_koBICInterfacecLcLkoBICSiPMData);
      instance.SetDeleteArray(&deleteArray_koBICInterfacecLcLkoBICSiPMData);
      instance.SetDestructor(&destruct_koBICInterfacecLcLkoBICSiPMData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::koBICInterface::koBICSiPMData*)
   {
      return GenerateInitInstanceLocal(static_cast<::koBICInterface::koBICSiPMData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICSiPMData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *koBICInterfacecLcLkoBICSiPMData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICSiPMData*>(nullptr))->GetClass();
      koBICInterfacecLcLkoBICSiPMData_TClassManip(theClass);
   return theClass;
   }

   static void koBICInterfacecLcLkoBICSiPMData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *koBICInterfacecLcLkoBICTowerData_Dictionary();
   static void koBICInterfacecLcLkoBICTowerData_TClassManip(TClass*);
   static void *new_koBICInterfacecLcLkoBICTowerData(void *p = nullptr);
   static void *newArray_koBICInterfacecLcLkoBICTowerData(Long_t size, void *p);
   static void delete_koBICInterfacecLcLkoBICTowerData(void *p);
   static void deleteArray_koBICInterfacecLcLkoBICTowerData(void *p);
   static void destruct_koBICInterfacecLcLkoBICTowerData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::koBICInterface::koBICTowerData*)
   {
      ::koBICInterface::koBICTowerData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::koBICInterface::koBICTowerData));
      static ::ROOT::TGenericClassInfo 
         instance("koBICInterface::koBICTowerData", "koBICInterface.h", 43,
                  typeid(::koBICInterface::koBICTowerData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &koBICInterfacecLcLkoBICTowerData_Dictionary, isa_proxy, 4,
                  sizeof(::koBICInterface::koBICTowerData) );
      instance.SetNew(&new_koBICInterfacecLcLkoBICTowerData);
      instance.SetNewArray(&newArray_koBICInterfacecLcLkoBICTowerData);
      instance.SetDelete(&delete_koBICInterfacecLcLkoBICTowerData);
      instance.SetDeleteArray(&deleteArray_koBICInterfacecLcLkoBICTowerData);
      instance.SetDestructor(&destruct_koBICInterfacecLcLkoBICTowerData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::koBICInterface::koBICTowerData*)
   {
      return GenerateInitInstanceLocal(static_cast<::koBICInterface::koBICTowerData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICTowerData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *koBICInterfacecLcLkoBICTowerData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICTowerData*>(nullptr))->GetClass();
      koBICInterfacecLcLkoBICTowerData_TClassManip(theClass);
   return theClass;
   }

   static void koBICInterfacecLcLkoBICTowerData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *koBICInterfacecLcLkoBICEdepData_Dictionary();
   static void koBICInterfacecLcLkoBICEdepData_TClassManip(TClass*);
   static void *new_koBICInterfacecLcLkoBICEdepData(void *p = nullptr);
   static void *newArray_koBICInterfacecLcLkoBICEdepData(Long_t size, void *p);
   static void delete_koBICInterfacecLcLkoBICEdepData(void *p);
   static void deleteArray_koBICInterfacecLcLkoBICEdepData(void *p);
   static void destruct_koBICInterfacecLcLkoBICEdepData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::koBICInterface::koBICEdepData*)
   {
      ::koBICInterface::koBICEdepData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::koBICInterface::koBICEdepData));
      static ::ROOT::TGenericClassInfo 
         instance("koBICInterface::koBICEdepData", "koBICInterface.h", 53,
                  typeid(::koBICInterface::koBICEdepData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &koBICInterfacecLcLkoBICEdepData_Dictionary, isa_proxy, 4,
                  sizeof(::koBICInterface::koBICEdepData) );
      instance.SetNew(&new_koBICInterfacecLcLkoBICEdepData);
      instance.SetNewArray(&newArray_koBICInterfacecLcLkoBICEdepData);
      instance.SetDelete(&delete_koBICInterfacecLcLkoBICEdepData);
      instance.SetDeleteArray(&deleteArray_koBICInterfacecLcLkoBICEdepData);
      instance.SetDestructor(&destruct_koBICInterfacecLcLkoBICEdepData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::koBICInterface::koBICEdepData*)
   {
      return GenerateInitInstanceLocal(static_cast<::koBICInterface::koBICEdepData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICEdepData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *koBICInterfacecLcLkoBICEdepData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICEdepData*>(nullptr))->GetClass();
      koBICInterfacecLcLkoBICEdepData_TClassManip(theClass);
   return theClass;
   }

   static void koBICInterfacecLcLkoBICEdepData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *koBICInterfacecLcLkoBICLeakageData_Dictionary();
   static void koBICInterfacecLcLkoBICLeakageData_TClassManip(TClass*);
   static void *new_koBICInterfacecLcLkoBICLeakageData(void *p = nullptr);
   static void *newArray_koBICInterfacecLcLkoBICLeakageData(Long_t size, void *p);
   static void delete_koBICInterfacecLcLkoBICLeakageData(void *p);
   static void deleteArray_koBICInterfacecLcLkoBICLeakageData(void *p);
   static void destruct_koBICInterfacecLcLkoBICLeakageData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::koBICInterface::koBICLeakageData*)
   {
      ::koBICInterface::koBICLeakageData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::koBICInterface::koBICLeakageData));
      static ::ROOT::TGenericClassInfo 
         instance("koBICInterface::koBICLeakageData", "koBICInterface.h", 65,
                  typeid(::koBICInterface::koBICLeakageData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &koBICInterfacecLcLkoBICLeakageData_Dictionary, isa_proxy, 4,
                  sizeof(::koBICInterface::koBICLeakageData) );
      instance.SetNew(&new_koBICInterfacecLcLkoBICLeakageData);
      instance.SetNewArray(&newArray_koBICInterfacecLcLkoBICLeakageData);
      instance.SetDelete(&delete_koBICInterfacecLcLkoBICLeakageData);
      instance.SetDeleteArray(&deleteArray_koBICInterfacecLcLkoBICLeakageData);
      instance.SetDestructor(&destruct_koBICInterfacecLcLkoBICLeakageData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::koBICInterface::koBICLeakageData*)
   {
      return GenerateInitInstanceLocal(static_cast<::koBICInterface::koBICLeakageData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICLeakageData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *koBICInterfacecLcLkoBICLeakageData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICLeakageData*>(nullptr))->GetClass();
      koBICInterfacecLcLkoBICLeakageData_TClassManip(theClass);
   return theClass;
   }

   static void koBICInterfacecLcLkoBICLeakageData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *koBICInterfacecLcLkoBICGenData_Dictionary();
   static void koBICInterfacecLcLkoBICGenData_TClassManip(TClass*);
   static void *new_koBICInterfacecLcLkoBICGenData(void *p = nullptr);
   static void *newArray_koBICInterfacecLcLkoBICGenData(Long_t size, void *p);
   static void delete_koBICInterfacecLcLkoBICGenData(void *p);
   static void deleteArray_koBICInterfacecLcLkoBICGenData(void *p);
   static void destruct_koBICInterfacecLcLkoBICGenData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::koBICInterface::koBICGenData*)
   {
      ::koBICInterface::koBICGenData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::koBICInterface::koBICGenData));
      static ::ROOT::TGenericClassInfo 
         instance("koBICInterface::koBICGenData", "koBICInterface.h", 83,
                  typeid(::koBICInterface::koBICGenData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &koBICInterfacecLcLkoBICGenData_Dictionary, isa_proxy, 4,
                  sizeof(::koBICInterface::koBICGenData) );
      instance.SetNew(&new_koBICInterfacecLcLkoBICGenData);
      instance.SetNewArray(&newArray_koBICInterfacecLcLkoBICGenData);
      instance.SetDelete(&delete_koBICInterfacecLcLkoBICGenData);
      instance.SetDeleteArray(&deleteArray_koBICInterfacecLcLkoBICGenData);
      instance.SetDestructor(&destruct_koBICInterfacecLcLkoBICGenData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::koBICInterface::koBICGenData*)
   {
      return GenerateInitInstanceLocal(static_cast<::koBICInterface::koBICGenData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICGenData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *koBICInterfacecLcLkoBICGenData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICGenData*>(nullptr))->GetClass();
      koBICInterfacecLcLkoBICGenData_TClassManip(theClass);
   return theClass;
   }

   static void koBICInterfacecLcLkoBICGenData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *koBICInterfacecLcLkoBICEventData_Dictionary();
   static void koBICInterfacecLcLkoBICEventData_TClassManip(TClass*);
   static void *new_koBICInterfacecLcLkoBICEventData(void *p = nullptr);
   static void *newArray_koBICInterfacecLcLkoBICEventData(Long_t size, void *p);
   static void delete_koBICInterfacecLcLkoBICEventData(void *p);
   static void deleteArray_koBICInterfacecLcLkoBICEventData(void *p);
   static void destruct_koBICInterfacecLcLkoBICEventData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::koBICInterface::koBICEventData*)
   {
      ::koBICInterface::koBICEventData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::koBICInterface::koBICEventData));
      static ::ROOT::TGenericClassInfo 
         instance("koBICInterface::koBICEventData", "koBICInterface.h", 98,
                  typeid(::koBICInterface::koBICEventData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &koBICInterfacecLcLkoBICEventData_Dictionary, isa_proxy, 4,
                  sizeof(::koBICInterface::koBICEventData) );
      instance.SetNew(&new_koBICInterfacecLcLkoBICEventData);
      instance.SetNewArray(&newArray_koBICInterfacecLcLkoBICEventData);
      instance.SetDelete(&delete_koBICInterfacecLcLkoBICEventData);
      instance.SetDeleteArray(&deleteArray_koBICInterfacecLcLkoBICEventData);
      instance.SetDestructor(&destruct_koBICInterfacecLcLkoBICEventData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::koBICInterface::koBICEventData*)
   {
      return GenerateInitInstanceLocal(static_cast<::koBICInterface::koBICEventData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICEventData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *koBICInterfacecLcLkoBICEventData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::koBICInterface::koBICEventData*>(nullptr))->GetClass();
      koBICInterfacecLcLkoBICEventData_TClassManip(theClass);
   return theClass;
   }

   static void koBICInterfacecLcLkoBICEventData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_tuplelEfloatcOfloatcOfloatgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::tuple<float,float,float> : new ::tuple<float,float,float>;
   }
   static void *newArray_tuplelEfloatcOfloatcOfloatgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::tuple<float,float,float>[nElements] : new ::tuple<float,float,float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_tuplelEfloatcOfloatcOfloatgR(void *p) {
      delete (static_cast<::tuple<float,float,float>*>(p));
   }
   static void deleteArray_tuplelEfloatcOfloatcOfloatgR(void *p) {
      delete [] (static_cast<::tuple<float,float,float>*>(p));
   }
   static void destruct_tuplelEfloatcOfloatcOfloatgR(void *p) {
      typedef ::tuple<float,float,float> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::tuple<float,float,float>

namespace ROOT {
   // Wrappers around operator new
   static void *new_koBICInterfacecLcLkoBICSiPMData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICSiPMData : new ::koBICInterface::koBICSiPMData;
   }
   static void *newArray_koBICInterfacecLcLkoBICSiPMData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICSiPMData[nElements] : new ::koBICInterface::koBICSiPMData[nElements];
   }
   // Wrapper around operator delete
   static void delete_koBICInterfacecLcLkoBICSiPMData(void *p) {
      delete (static_cast<::koBICInterface::koBICSiPMData*>(p));
   }
   static void deleteArray_koBICInterfacecLcLkoBICSiPMData(void *p) {
      delete [] (static_cast<::koBICInterface::koBICSiPMData*>(p));
   }
   static void destruct_koBICInterfacecLcLkoBICSiPMData(void *p) {
      typedef ::koBICInterface::koBICSiPMData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::koBICInterface::koBICSiPMData

namespace ROOT {
   // Wrappers around operator new
   static void *new_koBICInterfacecLcLkoBICTowerData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICTowerData : new ::koBICInterface::koBICTowerData;
   }
   static void *newArray_koBICInterfacecLcLkoBICTowerData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICTowerData[nElements] : new ::koBICInterface::koBICTowerData[nElements];
   }
   // Wrapper around operator delete
   static void delete_koBICInterfacecLcLkoBICTowerData(void *p) {
      delete (static_cast<::koBICInterface::koBICTowerData*>(p));
   }
   static void deleteArray_koBICInterfacecLcLkoBICTowerData(void *p) {
      delete [] (static_cast<::koBICInterface::koBICTowerData*>(p));
   }
   static void destruct_koBICInterfacecLcLkoBICTowerData(void *p) {
      typedef ::koBICInterface::koBICTowerData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::koBICInterface::koBICTowerData

namespace ROOT {
   // Wrappers around operator new
   static void *new_koBICInterfacecLcLkoBICEdepData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICEdepData : new ::koBICInterface::koBICEdepData;
   }
   static void *newArray_koBICInterfacecLcLkoBICEdepData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICEdepData[nElements] : new ::koBICInterface::koBICEdepData[nElements];
   }
   // Wrapper around operator delete
   static void delete_koBICInterfacecLcLkoBICEdepData(void *p) {
      delete (static_cast<::koBICInterface::koBICEdepData*>(p));
   }
   static void deleteArray_koBICInterfacecLcLkoBICEdepData(void *p) {
      delete [] (static_cast<::koBICInterface::koBICEdepData*>(p));
   }
   static void destruct_koBICInterfacecLcLkoBICEdepData(void *p) {
      typedef ::koBICInterface::koBICEdepData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::koBICInterface::koBICEdepData

namespace ROOT {
   // Wrappers around operator new
   static void *new_koBICInterfacecLcLkoBICLeakageData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICLeakageData : new ::koBICInterface::koBICLeakageData;
   }
   static void *newArray_koBICInterfacecLcLkoBICLeakageData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICLeakageData[nElements] : new ::koBICInterface::koBICLeakageData[nElements];
   }
   // Wrapper around operator delete
   static void delete_koBICInterfacecLcLkoBICLeakageData(void *p) {
      delete (static_cast<::koBICInterface::koBICLeakageData*>(p));
   }
   static void deleteArray_koBICInterfacecLcLkoBICLeakageData(void *p) {
      delete [] (static_cast<::koBICInterface::koBICLeakageData*>(p));
   }
   static void destruct_koBICInterfacecLcLkoBICLeakageData(void *p) {
      typedef ::koBICInterface::koBICLeakageData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::koBICInterface::koBICLeakageData

namespace ROOT {
   // Wrappers around operator new
   static void *new_koBICInterfacecLcLkoBICGenData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICGenData : new ::koBICInterface::koBICGenData;
   }
   static void *newArray_koBICInterfacecLcLkoBICGenData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICGenData[nElements] : new ::koBICInterface::koBICGenData[nElements];
   }
   // Wrapper around operator delete
   static void delete_koBICInterfacecLcLkoBICGenData(void *p) {
      delete (static_cast<::koBICInterface::koBICGenData*>(p));
   }
   static void deleteArray_koBICInterfacecLcLkoBICGenData(void *p) {
      delete [] (static_cast<::koBICInterface::koBICGenData*>(p));
   }
   static void destruct_koBICInterfacecLcLkoBICGenData(void *p) {
      typedef ::koBICInterface::koBICGenData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::koBICInterface::koBICGenData

namespace ROOT {
   // Wrappers around operator new
   static void *new_koBICInterfacecLcLkoBICEventData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICEventData : new ::koBICInterface::koBICEventData;
   }
   static void *newArray_koBICInterfacecLcLkoBICEventData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::koBICInterface::koBICEventData[nElements] : new ::koBICInterface::koBICEventData[nElements];
   }
   // Wrapper around operator delete
   static void delete_koBICInterfacecLcLkoBICEventData(void *p) {
      delete (static_cast<::koBICInterface::koBICEventData*>(p));
   }
   static void deleteArray_koBICInterfacecLcLkoBICEventData(void *p) {
      delete [] (static_cast<::koBICInterface::koBICEventData*>(p));
   }
   static void destruct_koBICInterfacecLcLkoBICEventData(void *p) {
      typedef ::koBICInterface::koBICEventData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::koBICInterface::koBICEventData

namespace ROOT {
   static TClass *vectorlEkoBICInterfacecLcLkoBICTowerDatagR_Dictionary();
   static void vectorlEkoBICInterfacecLcLkoBICTowerDatagR_TClassManip(TClass*);
   static void *new_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p = nullptr);
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(Long_t size, void *p);
   static void delete_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p);
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p);
   static void destruct_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<koBICInterface::koBICTowerData>*)
   {
      vector<koBICInterface::koBICTowerData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<koBICInterface::koBICTowerData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<koBICInterface::koBICTowerData>", -2, "vector", 389,
                  typeid(vector<koBICInterface::koBICTowerData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEkoBICInterfacecLcLkoBICTowerDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<koBICInterface::koBICTowerData>) );
      instance.SetNew(&new_vectorlEkoBICInterfacecLcLkoBICTowerDatagR);
      instance.SetNewArray(&newArray_vectorlEkoBICInterfacecLcLkoBICTowerDatagR);
      instance.SetDelete(&delete_vectorlEkoBICInterfacecLcLkoBICTowerDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEkoBICInterfacecLcLkoBICTowerDatagR);
      instance.SetDestructor(&destruct_vectorlEkoBICInterfacecLcLkoBICTowerDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<koBICInterface::koBICTowerData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<koBICInterface::koBICTowerData>","std::vector<koBICInterface::koBICTowerData, std::allocator<koBICInterface::koBICTowerData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICTowerData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEkoBICInterfacecLcLkoBICTowerDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICTowerData>*>(nullptr))->GetClass();
      vectorlEkoBICInterfacecLcLkoBICTowerDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEkoBICInterfacecLcLkoBICTowerDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICTowerData> : new vector<koBICInterface::koBICTowerData>;
   }
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICTowerData>[nElements] : new vector<koBICInterface::koBICTowerData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p) {
      delete (static_cast<vector<koBICInterface::koBICTowerData>*>(p));
   }
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p) {
      delete [] (static_cast<vector<koBICInterface::koBICTowerData>*>(p));
   }
   static void destruct_vectorlEkoBICInterfacecLcLkoBICTowerDatagR(void *p) {
      typedef vector<koBICInterface::koBICTowerData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<koBICInterface::koBICTowerData>

namespace ROOT {
   static TClass *vectorlEkoBICInterfacecLcLkoBICSiPMDatagR_Dictionary();
   static void vectorlEkoBICInterfacecLcLkoBICSiPMDatagR_TClassManip(TClass*);
   static void *new_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p = nullptr);
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(Long_t size, void *p);
   static void delete_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p);
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p);
   static void destruct_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<koBICInterface::koBICSiPMData>*)
   {
      vector<koBICInterface::koBICSiPMData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<koBICInterface::koBICSiPMData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<koBICInterface::koBICSiPMData>", -2, "vector", 389,
                  typeid(vector<koBICInterface::koBICSiPMData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEkoBICInterfacecLcLkoBICSiPMDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<koBICInterface::koBICSiPMData>) );
      instance.SetNew(&new_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR);
      instance.SetNewArray(&newArray_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR);
      instance.SetDelete(&delete_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR);
      instance.SetDestructor(&destruct_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<koBICInterface::koBICSiPMData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<koBICInterface::koBICSiPMData>","std::vector<koBICInterface::koBICSiPMData, std::allocator<koBICInterface::koBICSiPMData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICSiPMData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEkoBICInterfacecLcLkoBICSiPMDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICSiPMData>*>(nullptr))->GetClass();
      vectorlEkoBICInterfacecLcLkoBICSiPMDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEkoBICInterfacecLcLkoBICSiPMDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICSiPMData> : new vector<koBICInterface::koBICSiPMData>;
   }
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICSiPMData>[nElements] : new vector<koBICInterface::koBICSiPMData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p) {
      delete (static_cast<vector<koBICInterface::koBICSiPMData>*>(p));
   }
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p) {
      delete [] (static_cast<vector<koBICInterface::koBICSiPMData>*>(p));
   }
   static void destruct_vectorlEkoBICInterfacecLcLkoBICSiPMDatagR(void *p) {
      typedef vector<koBICInterface::koBICSiPMData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<koBICInterface::koBICSiPMData>

namespace ROOT {
   static TClass *vectorlEkoBICInterfacecLcLkoBICLeakageDatagR_Dictionary();
   static void vectorlEkoBICInterfacecLcLkoBICLeakageDatagR_TClassManip(TClass*);
   static void *new_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p = nullptr);
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(Long_t size, void *p);
   static void delete_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p);
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p);
   static void destruct_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<koBICInterface::koBICLeakageData>*)
   {
      vector<koBICInterface::koBICLeakageData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<koBICInterface::koBICLeakageData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<koBICInterface::koBICLeakageData>", -2, "vector", 389,
                  typeid(vector<koBICInterface::koBICLeakageData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEkoBICInterfacecLcLkoBICLeakageDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<koBICInterface::koBICLeakageData>) );
      instance.SetNew(&new_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR);
      instance.SetNewArray(&newArray_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR);
      instance.SetDelete(&delete_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR);
      instance.SetDestructor(&destruct_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<koBICInterface::koBICLeakageData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<koBICInterface::koBICLeakageData>","std::vector<koBICInterface::koBICLeakageData, std::allocator<koBICInterface::koBICLeakageData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICLeakageData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEkoBICInterfacecLcLkoBICLeakageDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICLeakageData>*>(nullptr))->GetClass();
      vectorlEkoBICInterfacecLcLkoBICLeakageDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEkoBICInterfacecLcLkoBICLeakageDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICLeakageData> : new vector<koBICInterface::koBICLeakageData>;
   }
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICLeakageData>[nElements] : new vector<koBICInterface::koBICLeakageData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p) {
      delete (static_cast<vector<koBICInterface::koBICLeakageData>*>(p));
   }
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p) {
      delete [] (static_cast<vector<koBICInterface::koBICLeakageData>*>(p));
   }
   static void destruct_vectorlEkoBICInterfacecLcLkoBICLeakageDatagR(void *p) {
      typedef vector<koBICInterface::koBICLeakageData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<koBICInterface::koBICLeakageData>

namespace ROOT {
   static TClass *vectorlEkoBICInterfacecLcLkoBICGenDatagR_Dictionary();
   static void vectorlEkoBICInterfacecLcLkoBICGenDatagR_TClassManip(TClass*);
   static void *new_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p = nullptr);
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICGenDatagR(Long_t size, void *p);
   static void delete_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p);
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p);
   static void destruct_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<koBICInterface::koBICGenData>*)
   {
      vector<koBICInterface::koBICGenData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<koBICInterface::koBICGenData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<koBICInterface::koBICGenData>", -2, "vector", 389,
                  typeid(vector<koBICInterface::koBICGenData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEkoBICInterfacecLcLkoBICGenDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<koBICInterface::koBICGenData>) );
      instance.SetNew(&new_vectorlEkoBICInterfacecLcLkoBICGenDatagR);
      instance.SetNewArray(&newArray_vectorlEkoBICInterfacecLcLkoBICGenDatagR);
      instance.SetDelete(&delete_vectorlEkoBICInterfacecLcLkoBICGenDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEkoBICInterfacecLcLkoBICGenDatagR);
      instance.SetDestructor(&destruct_vectorlEkoBICInterfacecLcLkoBICGenDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<koBICInterface::koBICGenData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<koBICInterface::koBICGenData>","std::vector<koBICInterface::koBICGenData, std::allocator<koBICInterface::koBICGenData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICGenData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEkoBICInterfacecLcLkoBICGenDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICGenData>*>(nullptr))->GetClass();
      vectorlEkoBICInterfacecLcLkoBICGenDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEkoBICInterfacecLcLkoBICGenDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICGenData> : new vector<koBICInterface::koBICGenData>;
   }
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICGenDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICGenData>[nElements] : new vector<koBICInterface::koBICGenData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p) {
      delete (static_cast<vector<koBICInterface::koBICGenData>*>(p));
   }
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p) {
      delete [] (static_cast<vector<koBICInterface::koBICGenData>*>(p));
   }
   static void destruct_vectorlEkoBICInterfacecLcLkoBICGenDatagR(void *p) {
      typedef vector<koBICInterface::koBICGenData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<koBICInterface::koBICGenData>

namespace ROOT {
   static TClass *vectorlEkoBICInterfacecLcLkoBICEdepDatagR_Dictionary();
   static void vectorlEkoBICInterfacecLcLkoBICEdepDatagR_TClassManip(TClass*);
   static void *new_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p = nullptr);
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(Long_t size, void *p);
   static void delete_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p);
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p);
   static void destruct_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<koBICInterface::koBICEdepData>*)
   {
      vector<koBICInterface::koBICEdepData> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<koBICInterface::koBICEdepData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<koBICInterface::koBICEdepData>", -2, "vector", 389,
                  typeid(vector<koBICInterface::koBICEdepData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEkoBICInterfacecLcLkoBICEdepDatagR_Dictionary, isa_proxy, 0,
                  sizeof(vector<koBICInterface::koBICEdepData>) );
      instance.SetNew(&new_vectorlEkoBICInterfacecLcLkoBICEdepDatagR);
      instance.SetNewArray(&newArray_vectorlEkoBICInterfacecLcLkoBICEdepDatagR);
      instance.SetDelete(&delete_vectorlEkoBICInterfacecLcLkoBICEdepDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlEkoBICInterfacecLcLkoBICEdepDatagR);
      instance.SetDestructor(&destruct_vectorlEkoBICInterfacecLcLkoBICEdepDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<koBICInterface::koBICEdepData> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<koBICInterface::koBICEdepData>","std::vector<koBICInterface::koBICEdepData, std::allocator<koBICInterface::koBICEdepData> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICEdepData>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEkoBICInterfacecLcLkoBICEdepDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<koBICInterface::koBICEdepData>*>(nullptr))->GetClass();
      vectorlEkoBICInterfacecLcLkoBICEdepDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEkoBICInterfacecLcLkoBICEdepDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICEdepData> : new vector<koBICInterface::koBICEdepData>;
   }
   static void *newArray_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<koBICInterface::koBICEdepData>[nElements] : new vector<koBICInterface::koBICEdepData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p) {
      delete (static_cast<vector<koBICInterface::koBICEdepData>*>(p));
   }
   static void deleteArray_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p) {
      delete [] (static_cast<vector<koBICInterface::koBICEdepData>*>(p));
   }
   static void destruct_vectorlEkoBICInterfacecLcLkoBICEdepDatagR(void *p) {
      typedef vector<koBICInterface::koBICEdepData> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<koBICInterface::koBICEdepData>

namespace ROOT {
   static TClass *maplEpairlEfloatcOfloatgRcOintgR_Dictionary();
   static void maplEpairlEfloatcOfloatgRcOintgR_TClassManip(TClass*);
   static void *new_maplEpairlEfloatcOfloatgRcOintgR(void *p = nullptr);
   static void *newArray_maplEpairlEfloatcOfloatgRcOintgR(Long_t size, void *p);
   static void delete_maplEpairlEfloatcOfloatgRcOintgR(void *p);
   static void deleteArray_maplEpairlEfloatcOfloatgRcOintgR(void *p);
   static void destruct_maplEpairlEfloatcOfloatgRcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<pair<float,float>,int>*)
   {
      map<pair<float,float>,int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<pair<float,float>,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<pair<float,float>,int>", -2, "map", 100,
                  typeid(map<pair<float,float>,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEpairlEfloatcOfloatgRcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<pair<float,float>,int>) );
      instance.SetNew(&new_maplEpairlEfloatcOfloatgRcOintgR);
      instance.SetNewArray(&newArray_maplEpairlEfloatcOfloatgRcOintgR);
      instance.SetDelete(&delete_maplEpairlEfloatcOfloatgRcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEpairlEfloatcOfloatgRcOintgR);
      instance.SetDestructor(&destruct_maplEpairlEfloatcOfloatgRcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<pair<float,float>,int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<pair<float,float>,int>","std::map<std::pair<float, float>, int, std::less<std::pair<float, float> >, std::allocator<std::pair<std::pair<float, float> const, int> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<pair<float,float>,int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEpairlEfloatcOfloatgRcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<pair<float,float>,int>*>(nullptr))->GetClass();
      maplEpairlEfloatcOfloatgRcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEpairlEfloatcOfloatgRcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEpairlEfloatcOfloatgRcOintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<pair<float,float>,int> : new map<pair<float,float>,int>;
   }
   static void *newArray_maplEpairlEfloatcOfloatgRcOintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<pair<float,float>,int>[nElements] : new map<pair<float,float>,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEpairlEfloatcOfloatgRcOintgR(void *p) {
      delete (static_cast<map<pair<float,float>,int>*>(p));
   }
   static void deleteArray_maplEpairlEfloatcOfloatgRcOintgR(void *p) {
      delete [] (static_cast<map<pair<float,float>,int>*>(p));
   }
   static void destruct_maplEpairlEfloatcOfloatgRcOintgR(void *p) {
      typedef map<pair<float,float>,int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<pair<float,float>,int>

namespace {
  void TriggerDictionaryInitialization_librootIO_Impl() {
    static const char* headers[] = {
"/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include/RootInterface.h",
"/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include/koBICInterface.h",
nullptr
    };
    static const char* includePaths[] = {
"/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include",
"/u/user/changhui/software/root/include",
"/u/user/changhui/koBIC2026/TBCERN2026b/rootIO",
"/u/user/changhui/software/root/include/",
"/u/user/changhui/koBIC2026/TBCERN2026b/build/rootIO/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "librootIO dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <typename ..._Elements> class __attribute__((annotate("$clingAutoload$type_traits")))  __attribute__((annotate("$clingAutoload$string")))  tuple;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "librootIO dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include/RootInterface.h"
#include "/u/user/changhui/koBIC2026/TBCERN2026b/rootIO/include/koBICInterface.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"koBICInterface", payloadCode, "@",
"koBICInterface::koBICEdepData", payloadCode, "@",
"koBICInterface::koBICEventData", payloadCode, "@",
"koBICInterface::koBICGenData", payloadCode, "@",
"koBICInterface::koBICLeakageData", payloadCode, "@",
"koBICInterface::koBICSiPMData", payloadCode, "@",
"koBICInterface::koBICTowerData", payloadCode, "@",
"koBICInterface::threeVector", payloadCode, "@",
"tuple<float,float,float>", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("librootIO",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_librootIO_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_librootIO_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_librootIO() {
  TriggerDictionaryInitialization_librootIO_Impl();
}
