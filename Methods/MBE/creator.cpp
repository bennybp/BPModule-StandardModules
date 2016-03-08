#include "MBE.hpp"

#define RegisterSuperModule(SMType,SMName)\
   using bpmodule::modulemanager::ModuleCreationFuncs;\
   extern "C" {\
      ModuleCreationFuncs InsertSupermodule(void){\
         ModuleCreationFuncs cf;\
         cf.AddCppCreator<SMType>(SMName);\
         return cf;\
      }\
    }

RegisterSuperModule(LibMBE::MBE,"MBE")

