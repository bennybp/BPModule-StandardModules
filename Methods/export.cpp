/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <pybind11/pybind11.h>
#include "Methods/MethodHelpers/MethodHelpers.hpp"

PYBIND11_PLUGIN(Methods){
pybind11::module mtop("Methods", "Python bindings for the reference implementations");

//pybind11::module m=mtop.def_submodule("pulsarmethods");
mtop.def("run_series_of_methods",pulsarmethods::RunSeriesOfMethods);
mtop.def("fill_deriv",pulsarmethods::FillDeriv);
 
return mtop.ptr();
}
