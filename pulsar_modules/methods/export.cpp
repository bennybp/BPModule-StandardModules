/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <pulsar/util/Pybind11.hpp>
#include "pulsar_modules/methods/method_helpers/MethodHelpers.hpp"

PYBIND11_PLUGIN(pulsar_modules){
pybind11::module mtop("pulsar_modules", "Default Pulsar Implementations");

//pybind11::module m=mtop.def_submodule("pulsarmethods");
mtop.def("run_series_of_methods",RunSeriesOfMethods);
mtop.def("fill_deriv",FillDeriv);
 
return mtop.ptr();
}
