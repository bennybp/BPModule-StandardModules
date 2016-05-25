/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "SystemFragmenters/Null.hpp"

pulsar::system::SystemMap NullFragmenter::Fragmentize_(
    const pulsar::system::System & mol){return {{"SYSTEM",mol}};}

