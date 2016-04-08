/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Null.hpp
 * Author: richard
 *
 * Created on April 6, 2016, 6:09 PM
 */

#ifndef NULL_HPP
#define NULL_HPP

#include <bpmodule/modulebase/SystemFragmenter.hpp>

/** \brief A class that generates a null fragmentation scheme
 *
 *  In order to use methods that assume fragments exist, such as MIM,
 *  with the entire system we need a way of telling MIM not to fragment.
 *  We do this with this class, which just returns a SystemMap that is
 *  the input system.
 *  
 */
class NullFragmenter: public bpmodule::modulebase::SystemFragmenter{
    public:
    ///Import the constructor from the base class
    using bpmodule::modulebase::SystemFragmenter::SystemFragmenter;

    ///Returns a SystemMap with the input system
    virtual bpmodule::system::SystemMap Fragmentize_
           (const bpmodule::system::System & mol);
};


#endif /* NULL_HPP */

