/**
 * @file ParticleRRHO.cpp
 *
 * @brief Implementation of the ParticleRRHO class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


#include "ParticleRRHO.h"
#include "Utilities.h"
#include "Constants.h"

#include <cassert>
#include <iostream>
using std::cout; // debug
using std::endl; // debug

using namespace Mutation::Utilities;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

ParticleRRHO::ParticleRRHO(const IO::XmlElement& xml_element)
    : m_hform(0), m_steric(0), m_linearity(0), m_rotational_t(0)
{
    // Load information stored in child elements
    IO::XmlElement::const_iterator iter = xml_element.begin();
    
    for ( ; iter != xml_element.end(); ++iter) {
        if (iter->tag() == "formation_enthalpy") {
            m_hform = atof(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "steric_factor") {
            m_steric = atoi(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "linear") {
            std::string yesno = String::trim(iter->text());
            if (yesno == "yes")
                m_linearity = 2;
            else if (yesno == "no")
                m_linearity = 3;
            else {
                iter->parseError(
                    "Values for linear can only be \"yes\" or \"no\"!");
            }
        }
        else if (iter->tag() == "rotational_temperature") {
            m_rotational_t = atof(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "vibrational_temperatures") {
            std::vector<std::string> tokens;
            String::tokenize(iter->text(), tokens, " ,\t\n\r\f\v");
            
            std::vector<std::string>::const_iterator t_iter = tokens.begin();
            for ( ; t_iter != tokens.end(); ++t_iter)
                m_vibrational_energies.push_back(atof(t_iter->c_str()));
        }
        else if (iter->tag() == "electronic_levels") {
            IO::XmlElement::const_iterator level_iter = iter->begin();
            
            int    degeneracy;
            double temperature;
            
            for ( ; level_iter != iter->end(); ++level_iter) {
                if (level_iter->tag() == "level") {
                    level_iter->getAttribute("degeneracy",  degeneracy);
                    level_iter->getAttribute("energy", temperature);
                    
                    // convert from 1/cm to K
                    m_electronic_energies.push_back(
                        std::make_pair(degeneracy, temperature * 1.4387));
                }
            }
        }
    }
}

//==============================================================================

ParticleRRHO::ParticleRRHO(const ParticleRRHO& rrho, const size_t level)
    : m_hform(rrho.m_hform + RU*rrho.electronicEnergy(level).second),
      m_steric(rrho.m_steric),
      m_linearity(rrho.m_linearity),
      m_rotational_t(rrho.m_rotational_t),
      m_electronic_energies(1, std::make_pair(rrho.electronicEnergy(level).first, 0.0)),
      m_vibrational_energies(rrho.m_vibrational_energies)
{
    // Make sure the level used is actually present in the given RRHO parameters
    assert(level < rrho.nElectronicLevels());
}

//==============================================================================

ParticleRRHO::ParticleRRHO(const IO::XmlElement& xml_element,
                           const std::vector<size_t>& indices)
    : m_hform(0), m_steric(0), m_linearity(0), m_rotational_t(0)
{
    /// Load basic data -------------------------------------------------------
    IO::XmlElement::const_iterator iter = xml_element.begin();
    
    for ( ; iter != xml_element.end(); ++iter) {
        if (iter->tag() == "formation_enthalpy") {
            m_hform = atof(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "steric_factor") {
            m_steric = atoi(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "linear") {
            std::string yesno = String::trim(iter->text());
            if (yesno == "yes")
                m_linearity = 2;
            else if (yesno == "no")
                m_linearity = 3;
            else {
                iter->parseError(
                    "Values for linear can only be \"yes\" or \"no\"!");
            }
        }
        else if (iter->tag() == "rotational_temperature") {
            m_rotational_t = atof(String::trim(iter->text()).c_str());
        }
        else if (iter->tag() == "vibrational_temperatures") {
            std::vector<std::string> tokens;
            String::tokenize(iter->text(), tokens, " ,\t\n\r\f\v");
            
            std::vector<std::string>::const_iterator t_iter = tokens.begin();
            for ( ; t_iter != tokens.end(); ++t_iter)
                m_vibrational_energies.push_back(atof(t_iter->c_str()));
        }
    }
    
    // Load data representing a single energy level ---------------------------
    
    m_electronic_energies.push_back(std::make_pair(1, 0.0));
    
    std::vector<std::string> modeTag = {"electronic_levels",
        "vibrational_levels", "rotational_levels"};
    IO::XmlElement::const_iterator level_it = xml_element.findTag(modeTag[0]);
    iter = level_it;
    
    for (int i = 0; i < indices.size(); ++i) {
        
        if (i > 0)
            iter = level_it->findTag(modeTag[i]);
        
        if (iter == level_it->end()) {
            iter->parseError(
                "Requested excited level data not found !");
        }
            
        level_it = iter->begin();
        size_t level = 0;
        for ( ; level_it != iter->end(); ++level_it) {
            if (level_it->tag() == "level") {
                
                if (level == indices[i]) {
                    
                    int    degeneracy;
                    double temperature;
                    level_it->getAttribute("degeneracy",  degeneracy);
                    level_it->getAttribute("energy", temperature);
                    // convert from 1/cm to J/mol
                    m_hform += RU*temperature*1.4387;
                    m_electronic_energies[0].first *= degeneracy;
                    
                    std::cout << "  degeneracy = " << degeneracy; // debug
                    std::cout << " energy = " << temperature << endl; // debug
                    
                    if (i < 1) {
                        IO::XmlElement::const_iterator t_it = 
                            level_it->findTag("vibrational_temperature");
                        if (t_it != level_it->end()) {
                            m_vibrational_energies.clear();
                            m_vibrational_energies.push_back(
                                atof(String::trim(t_it->text()).c_str()) );
                        }
                    } else
                        m_vibrational_energies.clear();
                        
                    if (i < 2){
                        IO::XmlElement::const_iterator t_it = 
                            level_it->findTag("rotational_temperature");
                        if (t_it != level_it->end())
                            m_rotational_t = atof(String::trim(t_it->text()).c_str());
                    } else {
                        m_rotational_t = 0.0;
                        m_linearity = 0;
                    }
                    
                    break;
                }
                
                level++;
            }
        }
            
    }
    
    for (size_t i = 0; i < m_vibrational_energies.size(); ++i)
        std::cout << "Tvib = " << m_vibrational_energies[i] << endl; // debug
    std::cout << "Trot = " << m_rotational_t << endl; // debug
    
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation

