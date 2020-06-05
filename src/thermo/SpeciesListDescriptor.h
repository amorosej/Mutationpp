/**
 * @file SpeciesListDescriptor.h
 *
 * @brief Declaration of the SpeciesListDescriptor class.
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

#ifndef THERMO_SPECIES_DESCRIPTOR_H
#define THERMO_SPECIES_DESCRIPTOR_H

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include "Species.h"
#include <sstream>

namespace Mutation {
    namespace Thermodynamics {

class Species;

/**
 * This is used upon initialization of the thermodynamic database
 * to parse and sort species names representing individual energy levels.
 */
class energyLevel
{
public:

    // Constructors.
    energyLevel(const std::string& name, const std::vector<size_t>& indices)
        : m_groundStateName(name), m_indices(indices) { }
    
    energyLevel(const Species& species)
        : m_groundStateName(species.groundStateName())
    {
        if (species.levelType() >= ELECTRONIC)
            m_indices.push_back(species.level());
        
        if (species.levelType() >= VIBRATIONAL)
            m_indices.push_back(species.vibLevel());
        
        if (species.levelType() >= ROTATIONAL)
            m_indices.push_back(species.rotLevel());
    }
    
    // Returns the full name of this level.
    std::string name() const {
        if (m_indices.size() > 0) {
            std::stringstream ind;
            ind << "(";
            for (size_t i = 0; i < m_indices.size()-1; ++i)
                ind << m_indices[i] << ",";
            ind << m_indices.back() << ")";
            return m_groundStateName+ind.str();
        } else
            return m_groundStateName;
    }
    
    // Returns the name of the ground state.
    const std::string& groundStateName() const {
        return m_groundStateName;
    }
    
    // Returns the indices representing this level.
    const std::vector<size_t>& indices() const { 
        return m_indices; 
    }
    
    bool operator<(const energyLevel& lvl) const
    {
        if (this->groundStateName() != lvl.groundStateName())
           return (this->groundStateName() < lvl.groundStateName());
        
        //if (this->indices().size() != lvl.indices().size())
        //   return (this->indices().size() < lvl.indices().size());
        
        return (this->indices() < lvl.indices());
    }
    
private:
    
    std::string         m_groundStateName;
    std::vector<size_t> m_indices;
};


/**
 * This class is used by thermodynamic databases to decide which species they 
 * should load upon initialization.  The species list could be determined from
 * a simple list of species names, or something more complex such as all gases
 * containing certain elements.
 */
class SpeciesListDescriptor
{
public:
    
    typedef std::map<energyLevel, int> levelMap;
    
    /**
     * Constructor which takes a string representation of the descriptor.
     */
    SpeciesListDescriptor(std::string descriptor);
    
    /**
     * Destructor.
     */
    ~SpeciesListDescriptor() { }
    
    /**
     * Tests if a species object is described by this descriptor and keeps track
     * of whether or not all explicitly defined species are matched.
     */
    bool matches(const Species& species) const;
    
    /**
     * Orders the species given as input in the output array.  Ensures that 
     * species explicitly listed by the user maintain the same order and that 
     * the electron is always at the beginning if it is present as a species.
     * Condensed phase species are listed at the end of the array.
     */
    void order(
        std::list<Species>& input, std::vector<Species>& output,
        std::vector<std::string>& missing) const;
    
private:

    /// Separates species names in a list (initializes m_species_names).
    void separateSpeciesNames(std::string descriptor);


private:

    /// Explicitly defined species names
    std::vector<std::string> m_species_names;
    
    /// List of allowed elements
    std::set<std::string> m_element_names;
    
    /// List of species who should be expanded into excited electronic states
    //std::set<std::string> m_expand_states;
    levelMap m_expand_states;
    
    /// True if gases are allowed
    bool m_gases;
    
    /// True if solids are allowed
    bool m_solids;
    
    /// True if liquids are allowed
    bool m_liquids;
    
};

    } // namespace Thermodynamics
} // namespace Mutation

#endif

