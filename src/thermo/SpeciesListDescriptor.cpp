/**
 * @file SpeciesListDescriptor.cpp
 *
 * @brief Implementation of the SpeciesListDescriptor class.
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

#include "SpeciesListDescriptor.h"
#include "Species.h"
#include "Utilities.h"

using std::cout;
using std::endl;

#include <iostream>
#include <iterator>
using namespace std;

namespace Mutation {
    namespace Thermodynamics {

//==============================================================================

SpeciesListDescriptor::SpeciesListDescriptor(std::string descriptor)
    : m_gases(false), m_solids(false), m_liquids(false)
{
    // First step is to look for implicit species defined by "{ rules }"
    size_t is = descriptor.find_first_of('{');
    if (is != std::string::npos) {
        size_t es = descriptor.find_first_of('}', is);
        std::string rule = descriptor.substr(is+1,es-is-1);
        descriptor.erase(is, es-is+1);

        // Split the rules up
        std::vector<std::string> rules;
        Utilities::String::tokenize(rule, rules, "with", false);
        
        // First check phase rule
        std::vector<std::string> phases;
        Utilities::String::tokenize(rules[0], phases, ", ");
        
        for (size_t i = 0; i < phases.size(); ++i) {
            if (phases[i] == "gases")
                m_gases = true;
            else if (phases[i] == "liquids")
                m_liquids = true;
            else if (phases[i] == "solids")
                m_solids = true;
            else if (phases[i] == "condensed") {
                m_solids = true;
                m_liquids = true;
            } else if (phases[i] == "all") {
                m_gases = true;
                m_solids = true;
                m_liquids = true;
            } else {
                throw InvalidInputError("species descriptor", descriptor)
                    << "Unknown phase keyword in implicit species rule.  "
                    << "Possible phase descriptors are 'gases', 'liquids', "
                    << "'solids', 'condensed', and 'all'.";
            }
        }
        
        // Second the elements
        std::vector<std::string> elements;
        Utilities::String::tokenize(rules[1], elements, ", ");
        m_element_names.insert(elements.begin(), elements.end());
    }
    
    // Separate out the species names
    separateSpeciesNames(descriptor);
    
    
    
    // Check on any species whose excited states should be expanded
    std::set<std::string> gsNames;
    std::vector<std::string>::iterator it = m_species_names.begin();
    
    while (it != m_species_names.end()) {
        
        std::string gsName;
        std::vector<size_t> indices;
        int expd = 0;
        
        size_t spos = it->find('(');
        
        if (spos != std::string::npos) {
            
            size_t epos = it->find(')', spos);
            
            std::vector<std::string> tokens;
            Utilities::String::tokenize(it->substr(spos+1, epos-spos-1), tokens, ",");
            
            for (int i = 0; i < tokens.size(); ++i) {
                if (tokens[i] == "*")
                    expd++;
                else if (expd == 0)
                    indices.push_back(atoi(tokens[i].c_str()));
                else
                    throw InvalidInputError("species descriptor", descriptor)
                        << "Star token can only be followed by another star token. \n"
                        << "    " << *it << " <--";
            }
            
            gsName = it->substr(0, spos);
            
        } else {
            gsName = *it;
        }
        
        // insert species into the map with its state expansion indicator
        energyLevel key(gsName, indices);
        std::pair<levelMap::iterator, bool> ret;
        ret = m_expand_states.insert( std::pair<energyLevel, int>(key, expd) );
        if (!ret.second)
            ret.first->second = std::max(ret.first->second, expd);
        
        *it = gsName;
        if (gsNames.count(gsName) > 0)
            (it = m_species_names.erase(it))--;
        else
            gsNames.insert(gsName);
        
        it++;
    }
    
    
    // look for duplicate or conflicting definitions of excited species
    levelMap::const_iterator rit1 = m_expand_states.end();
    levelMap::const_iterator rit2;
    std::string gsName;
    
    cout << "map ------ " << endl; // Debug ------
    for (rit1 = m_expand_states.begin() ; rit1 != m_expand_states.end(); ++rit1) {
        cout << rit1->first.groundStateName() << " indices:";
        for (int i = 0; i < rit1->first.indices().size(); ++i)
            cout << " " << rit1->first.indices()[i];
        cout << "  exp = " << rit1->second << endl;
        cout << endl;
    }
    rit1 = m_expand_states.end(); // ------ Debug
    
    rit1--;
    while (rit1 != m_expand_states.begin()) {
        
        (rit2 = rit1)--;
        gsName = rit1->first.groundStateName();
        
        while (rit2->first.groundStateName() == gsName) {
            
            size_t nsub = rit2->first.indices().size();
            std::vector<size_t> subind(rit1->first.indices().begin(),
                                        rit1->first.indices().begin()+nsub);
            
            if (subind == rit2->first.indices()) {
                if ( rit2->second == (rit1->first.indices().size()-nsub) && rit1->second == 0 ) {
                    // redundant definition
                    rit1 = m_expand_states.erase(rit1);
                    break;
                } else
                    // conflicting definition
                    throw InvalidInputError("species descriptor", descriptor)
                        << "Conflicting definitions of excited states: "
                        << rit1->first.name().append(rit1->second, '*')
                        << " <--> " << rit2->first.name().append(rit2->second, '*');
            }
            
            if (rit2 == m_expand_states.begin())
                break;
            rit2--;
        }
        
        rit1--;
    }
    
    cout << "final map ----------- " << endl; // Debug ------
    for (rit1 = m_expand_states.begin() ; rit1 != m_expand_states.end(); ++rit1) {
        cout << rit1->first.groundStateName() << " indices:";
        for (int i = 0; i < rit1->first.indices().size(); ++i)
            cout << " " << rit1->first.indices()[i];
        cout << "  exp = " << rit1->second << endl;
        cout << endl; // ------ Debug
    }
    
}

//==============================================================================

void SpeciesListDescriptor::separateSpeciesNames(std::string descriptor)
{
    // First trim the descriptor string
    descriptor = Utilities::String::trim(descriptor);

    // State-machine with 2 states (either in quotes or not)
    bool in_quotes = (descriptor[0] == '\"');
    std::string name = "";

    for (int i = (in_quotes ? 1 : 0); i < descriptor.length(); ++i) {
        char c = descriptor[i];

        if (in_quotes) {
            // Add everything to the name until we escape out of the quotes
            switch(c) {
            case '\"':
                // Escape and add name to the list
                in_quotes = false;
                if (name.length() > 0) {
                    m_species_names.push_back(name);
                    name = "";
                }
                break;
            default:
                // Otherwise just add the character to the name
                name += c;
            }
        } else {
            switch(c) {
            // All white-space should be ignored
            case ' ': case '\n': case '\r': case '\t': case '\f': case '\v':
                if (name.length() > 0) {
                    m_species_names.push_back(name);
                    name = "";
                }
                break;
            // Cannot include quotation mark in a name
            case '\"':
                if (name.length() > 0) {
                    throw InvalidInputError("species name", name)
                        << "Cannot include quotation mark in species name.\n"
                        << "    " << descriptor.substr(0, i+1) << " <--";
                }
                in_quotes = true;
                break;
            default:
                name += c;
            }
        }
    }

    // Push back the last name
    if (name.length() > 0)
        m_species_names.push_back(name);
}

//==============================================================================

bool SpeciesListDescriptor::matches(const Species& species) const
{
    // Check if this species is present in the explicit list (includes excited
    // state species which were implicitly defined using the '*' character)
    for (int i = 0; i < m_species_names.size(); ++i) {
        
        const std::string& name = m_species_names[i];
        
        if (species.groundStateName() == name) {
            
            levelMap::const_iterator it;
            std::vector<size_t> indices;
            
            indices = energyLevel(species).indices();       
            
            for (int j = 0; j < species.levelType(); ++j) {
                
                energyLevel key(species.groundStateName(), indices);
                it = m_expand_states.find(key);
                
                if ( (it != m_expand_states.end()) && (it->second == j) )
                    return true;
                
                indices.pop_back();
            }
            
            energyLevel key(species.groundStateName(), indices);
            it = m_expand_states.find(key);
            
            if ( (it != m_expand_states.end()) && 
                (it->second == species.levelType()) )
                return true;
            else
                return false;
            
        }
    }
    
    // do not apply implicit rules to excited states species
    if (species.levelType() != NONE)
        return false;
    
    // Check if the species is described by an implicit rule
    // Check phase
    if (species.phase() == GAS    && !m_gases)   return false;
    if (species.phase() == SOLID  && !m_solids)  return false;
    if (species.phase() == LIQUID && !m_liquids) return false;
    
    // Check elements
    Species::StoichList::const_iterator iter =
        species.stoichiometry().begin();
    for ( ; iter != species.stoichiometry().end(); ++iter)
        if (m_element_names.count(iter->first) == 0) return false;
    
    return true;
}

//==============================================================================

/// Predicate returns true if species name equals given name.
struct NameEquals {
    NameEquals(const std::string& str) : name(str) { }
    bool operator()(const Species& s) { return s.name() == name; }
    std::string name;
};

/// Predicate returns true if species ground state name equals given name.
struct GroundStateNameEquals {
    GroundStateNameEquals(const std::string& str) : name(str) { }
    bool operator()(const Species& s) { return s.groundStateName() == name; }
    std::string name;
};

/// Predicate returns true if species type equals given type.
struct TypeEquals {
    TypeEquals(ParticleType t) : type(t) { }
    bool operator()(const Species& s) { return s.type() == type; }
    ParticleType type;
};

/// Predicate returns true if species is condensed.
struct IsCondensed {
    bool operator()(const Species& s) { return s.phase() != GAS; }
};

void SpeciesListDescriptor::order(
    std::list<Species>& input, std::vector<Species>& output,
    std::vector<std::string>& missing) const
{
    std::list<Species> ordered_list;
    
    int index = 0;
    std::list<Species>::iterator it;
    levelMap::const_iterator it1;
    
    // First order all of the species that are explicitly listed
    for (std::size_t i = 0; i < m_species_names.size(); ++i) {
        const std::string& name = m_species_names[i];
        
        // Don't expand name
        it1 = m_expand_states.find(energyLevel(name, std::vector<size_t>()));
        if ( it1 != m_expand_states.end() && it1->second == 0 ) {
            it = std::find_if(input.begin(), input.end(), NameEquals(name));
            
            if (it == input.end())
                missing.push_back(name);
            else {
                ordered_list.push_back(*it);
                input.erase(it);
            }
        // Expand name
        } else {
            it = std::find_if(
                input.begin(), input.end(), GroundStateNameEquals(name));
            
            std::list<Species>::iterator lowest_level = it;
            while (lowest_level != input.end()) {
                it++;
                it = std::find_if(it, input.end(), GroundStateNameEquals(name));
                while (it != input.end()) {
                    if ( energyLevel(*it) < energyLevel(*lowest_level) )
                        lowest_level = it;
                    it++;
                    it = std::find_if(
                        it, input.end(), GroundStateNameEquals(name));
                }
            
                ordered_list.push_back(*lowest_level);
                input.erase(lowest_level);
                
                it = std::find_if(
                    input.begin(), input.end(), GroundStateNameEquals(name));
                lowest_level = it;
            }
        }
    }
    
    // Check if any explicitly listed excited state is missing
    for (it1 = m_expand_states.begin() ; it1 != m_expand_states.end(); ++it1) {
        if (it1->second == 0) {
            std::string name = it1->first.name(); 
            it = std::find_if(ordered_list.begin(), ordered_list.end(), NameEquals(name));
            if (it == ordered_list.end())
                missing.push_back(name);
        }
    }
    
    // Return early if we are missing some species as this is an error
    if (missing.size() > 0)
        return;
    
    // Now all that remains are the species that were implicitly defined
    ordered_list.insert(ordered_list.end(), input.begin(), input.end());
    input.clear();
    
    // Place the electron first
    it = find_if(
         ordered_list.begin(), ordered_list.end(), TypeEquals(ELECTRON));

    if (it != ordered_list.begin() && it != ordered_list.end())
        ordered_list.splice(ordered_list.begin(), ordered_list, it);

    // Move all condensed species to the end of the list
    int ncond = count_if(
        ordered_list.begin(), ordered_list.end(), IsCondensed());

    while (ncond-- > 0) {
        it = std::find_if(
            ordered_list.begin(), ordered_list.end(), IsCondensed());
        ordered_list.splice(ordered_list.end(), ordered_list, it);
    }

    // Finally copy the list to the output vector
    output.assign(ordered_list.begin(), ordered_list.end());
}

//==============================================================================

    } // namespace Thermodynamics
} // namespace Mutation
