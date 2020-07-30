/**
 * @file OmegaSEv.cpp
 *
 * @brief Implementation of OmegaSEv.
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

#include "Mixture.h"
#include "TransferModel.h"

#include <cmath>
#include <vector>

#include <iostream>

namespace Mutation {
    namespace Transfer {
      
/**
 * Gives the average vibrational energy gain for a single spontaneous
 * emission process, as a function of vibrational temperature,
 * fitted to a rational function. Units : J/mol
 */
class FittedSEVibSourceTerm
{
public:

    FittedSEVibSourceTerm(const Mutation::Utilities::IO::XmlElement& node) {
        node.getAttribute("a0", m_a0, 0.0);
        node.getAttribute("a1", m_a1, 0.0);
        node.getAttribute("a2", m_a2, 0.0);
        node.getAttribute("a3", m_a3, 0.0);
        node.getAttribute("b0", m_b0, 0.0);
        node.getAttribute("b1", m_b1, 0.0);
        node.getAttribute("b2", m_b2, 0.0);
        
        // debug
        std::cout << "a0 " << m_a0 << " a1 " << m_a1 << " a2 " << m_a2 
            << " a3 " << m_a3 << "b0" << m_b0 << " b1 " << m_b1
            << " b2 " << m_b2 << std::endl;
        
        double Tv, sv;
        for (int i=0 ; i < 100 ; i++) {
            Tv = i*100000./100.;
            sv = (m_a0 + (m_a1 + (m_a2 + m_a3*Tv)*Tv)*Tv) /
                 (m_b0 + (m_b1 + m_b2*Tv)*Tv;
            std::cout << Tv << " " << sv << std::endl;
        }
    }
       
    double rate(const double Tv) {
        return ( (m_a0 + (m_a1 + (m_a2 + m_a3*Tv)*Tv)*Tv) /
            (m_b0 + (m_b1 + m_b2*Tv)*Tv) );
    }

private:

    double m_a0, m_a1, m_a2, m_a3;
    double m_b0, m_b1, m_b2;
};

//=============================================================================
        
class OmegaSEv : public TransferModel
{

public:
    OmegaSEv(Mutation::Mixture& mix)
        : TransferModel(mix)
    {
        m_ns = m_mixture.nSpecies();
        m_nr = m_mixture.nReactions();
        mp_hv = new double [m_ns];
        mp_rate = new double [m_nr];
        mp_delta = new double [m_nr];
        
        for(int i = 0; i < m_nr; ++i) {
            // Store any SE reaction involving molecules
            if ( (mix.reactions()[i].type() == Kinetics::BND_BND_EMISSION)
                && (mix.species(mix.reactions()[i].reactants()[0]).type() 
                    == Mutation::Thermodynamics::MOLECULE) )
                m_rId.push_back(i);
        }
        
        if (m_rId.empty()) return;
        
        // If there is any SE reaction in the mixture, look for data
        std::string filename = databaseFileName("SEvibSource.xml", "transfer");
        Mutation::Utilities::IO::XmlDocument doc(filename);
        
        Mutation::Utilities::IO::XmlElement::const_iterator iter;
        
        // Loop over all SE reactions and load the data for the fitted source term
        for(int i = 0; i < m_rId.size(); ++i) {
            
            iter = doc.root().findTagWithAttribute("transition",
                "formula", mix.reactions()[ m_rId[i] ].formula())
                
            if (iter == doc.root().end())
                doc.root().parseError("Could not find requested transition.");
                
            m_rSourceTerm.push_back(FittedSEVibSourceTerm(*iter));
        }
    };

    ~OmegaSEv() {
        delete [] mp_hv;
        delete [] mp_rate;
        delete [] mp_delta;
    };

    /**
      * Computes the vibrational energy removed by spontaneous 
      * emission processes in \f$ [J/(m^3\cdot s)] \f$.
      */
    double source()
    {
        // Get vibrational energies
        m_mixture.speciesHOverRT(NULL, NULL, NULL, mp_hv, NULL, NULL);
        
        // Get reaction non-preferential vibrational energy change
        std::fill(mp_delta, mp_delta+m_nr, 0.0);
        m_mixture.getReactionDelta(mp_hv, mp_delta);

        // Get molar rates of progress
        m_mixture.netRatesOfProgress(mp_rate);

        double src = 0.0;
        double numdens;
        int j;
        for (int i = 0; i < m_rId.size(); ++i) {
            j = m_rId[i];
            
            // Compute contribution from OmegaCV model and remove it
            src -= mp_delta[j]*mp_rate[j]*RU*m_mixture.T();
            
            // Add vibrational energy gain for this process
            src += m_rSourceTerm[i].rate(m_mixture.Tv())*mp_rate[j];
        }

        return src;

    }

private:
    int m_ns;
    int m_nr;
    std::vector<int> m_rId;
    std::vector<FittedSEVibSourceTerm> m_rSourceTerm;
    double* mp_hv;
    double* mp_rate;
    double* mp_delta;
};
  
// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaSEv, TransferModel> omegaSEv("OmegaSEv");

    } // namespace Transfer
} // namespace Mutation
