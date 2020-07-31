/**
 * @file RateLaws.h
 *
 * @brief Declaration of various RateLaw classes.
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

#ifndef RATELAW_H
#define RATELAW_H

#include <vector>
#include <cmath>
#include <cstdlib>

#include "Utilities.h"

namespace Mutation {
    namespace Kinetics {

/**
 * Abstract base class for all rate laws which allows owners such as class 
 * Reaction to store any rate law polymorphically.
 */
class RateLaw
{
public:

    virtual ~RateLaw() { };
    virtual RateLaw* clone() const = 0;
};

/**
 * Arrhenius rate law \f$ k_f(T) = A T^\eta exp(-E_a / (R_u T)) \f$.
 */
class Arrhenius : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    Arrhenius(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    Arrhenius(const Arrhenius& to_copy)
        : m_lnA(to_copy.m_lnA), m_n(to_copy.m_n), m_temp(to_copy.m_temp)
    { }
    
    virtual ~Arrhenius() { };
    
    Arrhenius* clone() const {
        return new Arrhenius(*this);
    }
    
    inline double getLnRate(const double lnT, const double invT) const {
        return (m_lnA + m_n * lnT - m_temp * invT);
    }
    
    inline double derivative(const double k, const double lnT, const double invT) const {
        return (k*invT*(m_n + m_temp*invT));
    }

    double A() const { 
        return std::exp(m_lnA);
    }
    
    double n() const {
        return m_n;
    }
    
    double T() const {
        return m_temp;
    }
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;    
    static std::vector<Mutation::Utilities::Units> sm_eunits;

    double m_lnA;
    double m_n;
    double m_temp;
};


/**
 * Arrhenius-like rate law with a rational pre-exponential term.
 */
class rationalExp : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    rationalExp(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    rationalExp(const rationalExp& to_copy)
        : m_n(to_copy.m_n), m_temp(to_copy.m_temp), 
          m_a0(to_copy.m_a0), m_a1(to_copy.m_a1), m_a2(to_copy.m_a2),
          m_b0(to_copy.m_b0), m_b1(to_copy.m_b1), m_b2(to_copy.m_b2), m_b3(to_copy.m_b3)
    { }
    
    virtual ~rationalExp() { };
    
    rationalExp* clone() const {
        return new rationalExp(*this);
    }
    
    inline double getLnRate(const double lnT, const double invT, const double T, const double sqT) const {
        return (m_n*lnT - m_temp*invT + std::log( (m_a0+m_a1*T+m_a2*sqT)/(m_b0+m_b1*T+m_b2*sqT+m_b3*sqT*T) ));
    }
    
    inline double derivative(const double k, const double invT, const double T, const double sqT) const {
        return ( k*( invT*(m_n + m_temp*invT) + (m_a1+m_a2*2.0*T)/(m_a0+m_a1*T+m_a2*sqT) - (m_b1+m_b2*2.0*T+m_b3*3.0*sqT)/(m_b0+m_b1*T+m_b2*sqT+m_b3*sqT*T) ) );
    }
   
    double n() const {return m_n;}
    
    double T() const {return m_temp;}
    
    double a0() const {return m_a0;}
    
    double a1() const {return m_a1;}
    
    double a2() const {return m_a2;}
    
    double b0() const {return m_b0;}
    
    double b1() const {return m_b1;}
    
    double b2() const {return m_b2;}
    
    double b3() const {return m_b3;}
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;    
    static std::vector<Mutation::Utilities::Units> sm_eunits;

    double m_n;
    double m_temp;
    double m_a0;
    double m_a1;
    double m_a2;
    double m_b0;
    double m_b1;
    double m_b2;
    double m_b3;
    
};

/**
 * Constant rate law (independant of temperature).
 */
class constRate : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    constRate(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    constRate(const constRate& to_copy) : m_lnA(to_copy.m_lnA) { };
    
    virtual ~constRate() { };
    
    constRate* clone() const {
        return new constRate(*this);
    }
    
    inline double getLnRate() const {
        return (m_lnA);
    }
    
    inline double derivative() const {
        return (0.0);
    }
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;
    double m_lnA;
};


/**
 * Exponential of a rational function of T.
 */
class expRat33 : public RateLaw
{
public:

    static void setUnits(const Mutation::Utilities::IO::XmlElement& node);
    
    expRat33(const Mutation::Utilities::IO::XmlElement& node, const int order);
    
    expRat33(const expRat33& to_copy)
        : m_a0(to_copy.m_a0), m_a1(to_copy.m_a1), m_a2(to_copy.m_a2), m_a3(to_copy.m_a3),
          m_b0(to_copy.m_b0), m_b1(to_copy.m_b1), m_b2(to_copy.m_b2)
    { }
    
    virtual ~expRat33() { };
    
    expRat33* clone() const {
        return new expRat33(*this);
    }
    
    inline double getLnRate(const double T) const {
        return ( (m_a0 + (m_a1 + (m_a2 + m_a3*T)*T)*T)
            /(m_b0 + (m_b1 + (m_b2 + T)*T)*T) );
    }
    
    double a0() const {return m_a0;}
    double a1() const {return m_a1;}
    double a2() const {return m_a2;}
    double a3() const {return m_a3;}
    double b0() const {return m_b0;}
    double b1() const {return m_b1;}
    double b2() const {return m_b2;}
    
private:

    static std::vector<Mutation::Utilities::Units> sm_aunits;    

    double m_a0, m_a1, m_a2, m_a3;
    double m_b0, m_b1, m_b2;
    
};

    } // namespace Kinetics
} // namespace Mutation



#endif // RATELAW_HPP
