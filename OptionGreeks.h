//
//  OptionGreeks.h
//  A2.a
//  Copyright (c) 2015 Anees Attarwala. All rights reserved.
//

#ifndef A2_a_OptionGreeks_h
#define A2_a_OptionGreeks_h
#include <iostream>
#include <boost/math/distributions.hpp>
#include<cmath>
using namespace std;
using namespace boost::math;
class OptionGreeks  // base class for greeks
{
public:
    double K,r,b,T,sig;   // various option attributes.
    
    //constructor
    OptionGreeks(double _K,double _T,double _r,double _b, double _sig):K(_K),T(_T),b(_b),r(_r),sig(_sig){}
    
    //various functions
    virtual double price(double S) const{return 0;}
    virtual double delta(double S) const{return 0;}
    
    //returning vector of deltas for a vector of S.
    virtual vector <double> delta(const vector<double>& S) const{return vector<double> (0);}
    
    virtual double gamma(double S)const{return 0;}
    virtual double vega(double S) const{return 0;}
    virtual double theta(double S) const{return 0;}
       
};

class call_greeks:public OptionGreeks    // derived class for call option greeks
{
public:
    //constructor
    call_greeks(double _K,double _T,double _r,double _b, double _sig):OptionGreeks(_K,_T,_r,_b,_sig){}
    
    
    //various modified functions
    
    double price(double S) const           // Black Scholes call price formula
    {
        normal_distribution<> myNormal(0,1);
        double d1 = (log(S/K) + (b + (sig*sig)/2)*T)/sig/sqrt(T);   // BS formula
        double d2 = d1 - sig*sqrt(T);
        double Nd1 = cdf(myNormal,d1);
        double Nd2 = cdf(myNormal,d2);
        return (S*exp((b-r)*T)*Nd1-K*exp(-r*T)*Nd2);
    }
    
    double delta(double S) const{                               //call delta for single S
        normal_distribution<> myNormal(0,1);
        double nd1, d1;
        d1= (log(S/K) + (b+sig*sig/2)*T)/sig/sqrt(T);
        nd1=cdf(myNormal,d1);
        return ((exp((b-r)*T))*(nd1));
    }
                                                // function overloading
    vector<double> delta (const vector<double>& S1) const{     // call delta vector for a vector of S
        vector<double> result;
        for (int i=0;i<S1.size();++i){
            result.push_back(delta(S1[i]));}
        return result;
    }
    
    double gamma (double S) const{                        //call gamma
        normal_distribution<> myNormal(0,1);
        double n_d1, d1;
        d1= (log(S/K) + (b+sig*sig/2)*T)/sig/sqrt(T);
        n_d1=pdf(myNormal,d1);
        return (n_d1*exp((b-r)*T)/S/sig/sqrt(T));
    }
    double vega(double S) const{                          //call vega
        normal_distribution<> myNormal(0,1);
        double tmp = sig * sqrt(T);
        double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;
        double n_d1=pdf(myNormal,d1);
        return (S * exp((b-r)*T) * n_d1 * sqrt(T) );
    }
    double theta(double S)const{                         //call theta
        double tmp = sig * sqrt(T);
        normal_distribution<> myNormal(0,1);
        double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;
        double d2 = d1 - tmp;
        double n_d1= pdf(myNormal,d1);
        double Nd1=cdf(myNormal,d1);
        double Nd2 = cdf(myNormal,d2);
        double t1 =  (S * exp((b -r)*T ) * n_d1 * sig * 0.5 )/ sqrt(T);
        double t2 =  (b-r)*(S * exp((b-r)*T) * Nd1);
        double t3 =  r * K * exp(-r * T) * Nd2;
        return -(t1 + t2 + t3);
    }
 
};


class put_greeks:public OptionGreeks    // derived class for call option greeks
{
public:
    
    //constructor
    put_greeks(double _K,double _T,double _r,double _b, double _sig):OptionGreeks(_K,_T,_r,_b,_sig){}
    
    
    //various modifier methods
    double price(double S) const          // Black Scholes put price formula
    {
        normal_distribution<> myNormal(0,1);
        double d1 = (log(S/K) + (b + (sig*sig)/2)*T)/sig/sqrt(T);   // BS formula
        double d2 = d1 - sig*sqrt(T);
        double Nminusd1 = cdf(myNormal,-d1);
        double Nminusd2 = cdf(myNormal,-d2);
        return (K*exp(-r*T)*Nminusd2 -S*exp((b-r)*T)*Nminusd1);
    }
    
    double delta(double S) const{                     //put delta
        normal_distribution<> myNormal(0,1);
        double nd1, d1;
        d1= (log(S/K) + (b+sig*sig/2)*T)/sig/sqrt(T);
        nd1=cdf(myNormal,d1);
        return ((exp((b-r)*T))*(nd1-1));
    }
                                                // function overloading
    vector<double> delta (const vector<double>& S1) const{   // put delta vector for vector of underlying
        vector<double> result;
        for (int i=0;i<S1.size();++i){
            result.push_back(delta(S1[i]));}
        return result;
    }
    
    double gamma (double S) const{                    //put gamma
        normal_distribution<> myNormal(0,1);
        double n_d1, d1;
        d1= (log(S/K) + (b+sig*sig/2)*T)/sig/sqrt(T);
        n_d1=pdf(myNormal,d1);
        return (n_d1*exp((b-r)*T)/S/sig/sqrt(T));
    }
    
    double vega(double S) const{                          //put vega
        normal_distribution<> myNormal(0,1);
        double tmp = sig * sqrt(T);
        double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;
        double n_d1=pdf(myNormal,d1);
        return (S * exp((b-r)*T) * n_d1 * sqrt(T) );
    }
    double theta(double S) const{                         //put theta
        double tmp = sig * sqrt(T);
        normal_distribution<> myNormal(0,1);
        double d1 = ( log(S/K) + (b + (sig*sig)*0.5 ) * T )/ tmp;
        double d2 = d1 - tmp;
        double n_d1= pdf(myNormal,d1);
        double Nd1 = cdf(myNormal,d1);
        double Nd2 = cdf(myNormal,d2);
        double t1 =  (S * exp((b -r)*T ) * n_d1 * sig * 0.5 )/ sqrt(T);
        double t2 =  (b-r)*(S * exp((b-r)*T) * (1-Nd1));
        double t3 =  r * K * exp(-r * T) * (1-Nd2);
        return -(t1 + t2 + t3);
    }    
};

//Definitions of global functions for calculating delta, gamma, vega and theta based on divided differences

//Delta
vector<double> delta_DVD(const OptionGreeks& OG, double S,vector<double> h)
{// a, b are terms in divided difference  delta = (a-b)/2/h and gamma = (a+v(s)+b)/h^2
    vector<double> result;
    for (int i=0;i<h.size();++i){
    double a = OG.price(S-h[i]);
    double b = OG.price(S+h[i]);
    result.push_back((b-a)/2/h[i]);}
    return result;
}
//Delta
vector<double> delta_DVD(const OptionGreeks& OG, vector<double> S, double h)
{// a, b are terms in divided difference  delta = (a-b)/2/h and gamma = (a+v(s)+b)/h^2
    vector<double> result;
    for (int i=0;i<S.size();++i){
        double a = OG.price(S[i]-h);
        double b = OG.price(S[i]+h);
        result.push_back((b-a)/2/h);}
    return result;
}
//Gamma
vector<double> gamma_DVD(const OptionGreeks& OG, double S,vector <double> h)
{
    vector<double> result;
    for (int i=0;i<h.size();++i){
    double a = OG.price(S-h[i]);
    double b = OG.price(S+h[i]);
    result.push_back((b-2*OG.price(S)+a)/h[i]/h[i]);}
    return result;
}

// Vega and theta are templated on a class type - call_greeks class or a put_greeks class
// Call and put objects are passed by value as we want to change the object when computing sensitivities to volatility and time.

//Vega
template <typename T>
vector<double> vega_DVD(T OC, double S,vector<double> h,double sig1)
{
    vector<double> result;
    for (int i=0;i<h.size();++i){
    OC.sig=sig1-h[i];
    double a = OC.price(S);
    OC.sig=sig1+h[i];
    double b = OC.price(S);
    result.push_back((b-a)/2/h[i]);
    }
    return result;
}

//Theta
template <typename T1>
vector <double> theta_DVD(T1 OC1, double S1,vector<double> h, double time)
{
    vector<double> result;
    for (int i=0;i<h.size();++i){
    OC1.T=time-h[i];
    double a = OC1.price(S1);
    OC1.T=time+h[i];
    double b = OC1.price(S1);
    result.push_back(-(b-a)/2/h[i]);
    }
    return result;

}

#endif
