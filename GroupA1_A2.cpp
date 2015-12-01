//
//  main.cpp
//  GroupA

//  Copyright (c) 2015 Anees Attarwala. All rights reserved.
//  Exact call and put option prices on following batches
/*Batch 1: T = 0.25, K = 65, sig = 0.30, r = 0.08, S = 60 (then C = 2.13293, P = 5.84584).
 Batch 2: T = 1.0, K = 100, sig = 0.2, r = 0.0, S = 100 (then C = 7.96632, P = 7.96632).
 Batch 3: T = 1.0, K = 10, sig = 0.50, r = 0.12, S = 5 (C = 0.204121, P = 4.0733).
 Batch 4: T = 30.0, K = 100.0, sig = 0.30, r = 0.08, S = 100.0 (C = 92.1749, P = 1.24651).*/

#include <iostream>
#include <cmath>
#include"OptionGreeks.h"
#include <boost/math/distributions.hpp> // For non-member functions of distributions
#include <boost/math/distributions/normal.hpp>
using namespace std;
//Global Call option function definition
double Call (double& S, double& K, double& T, double& b, double& r, double& sigma)
{
    
    using namespace boost::math;
    normal_distribution<> myNormal(0,1);
    double d1 = (log(S/K) + (b + (sigma*sigma)/2)*T)/sigma/sqrt(T);   // BS formula
    double d2 = d1 - sigma*sqrt(T);
    double Nd1 = cdf(myNormal,d1);
    double Nd2 = cdf(myNormal,d2);
    return (S*exp((b-r)*T)*Nd1-K*exp(-r*T)*Nd2);
}
//Global Put option function definition
double Put (double& S, double& K, double& T, double& b, double& r, double& sigma)
{
    
    using namespace boost::math;
    normal_distribution<> myNormal(0,1);
    double d1 = (log(S/K) + (b + (sigma*sigma)/2)*T)/sigma/sqrt(T);  // BS Formula
    double d2 = d1 - sigma*sqrt(T);
    double Nminusd1 = cdf(myNormal,-d1);
    double Nminusd2 = cdf(myNormal,-d2);
    return (K*exp(-r*T)*Nminusd2-S*exp((b-r)*T)*Nminusd1);
}

// Encapsulate all data in one place
struct OptionData
{ // Option data + behaviour (all public!)
    double K;
    double T;
    double r;
    double sig;
    // Extra data later ...
};

//Global Call option function definition
double Call (const OptionData& OD, double& S)
{
    
    using namespace boost::math;
    normal_distribution<> myNormal(0,1);
    double d1 = (log(S/(OD.K)) + (OD.r + (OD.sig*OD.sig)/2)*OD.T)/(OD.sig)/sqrt(OD.T);
    double d2 = d1 - OD.sig*sqrt(OD.T);
    double Nd1 = cdf(myNormal,d1);
    double Nd2 = cdf(myNormal,d2);
    return (S*Nd1-OD.K*exp(-OD.r*OD.T)*Nd2);
}
//Global Put option function definition
double Put (const OptionData& OD, double& S)
{
    
    using namespace boost::math;
    normal_distribution<> myNormal(0,1);
    double d1 = (log(S/(OD.K)) + (OD.r + (OD.sig*OD.sig)/2)*OD.T)/(OD.sig)/sqrt(OD.T);
    double d2 = d1 - OD.sig*sqrt(OD.T);
    double Nminusd1 = cdf(myNormal,-d1);
    double Nminusd2 = cdf(myNormal,-d2);
    return (OD.K*exp(-OD.r*OD.T)*Nminusd2-S*Nminusd1);
}

//Global Mesh functions to return put and call option prices as functions of S.
//Demonstrating function overloading as function names are still same but they return vector of prices.
vector<double> Call(const OptionData& OP, vector<double> S)
{
    vector <double>optionprice;
    for (int i=0;i<S.size();++i)
    {
        optionprice.push_back(Call(OP,S[i]));
    }
    return optionprice;
}

vector<double> Put(const OptionData& OP, vector<double> S)
{
    vector <double>optionprice;
    for (int i=0;i<S.size();++i)
    {
        optionprice.push_back(Put(OP,S[i]));
        
    }
    return optionprice;
}


//Global functions that return put and call option prices as functions of volatility.
// passed by value as we need to change object data as we vary the volatility

vector<double> vol_variation(OptionData OP, vector<double> vol1, double S,int type)
{
    vector <double>optionprice;
    for (int i=0;i<vol1.size();++i)
    {
        OP.sig=vol1[i];
        if (type==1)
            optionprice.push_back(Call(OP,S));
        else if (type==-1)
            optionprice.push_back(Put(OP,S));
    }
    return optionprice;
}
// passed by value as we need to change object data as we vary the time to expiry
vector<double> exp_variation(OptionData OP, vector<double> exp1, double S, int type)
{
    vector <double>optionprice;
    for (int i=0;i<exp1.size();++i)
    {
        OP.T=exp1[i];
        if (type==1)
            optionprice.push_back(Call(OP,S));
        else if (type==-1)
            optionprice.push_back(Put(OP,S));
    }
    return optionprice;
}


// Global function to print a vector <double>
void print(const vector<double>& v)
{
    vector<double>::const_iterator it;
    for (it=v.begin();it!=v.end();++it)
    {
        cout<<*it<<" ";
    }
    
}
// Mesher returns a vector <double> array of mesh points
vector<double> mesher (double start, double end, double size)
{
    double n_steps = (end - start) / double (size);  //no of steps in mesh
    double temp=start;
    vector<double> result;
    result.push_back(start);
    for (int j = 1; j < n_steps+1; j++)
    {
        temp = temp+size;
        result.push_back(temp);
    }
    return result;
}




int main()
{//A1a) Computing prices per exact formula for call and put options
cout<<"***********BEGINNING GROUP A1******************************************"<<endl;
cout<<"A1a) Computing prices per exact formula for call and put options"<<endl;

   // creating 4 batches as vectors of T,K,Sig,r and S;
    vector<double> T={0.25,1.0,1.0,30.0};
    vector<double> K={65,100,10,100.0};
    vector<double> sig={0.30,0.2,0.50,0.30};
    vector<double> r={0.08,0.0,0.12,0.08};
    vector<double> S={60,100,5,100.0};
    // vectors for storing call and put prices for 4 batch runs
    vector <double>callprice;
    vector <double>putprice;
    for (int i=0;i<4;++i)
    {
        callprice.push_back(Call(S[i],K[i],T[i],r[i],r[i],sig[i]));
        putprice.push_back(Put(S[i],K[i],T[i],r[i],r[i],sig[i]));
    }
    cout<<"Call option prices below :"<<endl;
    for (int i=0; i<4; ++i)
    {
        cout<<callprice[i]<<endl;
    }
    cout<<"Put option prices below :"<<endl;
    for (int i=0; i<4; ++i)
    {
        cout<<putprice[i]<<endl;
    }
    
    //A1b) Testing put call parity. Given call price computing put price using put call parity and testing put price from formula implemented in a).
  
    cout<<"*********************************************************************"<<endl;
    cout<<"A1b) Testing put call parity. Given call price computing put price using put call parity and testing put price from formula implemented in a)."<<endl;
    
      //vector for storing put results computed from put call parity relation
    vector<double> PutFromParity;
    cout<<"Put option prices from put call parity"<<endl;
    for (int i=0;i<4;++i)
    {
        PutFromParity.push_back(callprice[i] + K[i]*exp(-r[i]*T[i])-S[i]);
        cout<<PutFromParity[i]<<endl;
    }
    cout<<"Note the Put price from put call parity relation is same as put price from formula!"<<endl<<endl;
    
    //A1c) Repeating part A1a) and A1b) using OD
    cout<<"*********************************************************************"<<endl;
    cout<<"A1c) Repeating part A1a) and A1b) using option structure to encapsulate option data "<<endl<<endl;
    
    vector<OptionData> OD_call;  //vector of call option structs
    vector<OptionData> OD_put;  //vector of put option structs
    OptionData calloption;    //structures to first set values, and then pass it onto
    OptionData putoption;     // vector of structures using push_back.
    
    
    for (int i=0;i<4;++i)         //filling in option data call and put vector of structures
    {
        calloption.K=K[i];
        calloption.T=T[i];
        calloption.r=r[i];
        calloption.sig=sig[i];
        OD_call.push_back(calloption);
        
        putoption.K=K[i];
        putoption.T=T[i];
        putoption.r=r[i];
        putoption.sig=sig[i];
        OD_put.push_back(putoption);
    }
    
    // vectors for storing callprice and putprice
    
    vector <double>callprice_c;
    vector <double>putprice_c;
    cout <<"CALL AND PUT PRICES USING STRUCTURE OPTION DATA (Repeating A1a) with option structure)"<<endl;
    
    for (int i=0;i<4;++i)
    {
        callprice_c.push_back(Call(OD_call[i],S[i]));
        putprice_c.push_back(Put(OD_put[i],S[i]));
    }
    cout<<"Call option prices below :"<<endl;
    for (int i=0; i<4; ++i)
    {
        cout<<callprice_c[i]<<endl;
    }
    cout<<"Put option prices below :"<<endl;
    for (int i=0; i<4; ++i)
    {
        cout<<putprice_c[i]<<endl;
    }
    
    // Testing put call parity. Given call price computing put price using put call parity and testing put price from formula implemented in a).
    //*********************************************************************
    //vector for storing put results computed from put call parity relation
    vector<double> PutFromParity_c;
    cout<<"Put option prices from put call parity (repeating b) with option structure)"<<endl;
    for (int i=0;i<4;++i)
    {
        PutFromParity_c.push_back(callprice_c[i] + OD_call[i].K*exp(-OD_call[i].r*OD_call[i].T)-S[i]);
        cout<<PutFromParity_c[i]<<endl;
    }
    cout<<"Note the Put price from put call parity relation is same as put price from formula!"<<endl<<endl;

    cout<<"*********************************************************************"<<endl<<endl;
    cout<<"A1d) Variation of Call and Put price for batch 4) with equity"<<endl;
    
    // Computing prices per exact formula for call and put options
    
    // Creating options data call and put objects
    
    OptionData OD_call_d;  // call option struct
    OptionData OD_put_d;  //put option struct
    
    //filling in Batch 4 option data for call and put options.
    OD_call_d.K=100;
    OD_call_d.T=30;
    OD_call_d.r=.08;
    OD_call_d.sig=0.3;
    
    OD_put_d.K=100;
    OD_put_d.T=30;
    OD_put_d.r=.08;
    OD_put_d.sig=0.3;
    
    // vectors for storing callprice and putprice

    vector<double> equity=mesher(75,125,5);  // creating equity mesh of underlying prices
    vector<double> callprice_d;
    vector <double>putprice_d;
    callprice_d = Call(OD_call_d,equity);  // calling Call and Put price functions
    putprice_d = Put(OD_put_d,equity);
    cout.precision(4);
    cout <<"VARIATION OF CALL AND PUT PRICES"<<endl;
    cout << "Equity varition "<<endl;
    print(equity);
    cout<<endl<<"*****CALL VARIATION*****************"<<endl;
    print(callprice);
    cout<<endl<<"*****PUT VARIATION******************"<<endl;
    print(putprice);
    cout<<endl<<endl;
    cout<<"We note that option price increases with increase in underlying equity price."<<endl;
    
    cout<<"*********************************************************************"<<endl<<endl;
    cout<<"A1e) Variation of Call and Put price for batch 4) with volatility and expiry"<<endl;

    vector<double> vol=mesher(.1,.5,.05);  // creating volatility mesh
    vector<double> expiry1=mesher(10,50.0,5);  // creating expiry mesh
    double S_e = 100; // Batch 4 S=100
    vector<double> call_vol;  // variation of call and put with volatility
    vector <double>put_vol;
    vector<double> call_exp; // variation of call and put with expiry
    vector <double>put_exp;
    
    call_vol = vol_variation(OD_call_d,vol,S_e,1); // 1 for call; calling Call and Put price functions
    put_vol = vol_variation(OD_put_d,vol,S_e,-1);  //-1 for put
    call_exp = exp_variation(OD_call_d,expiry1,S_e,1); // 1 for call; calling Call and Put price functions
    put_exp = exp_variation(OD_put_d,expiry1,S_e,-1); //-1 for put
    
    cout.precision(4);
    cout <<"VARIATION OF CALL AND PUT PRICES"<<endl<<endl;
    cout << "Volatility variation: "<<endl;
    print(vol);
    cout<<endl<<"*****CALL VARIATION WITH VOLATILITY*****************"<<endl;
    print(call_vol);
    cout<<endl<<"*****PUT VARIATION WITH VOLATILTIY******************"<<endl;
    print(put_vol);
    cout<<endl<<endl;
    cout << "Expiry variation: "<<endl;
    print(expiry1);
    cout<<endl<<"*****CALL VARIATION WITH EXPIRY*****************"<<endl;
    print(call_exp);
    cout<<endl<<"*****PUT VARIATION WITH EXPIRY******************"<<endl;
    print(put_exp);
    cout<<endl<<endl;
    cout<<"We note that call option prices increase with increasing volatility and increasing time to expiry. For put options the variation with volatility is monotonic increasing but variation with respect to time to expiry can be decreasing or increasing. "<<endl;

    
    
    cout<<"*********************************************************************"<<endl<<endl;
    cout<<"*****************BEGINNING GROUP A2*********************************"<<endl<<endl;
    double K1=100, T1=.5, S1=105,r1=0.1, b1=0, sig1=0.36;  //option paramaters
    
    //pointing derived class objects using base class pointers.
    
    boost::shared_ptr<OptionGreeks> call (new call_greeks(K1,T1,r1,b1,sig1));  //using shared pointers
    boost::shared_ptr<OptionGreeks> put (new put_greeks(K1,T1,r1,b1,sig1));
    cout<<"SOLUTION TO A2a)"<<endl<<endl;
    cout<<"Call delta is :"<<call->delta(S1)<<endl;
    cout<<"Call gamma is :"<<call->gamma(S1)<<endl;
    cout<<"Call vega is :"<<call->vega(S1)<<endl;
    cout<<"Call theta is :"<<call->theta(S1)<<endl;
    cout<<"****************************************"<<endl;
    
    cout<<"Put delta is :"<<put->delta(S1)<<endl;
    cout<<"Put gamma is :"<<put->gamma(S1)<<endl;
    cout<<"Put vega is :"<<put->vega(S1)<<endl;
    cout<<"Put theta is :"<<put->theta(S1)<<endl<<endl;
    cout<<"****************************************"<<endl<<endl;
    // part A2b) calcuation put and call variation with equity prices
    vector<double> equity2;          // vector to store equity mesh
    vector<double> CallDeltaValues;       // vector to store call deltas
    vector<double> PutDeltaValues;       // vector to store put deltas
    
    equity2=mesher(80,120,5);       // creating equity mesh
    
    cout<<"SOLUTION TO A2b)"<<endl<<endl;
    cout<<"Equity variation :"<<endl;
    print(equity2);
    cout<<endl<<endl;
    
    CallDeltaValues = call->delta(equity2);  //calculate call delta and store in vector.
    PutDeltaValues = put->delta(equity2);  //calculate put delta and store in vector.
    
    //print results
    cout<<"Call delta variation with underlying equity :"<<endl;
    cout<<setprecision(4);
    print(CallDeltaValues);
    cout<<endl<<endl;
    cout<<"Put delta variation with underlying equity :"<<endl;
    print(PutDeltaValues);
    cout<<endl;
    cout<<"*********************************************************************"<<endl<<endl;
    cout<<"SOLUTION TO A2c)"<<endl;
    call_greeks call1(K1,T1,r1,b1,sig1);  //call object to store vega and theta
    put_greeks put1(K1,T1,r1,b1,sig1);  //put object to store vega and theta
    vector<double> calldelta_DVD;             //call delta from divided difference
    vector<double> callgamma_DVD;             //call gamma from divided difference
    vector<double> callvega_DVD;             //call vega from divided difference
    vector<double> calltheta_DVD;             //call theta from divided difference
    vector<double> putdelta_DVD;             //put delta from divided difference
    vector<double> putgamma_DVD;             //put gamma from divided difference
    vector<double> putvega_DVD;             //put vega from divided difference
    vector<double> puttheta_DVD;             //put theta from divided difference
    vector<double> h{0.2,0.1,.01,.001};      // size paramater for divided difference
 
    // performing divided difference to calculate delta for various "h" sizes
    // Calling global functions for divided difference defined in OptionGreeks.h
    calldelta_DVD = delta_DVD(*call,S1,h);
    callgamma_DVD = gamma_DVD(*call,S1,h);
    callvega_DVD = vega_DVD<call_greeks>(call1,S1,h,sig1);
    calltheta_DVD = theta_DVD<call_greeks>(call1,S1,h,T1);
    putdelta_DVD = delta_DVD(*put,S1,h);
    putgamma_DVD = gamma_DVD(*put,S1,h);
    putvega_DVD = vega_DVD<put_greeks>(put1,S1,h,sig1);
    puttheta_DVD = theta_DVD<put_greeks>(put1,S1,h,T1);
    cout<<endl;
    cout<<"***REDOING PART A2a) by Divided differences for Call option***"<<endl<<endl;
    cout<<"Various values of h are :"<<endl;
    print(h);
    cout<<endl;
    cout<<"Divided Difference Call Deltas :"<<endl;
    print(calldelta_DVD);
    cout<<endl;
    cout<<"Divided Difference Call Gammas :"<<endl;
    print(callgamma_DVD);
    cout<<endl;
    cout<<"Divided Difference Call  Vega:"<<endl;
    print(callvega_DVD);
    cout<<endl;
    cout<<"Divided Difference Call Theta :"<<endl;
    print(calltheta_DVD);
    cout<<endl<<endl<<endl;
    cout<<"***REDOING PART A2a) by Divided differences for Put option***"<<endl;
    cout<<"Various values of h are :"<<endl;
    print(h);
    cout<<endl;
    cout<<"Divided Difference Put Deltas :"<<endl;
    print(putdelta_DVD);
    cout<<endl;
    cout<<"Divided Difference Put Gammas :"<<endl;
    print(putgamma_DVD);
    cout<<endl;
    cout<<"Divided Difference Put  Vega:"<<endl;
    print(putvega_DVD);
    cout<<endl;
    cout<<"Divided Difference Put Theta :"<<endl;
    print(puttheta_DVD);
    cout<<endl;
    cout<<"Clearly as h decreases divided difference value approaches the formulaic value"<<endl<<endl;
    //vector to store call and put delta variation vs equity for different 'h'
    vector<double>CallDeltaValues_DVD,PutDeltaValues_DVD;
    cout<<setprecision(4);
    cout<<"*********REDIONG PART A2b) by Divided differences *************"<<endl;
    cout<<endl;
    //call and put delta variation as a function of S using divided differences
    // Calling global functions for divided difference defined in OptionGreeks.h
    cout<<"Call and put option delta as a function of S by divided difference "<<endl;
    for (int i=0;i<h.size();++i)
    {
        cout<<"CALL OPTION VARIATION WITH EQUITY USING h = "<<h[i]<<endl;
        CallDeltaValues_DVD = delta_DVD(*call,equity,h[i]);
        print(CallDeltaValues_DVD);
        cout<<endl;
    }
    for (int i=0;i<h.size();++i)
    {
        cout<<"PUT OPTION VARIATION WITH EQUITY USING h = "<<h[i]<<endl;
        PutDeltaValues_DVD = delta_DVD(*put, equity,h[i]);
        print(PutDeltaValues_DVD);
        cout<<endl;
    }
   
    return 0;
}

