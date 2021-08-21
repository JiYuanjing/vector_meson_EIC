//This function returns the photon flux within Q2min<Q2<Q2max and k_low<photon energy<k_high, (unit is GeV), photon energy is in proton rest framework 
//The differential photon flux is written as the function of y and Q2, see reference P. Fleischmann, PhD thesis, DESY-THESIS-2004-013
//s is the center of mass energy for proton and electron
//MJpsi is the vector meson mass (default set is J/Psi), also unit is GeV/c^2
double calPhotonFluxY(double Q2min, double Q2max, double k_low, double k_high, double s, double MJpsi=3.0969)
{
  double me = 0.000511; // electron mass
  double m_I = 0.93827; // proton mass

  double dndE = 0;
  double alpha = 1./137.036, pi = TMath::Pi();

  double minQ2 = Q2min;
  double maxQ2 = Q2max;
  //Q2 and y nsteps
  int const nQ2 = 100;
  int const ny=100;
  double lnQ2ratio = std::log(maxQ2/minQ2)/(1.*nQ2);
  double lnQ2_min = std::log(minQ2);

  for (int iQ2=0;iQ2<nQ2;iQ2++)
  {
    double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
    double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
    double mQ2 = (q2_2+q2_1)/2.; //mean
    double dQ2 = q2_2-q2_1;

    double yH = (k_high*k_high + mQ2-m_I)/(s-m_I);
    double yL = (k_low*k_low+ mQ2-m_I)/(s-m_I);
    double lnyratio = std::log(yH/yL)/(1.*ny); //in log 
    double lny_min = std::log(yL);

    for (int iy=0;iy<ny;iy++){
      double y2_1 = exp( lny_min + iy*lnyratio);	    
      double y2_2 = exp( lny_min + (iy+1)*lnyratio);
      double dy = y2_2-y2_1;
      double my = (y2_2+y2_1)/2.; //mean
      //kenimatic limit
      double min_Q2 = me*me*my*my/(1-my);
      if (mQ2<min_Q2) continue;

      dndE+=dQ2*dy*( alpha/2./ pi/mQ2 ) * ( (1+(1-my)*(1-my))/my - 2*me*me*my/mQ2 );
    }
  }
  return dndE;
}

//This function returns the photon flux within Q2min<Q2<Q2max and k_low<photon energy<k_high, (unit is GeV), photon energy is in proton rest framework 
//The differential photon flux is written as the function of photon energy k and Q2, see reference arXiv:1803.06420 
//Ee is the electron beam energy at proton rest framework
double calPhotonFluxEg(double Q2min, double Q2max, double k_low, double k_high, double Ee )
{
  double me = 0.000511;  //electron mass

  double dndE = 0;
  double alpha = 1./137.036, pi = TMath::Pi();
  
  //Q2 and photon energy nsteps
  int const nk=100;
  int const nQ2 = 100;
  double kH = k_high;
  double kL = k_low;
  double lnkratio = std::log(kH/kL)/(1.*nk); //in log 
  double lnk_min = std::log(kL);
  // double test=0;
  for (int ik=0;ik<nk;ik++){
    double k2_1 = exp( lnk_min + ik*lnkratio);	    
    double k2_2 = exp( lnk_min + (ik+1)*lnkratio);
    double dk = k2_2-k2_1;
    double mk = (k2_2+k2_1)/2.; //mean

    double min_Q2 = me*me*mk*mk/(Ee*(Ee-mk));
    double minQ2 = min_Q2>Q2min? min_Q2:Q2min;
    double maxQ2 = Q2max; // note: the maximum Q2 at kenimatic limit is 4*Ee*(Ee-k)
    double lnQ2ratio = std::log(maxQ2/minQ2)/(1.*nQ2);
    double lnQ2_min = std::log(minQ2);

    for (int iQ2=0;iQ2<nQ2;iQ2++)
    {
      double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
      double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
      double mQ2 = (q2_2+q2_1)/2.; //mean
      double dQ2 = q2_2-q2_1;

      // Integrating the effective photon flux
      dndE+=dQ2*dk*( alpha/ pi/mk/mQ2 ) * ( 1 - mk/Ee + mk*mk/2./Ee/Ee - (1-mk/Ee)*fabs(min_Q2/mQ2));
    }
  }
  return dndE;
}

