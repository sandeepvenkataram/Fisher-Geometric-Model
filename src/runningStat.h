#ifndef RUNNINGSTAT_H
#define RUNNINGSTAT_H

#include <cmath>
/*
slightly modified class from John Cook:
http://www.johndcook.com/standard_deviation.html
implements an algorithm from Knuth 'The art of computer programing' attributted to B. P. Welford (1962)."Note on a method for calculating corrected sums of squares and products". Technometrics 4(3):419–420. 

 */
class runningStat{
public:
  runningStat(): m_n(0) {}
  void Clear(){
    m_n = 0;
  }
  void Push(double x);
  int NumDataValues() const{
    return m_n;
  }
  double Mean() const{
    return (m_n > 0) ? m_newM : 0.0;
  }
  double Variance() const{
    return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
  }
  double StandardDeviation() const{
    return sqrt( Variance() );
  }
 private:
  int m_n;
  double m_oldM, m_newM, m_oldS, m_newS;
};

#endif
