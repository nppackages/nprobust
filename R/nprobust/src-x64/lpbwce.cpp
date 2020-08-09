#include <RcppArmadillo.h>
#include <iostream> 

using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

List lpbwce(arma::vec y, arma::vec x, arma::vec K, arma::vec L, arma::vec res, double c, int p, int q, double h, double b, int deriv, int fact) {

  int N = y.n_rows;
  float rho = h / b;

  arma::vec Wp = K/h;
  arma::vec Wq = L/b;
  arma::vec Xh = (x-c)/h;
  arma::vec Xb = (x-c)/b;
  arma::mat Rq(N,q+1,arma::fill::zeros);
  arma::mat Rp(N,p+1,arma::fill::zeros);
  for (int i=0; i<q+1 ; i++) Rq.col(i) = pow(Xb,i);
  for (int i=0; i<p+1 ; i++) Rp.col(i) = pow(Xh,i);

  arma::mat dWp = diagmat(Wp);
  arma::mat dWq = diagmat(Wq);
  arma::mat dK  = diagmat(K);
  arma::mat dL  = diagmat(L);

  arma::mat Lp1 = Rp.t()*dWp*pow(Xh,p+1)/N;
  arma::mat Gp  = Rp.t()*dWp*Rp/N;
  arma::mat Gq  = Rq.t()*dWq*Rq/N;

  arma::mat iGp = Gp.i();
  arma::mat iGq = Gq.i();

  arma::mat ep1(q+1,1,arma::fill::zeros);
  ep1(p+1) = 1;
  arma::mat e0(p+1,1,arma::fill::zeros);
  e0(deriv) = fact;

  arma::mat lus0 = e0.t()*iGp*(dK*Rp).t();
  arma::mat lbc0 = lus0 - pow(rho,p+1)*(e0.t()*iGp)*Lp1*ep1.t()*iGq*(dL*Rq).t();
  arma::vec vx = pow(res,2);

  double sums2=0;
  for (int i = 0 ; i < N ; i++) {
    sums2 += pow(lbc0(i),2)*vx(i);
  }
  double s2=sums2/(N*h);

  arma::mat Krrp(p+1,p+1,arma::fill::zeros);
  arma::mat Krxip(1, p+1,arma::fill::zeros);
  arma::mat Krxp(1,  p+1,arma::fill::zeros);
  arma::mat Lrrq(q+1,q+1,arma::fill::zeros);

  for (int i = 0 ; i < N ; i++) {
    arma::mat Rpi = Rp.row(i);
    arma::mat Rqi = Rq.row(i);
    Krrp  += K(i)*Rpi.t()*Rpi;
    Lrrq  += L(i)*Rqi.t()*Rqi;
    Krxip += K(i)*Rpi*pow(Xh(i),p+1) ;
    for (int j = 0 ; j < N ; j++) {
      if (j != i) Krxp += K(i)*Rpi*pow(Xh(j),p+1);
    }
  }

  arma::mat EKrrp = Krrp/N;
  arma::mat EKrxp = Krxp/(N*(N-1));
  arma::mat EKrxip = Krxip/N;
  arma::mat ELrrq = Lrrq/N;

  arma::mat q1(1,1,arma::fill::zeros);
  arma::mat q2(1,1,arma::fill::zeros);
  arma::mat q3(1,1,arma::fill::zeros);
  arma::mat q4(1,1,arma::fill::zeros);
  arma::mat q5a(1,q+1,arma::fill::zeros);
  arma::mat q5b(q+1,1,arma::fill::zeros);
  arma::mat q6(1,1,arma::fill::zeros);
  arma::mat q7a(1,q+1,arma::fill::zeros);
  arma::mat q7b(q+1,q+1,arma::fill::zeros);
  arma::mat q7c(q+1,1,arma::fill::zeros);
  arma::mat q8(1,1,arma::fill::zeros);
  arma::mat q9(1,1,arma::fill::zeros);
  arma::mat q10(1,1,arma::fill::zeros);
  arma::mat q11(1,1,arma::fill::zeros);
  arma::mat q12(1,1,arma::fill::zeros);
  arma::mat q3a(1,1,arma::fill::zeros);

  for (int i = 0 ; i < N ; i++) {
    arma::mat Rpi = Rp.row(i);
    arma::mat Rqi = Rq.row(i);

    q1 += pow(lbc0(i)*res(i),3);

    arma::mat lus1 = fact * iGp.row(deriv) * (EKrrp - K(i)*Rpi.t()*Rpi)*iGp*K(i)*Rpi.t();
    arma::mat T1   = fact * iGp.row(deriv) * ((EKrrp - K(i)*Rpi.t()*Rpi)*iGp*Lp1*ep1.t())    *(iGq*L(i)*Rqi.t());
    arma::mat T2   = fact * iGp.row(deriv) * ((K(i)*Rpi*pow(Xh(i),p+1) - EKrxip).t())*ep1.t()*(iGq*L(i)*Rqi.t());
    arma::mat T3   = fact * iGp.row(deriv) * ((Lp1*ep1.t()*iGq)*(ELrrq - L(i)*Rqi.t()*Rqi))  *(iGq*L(i)*Rqi.t());
    arma::mat lbc1 = lus1 - pow(rho,p+1)*(T1 + T2 + T3);

    q2 += lbc1*lbc0(i)*pow(res(i),2);

    q3 += pow(lbc0(i),4)*(pow(res(i),4)-pow(vx(i),2));

    q4 += pow(lbc0(i),2)*(Rqi*iGq*L(i)*Rqi.t())*pow(res(i),2);

    q5a += pow(lbc0(i),3)*(Rqi*iGq)*pow(res(i),2);
    q5b += L(i)*Rqi.t()*lbc0(i)*pow(res(i),2);

    q7a += lbc0(i)*pow(res(i),2)*L(i)*Rqi*iGq;
    q7b += pow(lbc0(i),2)*Rqi.t()*Rqi*iGq;
    q7c += lbc0(i)*pow(res(i),2)*L(i)*Rqi.t();

    q8  += pow(lbc0(i)*res(i),4);
    q9  += (pow(lbc0(i),2)*vx(i)-h*s2)*pow(lbc0(i)*res(i),2);

    q12 += pow(pow(lbc0(i),2)*vx(i)-h*s2,2);

    q3a += pow(lbc0[i]*res[i],3);


    for (int j = 0 ; j < N ; j++) {
      if (j != i) {
        arma::mat Rpj = Rp.row(j);
        arma::mat Rqj = Rq.row(j);

        arma::mat lus1 = fact * iGp.row(deriv) *  (EKrrp - K(j)*Rpj.t()*Rpj)*iGp*K(i)*Rpi.t();
        arma::mat T1   = fact * iGp.row(deriv) * ((EKrrp - K(j)*Rpj.t()*Rpj)*iGp*Lp1*ep1.t())    *(iGq*L(i)*Rqi.t());
        arma::mat T2   = fact * iGp.row(deriv) * ((K(j)*Rpj*pow(Xh(i),p+1) - EKrxp).t())*ep1.t() *(iGq*L(i)*Rqi.t());
        arma::mat T3   = fact * iGp.row(deriv) * ((Lp1*ep1.t()*iGq)*(ELrrq - L(j)*Rqj.t()*Rqj))  *(iGq*L(i)*Rqi.t());
        arma::mat lbc1 = lus1 - pow(rho,p+1)*(T1 + T2 + T3);

        q10 += lbc1*lbc0(i)*pow(lbc0(j)*res(j),2)*vx(i);
        q11 += lbc1*lbc0(i)*(pow(lbc0(j),2)*vx(j)-h*s2)*pow(res(i),2);

        q6 += pow(lbc0(i),2)*pow(Rqi*iGq*L(j)*Rqj.t(),2)*pow(res(j),2) ;
       }
    }
  }

  arma::mat Eq1  = pow(q1/(N*h),2);
  arma::mat Eq2  = q2/(N*h);
  arma::mat Eq3  = q3/(N*h);
  arma::mat Eq4  = q4/(N*h);
  arma::mat Eq5  = (q5a/(N*h))*(q5b/(N*h));
  arma::mat Eq6  = q6/(N*(N-1)*pow(h,2));
  arma::mat Eq7  = (q7a/(N*h))*(q7b/(N*h))*(q7c/(N*h));
  arma::mat Eq8  = q8/(N*h);
  arma::mat Eq9  = q9/(N*h);
  arma::mat Eq10 = q10/(N*(N-1)*pow(h,2));
  arma::mat Eq11 = q11/(N*(N-1)*pow(h,2));
  arma::mat Eq12 = q12/(N*h);

  double z  = 1.959964;
  double pz = 0.05844507;

  arma::mat  q1bc =   pz*(
       Eq1*(pow(z,3)/3+7*z/4+s2*z*(pow(z,2)-3)/4)/pow(s2,3)
    +  Eq2*(-z*(pow(z,2)-3)/2)/s2
    +  Eq3*(z*(pow(z,2)-3)/8)/pow(s2,2)
    -  Eq4*(z*(pow(z,2)-1)/2)/s2
    -  Eq5*(z*(pow(z,2)-1))/pow(s2,2)
    +  Eq6*(z*(pow(z,2)-1)/4)/s2
    +  Eq7*(z*(pow(z,2)-1)/2)/pow(s2,2)
    +  Eq8*(-z*(pow(z,2)-3)/24)/pow(s2,2)
    +  Eq9*(z*(pow(z,2)-1)/4)/pow(s2,2)
    +  Eq10*(z*(pow(z,2)-3))/pow(s2,2)
    +  Eq11*(-z)/pow(s2,2)
    +  Eq12*(-z*(pow(z,2)+1)/8)/pow(s2,2)
    );

    double q2bc = -pz*z/(2*s2);

    arma::mat Eq3a  = q3a/(N*h);
    arma::mat q3bc = pz*Eq3a/(pow(s2,2))*(pow(z,3)/3);

 //std::cout << "The matrix q3a is:\n" << q3a << "\n\n"; PROBLEMATIC
 //std::cout << "The matrix N is:\n" << N << "\n\n"; DONE
 //std::cout << "The matrix h is:\n" << h << "\n\n"; DONE


 //std::cout << "The matrix pz is:\n" << pz << "\n\n"; DONE
 //std::cout << "The matrix Eq3a is:\n" << Eq3a << "\n\n"; PROBLEMATIC
 //std::cout << "The matrix s2 is:\n" << s2 << "\n\n"; DONE
 //std::cout << "The matrix z is:\n" << z << "\n\n"; DONE

 //std::cout << "The matrix q1bc is:\n" << q1bc << "\n\n"; PROBLEMATIC
 //std::cout << "The matrix q2bc is:\n" << q2bc << "\n\n"; DONE
 //std::cout << "The matrix q3bc is:\n" << q3bc << "\n\n"; PROBLEMATIC


    arma::mat q1rbc = 2*q1bc/pz;
    double    q2rbc = 2*q2bc/pz;
    arma::mat q3rbc = 2*q3bc/pz;


  return Rcpp::List::create(
    Rcpp::Named("q1rbc") = q1rbc,
    Rcpp::Named("q2rbc") = q2rbc,
    Rcpp::Named("q3rbc") = q3rbc
  ) ;

}

