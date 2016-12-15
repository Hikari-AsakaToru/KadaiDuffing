#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <gl/freeglut.h>
#include <limits>
#define PI 4*std::atan(1)
#define InitRange 0.7
#define PlotRange 1.0
#define InitOffset 0.0
#define Scale 10
struct Max
{
	double x;
	double p;
};

struct CoffData{
	double Alpha;	// x'
	double Beta;	// x
	double Gunma;	// x^3
	double Delta;	// cos
	double Omega;	// wtŠp‘¬“x
	double Phi;		// wtˆÊ‘Š
	double H;		// dt
	double Wn;		// ŒÅ—LU“®”
	double Zeta;	// Œ¸Š”ä
};
class Cofficient {
protected:
	CoffData CoffRaw;
public:
	inline void SetCoff(const double dt, const double dxdt, const double x, const double x3, const double cos, const double AddFreq, const double PhiIn) {
		CoffRaw.Alpha = dxdt;
		CoffRaw.Beta = x;
		CoffRaw.Gunma = x3;
		CoffRaw.Delta = cos;
		CoffRaw.Omega = 2 * 3.1415*AddFreq;
		CoffRaw.Wn = 0.0;
		CoffRaw.Zeta = 0.0;
		CoffRaw.Phi = PhiIn;
		CoffRaw.H = dt;
	}
	inline void SetRawData(const double m, const double k, const double c, const double AddFreq,
		const double f0, const double dt, const double coffX3, const double PhiIn)
	{
		CoffRaw.Zeta = c / (2 * std::sqrt(m*k));
		CoffRaw.Wn = 2 * 3.1415*std::sqrt(k / m);
		CoffRaw.Alpha = 2 * CoffRaw.Zeta *CoffRaw.Wn;
		CoffRaw.Beta = k / m;
		CoffRaw.Gunma = coffX3*CoffRaw.Wn;
		CoffRaw.Delta = f0 / m;
		CoffRaw.Omega = 2 * 3.1415*AddFreq;
		CoffRaw.H = dt;
		CoffRaw.Phi = PhiIn;
	}
	inline double GetCoffdydx() {
		return CoffRaw.Alpha;
	}
	inline double GetCoffx() {
		return CoffRaw.Beta;
	}
	inline double GetCoffx3() {
		return CoffRaw.Gunma;
	}
	inline double GetCoffcos() {
		return CoffRaw.Delta;
	}
	inline double GetOmega() {
		return CoffRaw.Omega;
	}
	inline double Getdt() {
		return CoffRaw.H;
	}
	inline double GetCoffPhi() {
		return CoffRaw.Phi;
	}
	inline double GetWn() {
		return CoffRaw.Wn;
	}
	inline double GetZeta() {
		return CoffRaw.Zeta;
	}

	Cofficient& operator =(const CoffData& Inp) {
		this->CoffRaw = Inp;
		return *this;
	}
	Cofficient& operator =(Cofficient& Inp) {
		this->CoffRaw.Alpha = Inp.GetCoffdydx();
		this->CoffRaw.Beta	= Inp.GetCoffx();
		this->CoffRaw.Gunma = Inp.GetCoffx3();
		this->CoffRaw.Delta = Inp.GetCoffcos();
		this->CoffRaw.H		= Inp.Getdt();
		this->CoffRaw.Omega = Inp.GetOmega();
		this->CoffRaw.Phi	= Inp.GetCoffPhi();
		this->CoffRaw.Wn	= Inp.GetWn();
		this->CoffRaw.Zeta	= Inp.GetZeta();
		return *this;
	}
};
struct InputRung {
	double t;
	double x;
	double p;
}typedef  RungState;


class RungeKuttaSolve {
	std::vector<double> Calcx;
	std::vector<double> Calcdxdt;
	Cofficient CoffList;
	double Getdydx(const InputRung Inp, const double t, const double x, const double p) {
		return -CoffList.GetCoffdydx()	*(Inp.p + p)
			- CoffList.GetCoffx()		*(Inp.x + x)
			- CoffList.GetCoffx3()		*std::pow((Inp.x + x), 3.0)
			+ CoffList.GetCoffcos()		*std::cos(CoffList.GetOmega()*(Inp.t + t));
	}
public:
	RungeKuttaSolve() {
		Calcx.resize(0);
		Calcdxdt.resize(0);
	}
	inline void SetCoff( Cofficient SetCoff){
		CoffList = SetCoff;
	}
	inline void SetCoff(const CoffData SetCoff) {
		CoffList = SetCoff;
	}
	template <class ALUAUTO> inline void Init(ALUAUTO NumSeg) {
		Calcx.resize(NumSeg);
		Calcdxdt.resize(NumSeg);
	}
	inline void Init() {
		Calcx.resize(4);
		Calcdxdt.resize(4);
	}
	template <class STLAUTO> inline void CalcRungeSolve(STLAUTO* Output) {
		double H = CoffList.Getdt();
		for (auto DataPos = 0; DataPos < (*Output).size() - 1; DataPos++) {
			InputRung Status = (*Output)[DataPos];
			Calcdxdt[0] = H*Getdydx(Status, 0.0, 0.0, 0.0);
			Calcx[0] = H*(Status.p + 0.5*Calcdxdt[0]);
			for (auto i = 1; i < 3; i++) {
				Calcdxdt[i] = H*Getdydx(Status, H / 2.0, Calcx[i - 1] / 2, Calcdxdt[i - 1] / 2);
				Calcx[i] = H*(Status.p + 0.5* Calcdxdt[i]);
			}
			Calcdxdt[3] = H*Getdydx(Status, H, Calcx[2], Calcdxdt[2]);
			Calcx[3] = H*(Status.p + 0.5*Calcdxdt[3]);
			(*Output)[DataPos + 1].t = CoffList.Getdt() + Status.t;
			(*Output)[DataPos + 1].x = (Calcx[0] + 2.0*Calcx[1] + 2.0*Calcx[2] + Calcx[3]) / 6.0 + Status.x;
			(*Output)[DataPos + 1].p = (Calcdxdt[0] + 2.0*Calcdxdt[1] + 2.0*Calcdxdt[2] + Calcdxdt[3]) / 6.0 + Status.p;

		}
	}
	template <class STLAUTO> inline void CalcRungeSolveOnce(STLAUTO* Output) {
		double H = CoffList.Getdt();
		InputRung Status = (*Output);
		Output++;
		
		Calcdxdt[0] = H*Getdydx(Status, 0.0, 0.0, 0.0);
		Calcx[0]	= H*(Status.p + 0.5*Calcdxdt[0]);
		
		for (auto i = 1; i < 3; i++) {
			Calcdxdt[i] = H*Getdydx(Status, H / 2.0, Calcx[i - 1] / 2, Calcdxdt[i - 1] / 2);
			Calcx[i]	= H*(Status.p + 0.5* Calcdxdt[i]);
		}
		
		Calcdxdt[3] = H*Getdydx(Status, H, Calcx[2], Calcdxdt[2]);
		Calcx[3]	= H*(Status.p + 0.5*Calcdxdt[3]);

		(*Output).t = CoffList.Getdt() + Status.t;
		(*Output).x = (Calcx[0] + 2.0*Calcx[1] + 2.0*Calcx[2] + Calcx[3]) / 6.0 + Status.x;
		(*Output).p = (Calcdxdt[0] + 2.0*Calcdxdt[1] + 2.0*Calcdxdt[2] + Calcdxdt[3]) / 6.0 + Status.p;
	}
};

class CoffcientFunc :public Cofficient {
	RungeKuttaSolve Solver;
	CoffData dCoff;
	std::vector<RungState>  State;
	inline void WriteBackCoff() {
		if (!isnan((State.back()).x)) {
			CoffRaw.Omega = (State.back()).x;
		}
		else {
			CoffRaw.Omega = 0.0;

		}
	}
public:
	inline void UpdateCoff() {
		(*State.begin()) = (State.back());
		Solver.SetCoff(dCoff);
		Solver.CalcRungeSolve(&State);
		WriteBackCoff();
	}
	inline void SelectCalcCoff() {
		for (auto itr = State.begin(); itr != State.end();itr++) {
			(*itr).t = 0.0;
			(*itr).x = CoffRaw.Omega;
			(*itr).p = dCoff.Alpha;
		}
	}
	inline void Init(const int Num, const RungState* InitState) {
		State.resize(Num,(*InitState));
		Solver.Init();
	}
	inline void SetdCoffdfdt(const CoffData Inp) {
		dCoff = Inp;
		Solver.SetCoff(Inp);
	}
	inline void SetdCoffdfdt(Cofficient Inp) {
		dCoff.Alpha = Inp.GetCoffdydx();
		dCoff.Beta = Inp.GetCoffx();
		dCoff.Delta = Inp.GetCoffcos();
		dCoff.Gunma = Inp.GetCoffx3();
		dCoff.H = Inp.Getdt();
		dCoff.Omega = Inp.GetOmega();
		dCoff.Phi = Inp.GetCoffPhi();
		dCoff.Wn = Inp.GetWn();
		dCoff.Zeta = Inp.GetZeta();
		Solver.SetCoff(Inp);
	}

	inline void SetCoff(Cofficient Inp) {
		this->CoffRaw.Alpha = Inp.GetCoffdydx();
		this->CoffRaw.Beta = Inp.GetCoffx();
		this->CoffRaw.Delta = Inp.GetCoffcos();
		this->CoffRaw.Gunma = Inp.GetCoffx3();
		this->CoffRaw.H = Inp.Getdt();
		this->CoffRaw.Omega = Inp.GetOmega();
		this->CoffRaw.Phi = Inp.GetCoffPhi();
		this->CoffRaw.Wn = Inp.GetWn();
		this->CoffRaw.Zeta = Inp.GetZeta();
	}
};