#include "Common.h"
#include "StringOpenGL.h"
class RungeKutta;
RungeKutta	*Rung_global_pointer;
GLFONT		*font;
int winID[3];
#define glAxisColor			glColor3d(1.0,0.0,0.0)
#define glProtColor			glColor3d(0.0,1.0,0.0)
#define glBackGroundColor	glColor3d(0.0,0.0,0.0)
#define glScaleColor		glColor3d(1.0,1.0,1.0)
#define glScalePlotColor	glColor3d(0.5,0.5,0.5)


struct Max
{
	double x;
	double p;
};
struct Cofficient {
protected:
	double Alpha;	// x'
	double Beta;	// x
	double Gunma;	// x^3
	double Delta;	// cos
	double Omega;	// wtäpë¨ìx
	double dOmegadt;	// wtäpë¨ìx
	double Phi;		// wtà ëä
	double H;		// dt
	double Wn;		// å≈óLêUìÆêî
	double Zeta;	// å∏êäî‰
public:
	inline void SetCoff(const double dt, const double dxdt, const double x, const double x3, const double cos, const double AddFreq, const double PhiIn) {
		Alpha = dxdt;
		Beta = x;
		Gunma = x3;
		Delta = cos;
		Omega = 2 * 3.1415*AddFreq;
		Wn = 0.0;
		Zeta = 0.0;
		Phi = PhiIn;
		H = dt;
		dOmegadt = 0.0;
	}
	inline void SetRawData(const double m, const double k, const double c, const double AddFreq,
				const double f0 , const double dt, const double coffX3, const double PhiIn)
	{
		Zeta = c/(2*std::sqrt(m*k));
		Wn = 2*3.1415*std::sqrt(k / m);
		Alpha = 2 * Zeta *Wn;
		Beta = k / m;
		Gunma = coffX3*Wn;
		Delta = f0 / m;
		Omega = 2*3.1415*AddFreq;
		H = dt;
		Phi = PhiIn;
		dOmegadt = 0.0;
	}
	inline double GetCoffdydx() {
		return Alpha;
	}
	inline double GetCoffx() {
		return Beta;
	}
	inline double GetCoffx3() {
		return Gunma;
	}
	inline double GetCoffcos() {
		return Delta;
	}
	inline double GetOmega() {
		return Omega;
	}
	inline double Getdt() {
		return H;
	}
	inline void UpdateParam() {
		this->Omega += dOmegadt*H;
	}
	inline void SetTimeOmega(double dt) {
		this->dOmegadt = dt;
	}
};
struct InputRung {
	double t;
	double x;
	double p;
}typedef  RungState;
#define PI 4*std::atan(1)
class RungeKutta {
	std::vector<RungState> Output;
	std::vector<double> Calcx;
	std::vector<double> Calcdxdt;
	Cofficient CoffList;
	InputRung InitStatus;
	Max MaxPara;
	double Getdydx(const InputRung Inp,const double t, const double x, const double p) {
		return -CoffList.GetCoffdydx()*(Inp.p + p) 
			+ CoffList.GetCoffx()*(Inp.x + x) 
			- CoffList.GetCoffx3()*std::pow((Inp.x + x),3.0)
			+ CoffList.GetCoffcos()*std::cos(CoffList.GetOmega()*(Inp.t+t));
	}
	inline void DrawScalePlotLineX() {
		for (auto i = 0; i < Scale; i++) {
			glScalePlotColor;
			glBegin(GL_LINES);
			glVertex2d((double)i / (Scale/2 ) - 1.0, -1.0);
			glVertex2d((double)i / (Scale/2 ) - 1.0, 1.0);
			glEnd();
		}
	}
	inline void DrawScalePlotLineY() {
		for (auto i = 0; i < Scale; i++) {
			glScalePlotColor;
			glBegin(GL_LINES);
			glVertex2d(-1.0, (double)i / (Scale / 2)-1.0);
			glVertex2d(1.0, (double)i / (Scale / 2) - 1.0);
			glEnd();
		}
	}
	inline void DrawScaleX(double (RungeKutta::*TypeScaleX)(int)) {
			std::string Tmp;
		for (auto i = 0; i < Scale; i++) {
			glScaleColor;
			double ValX = (this->*TypeScaleX)((int)i);
			Tmp = std::to_string(ValX);
			if (ValX < 0) {
				Tmp.erase(Tmp.begin() + 5, Tmp.end());
			}
			else
			{
				Tmp.erase(Tmp.begin() + 4, Tmp.end());
			}
			glRasterPos2d(-1.0+i*0.2, 0);
			glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, reinterpret_cast<const unsigned char*>(Tmp.c_str()));
		}
		DrawScalePlotLineX();
	}
	inline void DrawScaleY(double (RungeKutta::*TypeScaleY)(int )) {
		std::string Tmp;
		for (auto i = 0; i < Scale; i++) {
			if (i!=5) {
				glScaleColor;
				double ValY = (this->*TypeScaleY)((int)i);
				Tmp = std::to_string(ValY);
				glRasterPos2d(0.05, -1.0 + i*0.2);
				if (ValY < 0) {
					Tmp.erase(Tmp.begin() + 5, Tmp.end());
				}
				else
				{
					Tmp.erase(Tmp.begin() + 4, Tmp.end());
				}
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, reinterpret_cast<const unsigned char*>(Tmp.c_str()));
			}
		}
		DrawScalePlotLineY();


	}
	inline void DrawScale(double (RungeKutta::*TypeScaleX)(int i), double (RungeKutta::*TypeScaleY)(int i)) {
		DrawScaleX(TypeScaleX);
		DrawScaleY(TypeScaleY);
	}
	inline double GetScaleX(int i) {
		return MaxPara.x / Range*(i - Scale / 2);
	}
	inline double GetScaleT(int i) {
		return (Output.size()*CoffList.Getdt() / 2.0)*i*2 ;
	}
	inline double GetScaleP(int i) {
		return MaxPara.p / Range*(i - Scale / 2);
	}
public:
	template<class TYPE,class NUM> RungeKutta(const TYPE dt, const TYPE TRange, const NUM NumSeg) {
		Output.resize((int)(TRange/dt+1)) ;
		Calcx.resize(NumSeg);
		Calcdxdt.resize(NumSeg);
		
	}
	inline void SetParam(const Cofficient SetCoff,  const InputRung Inp) {
		CoffList = SetCoff;
		InitStatus = Inp;
	}
	inline void CalcDffingEqudydx() {
		Output[0] = InitStatus;
		double H = CoffList.Getdt();
		for (auto DataPos = 0; DataPos < Output.size()-1;DataPos++) {
			InputRung Status = Output[DataPos];
			Calcdxdt[0] = H*Getdydx(Status, 0.0, 0.0, 0.0);
			Calcx[0] = H*(Status.p + 0.5*Calcdxdt[0]);
			for (auto i = 1; i < 3; i++) {
				Calcdxdt[i] = H*Getdydx(Status, H / 2.0, Calcx[i - 1] / 2, Calcdxdt[i - 1] / 2);
				Calcx[i] = H*(Status.p + 0.5* Calcdxdt[i]);
			}
			Calcdxdt[3] = H*Getdydx(Status, H, Calcx[2], Calcdxdt[2]);
			Calcx[3] = H*(Status.p + 0.5*Calcdxdt[3]);
			Output[DataPos + 1].t = CoffList.Getdt() + Status.t;
			Output[DataPos + 1].x = (Calcx[0] + 2.0*Calcx[1] + 2.0*Calcx[2] + Calcx[3]) / 6.0 + Status.x;
			Output[DataPos + 1].p = (Calcdxdt[0] + 2.0*Calcdxdt[1] + 2.0*Calcdxdt[2] + Calcdxdt[3]) / 6.0 + Status.p;
			
		}
	}
	inline void OutputCSV() {
		std::ofstream OutputFile("RungKutta001.csv",std::ofstream::trunc|std::ofstream::out);
		OutputFile << "t,x,p" << std::endl;
		//for (auto Rungitr = Output.begin(); Rungitr != Output.end(); Rungitr++) {
		//	OutputFile << (*Rungitr).t << "," << (*Rungitr).x << "," << (*Rungitr).p << std::endl;
		//}
		for (auto Rungitr =0; Rungitr < Output.size(); Rungitr++) {
			OutputFile << Output[Rungitr].t << "," << Output[Rungitr].x << "," << Output[Rungitr].p << std::endl;
		}


	}
	inline void MakeGraphData() {
		auto Maxx = Output[0].x;
		auto Maxp = Output[0].p;

		for (auto i = 1; i < Output.size(); i++) {
			if (Maxx < Output[i].x) 
			{
				Maxx = Output[i].x;
			}
			if (Maxp < Output[i].p)
			{
				Maxp = Output[i].p;
			}
		}
		for (auto i = 0; i < Output.size(); i++) {
			Output[i].x = Output[i].x / (Maxx / Range);
			Output[i].p = Output[i].p / (Maxp / Range);
			Output[i].t = Output[i].t / (Output.size()*CoffList.Getdt()/2.0)- 1.0;
		}
		this->MaxPara.x = Maxx;
		this->MaxPara.p = Maxp;

	}
	inline void DrawOutputXP() {
		glProtColor;
		glBegin(GL_LINE_STRIP);
		for (auto i = 0; i < Rung_global_pointer->Output.size(); i++) {
			glVertex2d(Rung_global_pointer->Output[i].x / Range, Rung_global_pointer->Output[i].p / Range);
		}
		glEnd();
	}
	inline void DrawOutputTX() {
		glProtColor;
		// ê¸ÇÃï`âÊ
		glBegin(GL_LINE_STRIP);
		for (auto i = 0; i < Rung_global_pointer->Output.size(); i++) {
			glVertex2d(Rung_global_pointer->Output[i].t, Rung_global_pointer->Output[i].x);
		}
		glEnd();
	}
	inline void DrawOutputTP() {
		glProtColor;
		// ê¸ÇÃï`âÊ
		glBegin(GL_LINE_STRIP);
		for (auto i = 0; i < Rung_global_pointer->Output.size(); i++) {
			glVertex2d(Rung_global_pointer->Output[i].t, Rung_global_pointer->Output[i].p);
		}
		glEnd();
	}
	inline void DrawAxis(double (RungeKutta::*TypeScaleX)(int), double (RungeKutta::*TypeScaleY)(int)) {
		DrawScale(TypeScaleX, TypeScaleY);
		glAxisColor;
		glBegin(GL_LINES);
		glVertex2d(-1.0, 0.0);
		glVertex2d(1.0, 0.0);
		glVertex2d(0.0, 1.0);
		glVertex2d(0.0, -1.0);
		glEnd();

	}

	inline static void dispTX() {

		glClear(GL_COLOR_BUFFER_BIT);
		Rung_global_pointer->DrawOutputTX();
		Rung_global_pointer->DrawAxis(&RungeKutta::GetScaleT, &RungeKutta::GetScaleX);
		glutSwapBuffers();
	}
	inline static void dispXP() {

		glClear(GL_COLOR_BUFFER_BIT);
		Rung_global_pointer->DrawOutputXP();
		Rung_global_pointer->DrawAxis(&RungeKutta::GetScaleX, &RungeKutta::GetScaleP);
		glutSwapBuffers();
	}
	inline static void dispTP() {

		glClear(GL_COLOR_BUFFER_BIT);
		Rung_global_pointer->DrawOutputTP();
		Rung_global_pointer->DrawAxis(&RungeKutta::GetScaleT, &RungeKutta::GetScaleP);
		glutSwapBuffers();
	}
	inline void Draw() {
		int argc = 0;
		glutInit(&argc,NULL);
		glutInitWindowPosition(200, 100);
		glutInitWindowSize(400, 300);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
		winID[0] = glutCreateWindow("Duffing T-X");
		font = new GLFONT(L"ÇlÇrñæí©", 24);
		glutDisplayFunc(dispTX);
		glutInitWindowPosition(600, 100);
		glutInitWindowSize(400, 300);
		winID[1] = glutCreateWindow("Duffing T-P");
		glutDisplayFunc(dispTP);
		glutInitWindowPosition(200, 500);
		glutInitWindowSize(400, 300);
		winID[1] = glutCreateWindow("Duffing X-P");
		glutDisplayFunc(dispXP);
		glutMainLoop();
	}
};


void main() {
	RungeKutta Rung(0.01,1000.0,4);
	Rung_global_pointer = &Rung;
	InputRung Input;
	Input.p = 0.0;
	Input.x = 0.3;
	Input.t = 0.0;
	Cofficient InputCof;
	
	InputCof.SetCoff(0.01,0.2,1.0 ,1.00, 0.3, 1/(2*PI), 0);
	Rung.SetParam(InputCof,Input);
	Rung.CalcDffingEqudydx();
	Rung.OutputCSV();
	Rung.MakeGraphData();
	Rung.Draw();
}