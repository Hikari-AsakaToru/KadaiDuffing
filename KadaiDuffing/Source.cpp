#include "Common.h"
#include "StringOpenGL.h"
class DuffingSolv;
DuffingSolv	*Rung_global_pointer;
GLFONT		*font;
int winID[3];
#define glAxisColor			glColor3d(1.0,0.0,0.0)
#define glProtColor			glColor3d(0.0,1.0,0.0)
#define glBackGroundColor	glColor3d(0.0,0.0,0.0)
#define glScaleColor		glColor3d(1.0,1.0,1.0)
#define glScalePlotColor	glColor3d(0.5,0.5,0.5)
void Idle() {
	glutPostRedisplay();
}


class DuffingSolv {
	RungeKuttaSolve Solver;
	std::vector<RungState> Output;
	std::vector<RungState> Analyze;
	CoffcientFunc CoffList;
	InputRung InitStatus;
	Max MaxPara;
	bool SetCoffdfdt;
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
	inline void DrawScaleX(double (DuffingSolv::*TypeScaleX)(int)) {
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
	inline void DrawScaleY(double (DuffingSolv::*TypeScaleY)(int )) {
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
					Tmp = "    " + Tmp;
				}
				glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10, reinterpret_cast<const unsigned char*>(Tmp.c_str()));
			}
		}
		DrawScalePlotLineY();


	}
	inline void DrawScale(double (DuffingSolv::*TypeScaleX)(int i), double (DuffingSolv::*TypeScaleY)(int i)) {
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
	double Range;
	template<class TYPE,class NUM> DuffingSolv(const TYPE dt, const TYPE TRange, const NUM NumSeg) {
		Output.resize((int)(TRange / dt + 1));
		Analyze.resize((int)(TRange / dt + 1));
		Solver.Init();
		SetCoffdfdt = FALSE;
		Range = InitRange;
	}
	inline void SetParam(Cofficient SetCoff, const InputRung Inp) {
		CoffList.SetCoff(SetCoff);
		InitStatus = Inp;
		Solver.SetCoff(SetCoff);
	}
	inline void SetParamdfdt(Cofficient SetCoff, const InputRung Inp) {
		 auto Div = (int)(CoffList.Getdt()/SetCoff.Getdt());
		CoffList.Init(Div, &InitStatus);
		CoffList.SetdCoffdfdt(SetCoff);
		CoffList.SelectCalcCoff();
		SetCoffdfdt = TRUE;
	}
	inline void CalcDffingEqudydx() {
		Output[0] = InitStatus;
		if (SetCoffdfdt) {
			for (auto Itr = Output.begin(); Itr != Output.end(); Itr++) {
				Solver.CalcRungeSolveOnce(&(*Itr));
				CoffList.UpdateCoff();
				Solver.SetCoff(CoffList);
			}
		}
		else {
			Solver.CalcRungeSolve(&Output);
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
		SearchMax();
		Quantum();
	}
	inline void Quantum() {
		for (auto i = 0; i < Output.size(); i++) {
			Output[i].x = Output[i].x / (MaxPara.x / Range);
			Output[i].p = Output[i].p / (MaxPara.p / Range);
			Output[i].t = Output[i].t / (Output.size()*CoffList.Getdt() / 2.0) - 1.0;
		}

	}
	inline void SearchMax() {
		MaxPara.x = Output[0].x;
		MaxPara.p = Output[0].p;

		for (auto i = 1; i < Output.size(); i++) {
			if (MaxPara.x < Output[i].x)
			{
				MaxPara.x = Output[i].x;
			}
			if (MaxPara.p < Output[i].p)
			{
				MaxPara.p = Output[i].p;
			}
		}
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
	inline void DrawAxis(double (DuffingSolv::*TypeScaleX)(int), double (DuffingSolv::*TypeScaleY)(int)) {
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
		Rung_global_pointer->DrawAxis(&DuffingSolv::GetScaleT, &DuffingSolv::GetScaleX);
		glutSwapBuffers();
	}
	inline static void dispXP() {

		glClear(GL_COLOR_BUFFER_BIT);
		Rung_global_pointer->DrawOutputXP();
		Rung_global_pointer->DrawAxis(&DuffingSolv::GetScaleX, &DuffingSolv::GetScaleP);
		glutSwapBuffers();
	}
	inline static void dispTP() {

		glClear(GL_COLOR_BUFFER_BIT);
		Rung_global_pointer->DrawOutputTP();
		Rung_global_pointer->DrawAxis(&DuffingSolv::GetScaleT, &DuffingSolv::GetScaleP);
		glutSwapBuffers();
	}
	inline static void MouseWheel(int wheel_number, int direction, int x, int y){
		if (direction == 1) { Rung_global_pointer->Range += 0.01; }
		else { Rung_global_pointer->Range -= 0.01; }
		Rung_global_pointer->Quantum();
		glutPostRedisplay();
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
		glutMouseWheelFunc(MouseWheel);
		glutInitWindowPosition(600, 100);
		glutInitWindowSize(400, 300);
		winID[1] = glutCreateWindow("Duffing T-P");
		glutDisplayFunc(dispTP);
		glutMouseWheelFunc(MouseWheel);
		glutInitWindowPosition(200, 500);
		glutInitWindowSize(400, 300);
		winID[1] = glutCreateWindow("Duffing X-P");
		glutDisplayFunc(dispXP);
		glutMouseWheelFunc(MouseWheel);
		glutMainLoop();
	}
};



void main() {
	DuffingSolv Duff(0.01, 1000.0, 4);
	Rung_global_pointer = &Duff;
	InputRung Input;
	Input.p = 0.0;
	Input.x = 0.3;
	Input.t = 0.0;
	Cofficient InputCof;
	Cofficient InputCofdfdt;

//	InputCof.SetCoff(0.01, 0.2, 1.0, 1.00, 0.3, 1 / (2 * PI), 0);
	InputCof.SetCoff(0.01, 0.2, 1.0, 1.00, 0.3, 0.1 / (2 * PI), 0);
	Duff.SetParam(InputCof, Input);
	InputCofdfdt.SetCoff(0.001,-0.001, 0.0, 0.0, 0.0, 0.0		,0);
	Duff.SetParamdfdt(InputCofdfdt, Input);
	Duff.CalcDffingEqudydx();
	Duff.OutputCSV();
	Duff.MakeGraphData();
	Duff.Draw();
}