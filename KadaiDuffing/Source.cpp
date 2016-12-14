#include "Common.h"
#include "StringOpenGL.h"
class DuffingSolv;
constexpr double NAN_DEF = std::numeric_limits<double>::has_quiet_NaN;

DuffingSolv	*Rung_global_pointer;
int winID[4];
#define glAxisColor			glColor3d(1.0,0.0,0.0)
#define glProtColor			glColor3d(0.0,1.0,0.0)
#define glBackGroundColor	glColor3d(0.0,0.0,0.0)
#define glScaleColor		glColor3d(1.0,1.0,1.0)
#define glScalePlotColor	glColor3d(0.5,0.5,0.5)
#define glMassColor			glColor3d(0.5,0.5,0.5)
void Idle() {
	glutPostRedisplay();
}
GLfloat top = -0.9;
void timer(int value) {


	glutTimerFunc(10, timer, 0);
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
	unsigned int Loop;
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
			double ValX;
			ValX = (this->*TypeScaleX)((int)i);
			if (ValX == NAN_DEF){
				return;
			}

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
		return (Output.size()*CoffList.Getdt())*i/10 ;
	}
	inline double GetScaleP(int i) {
		return MaxPara.p / Range*(i - Scale / 2);
	}
	inline double NotScale(int i) {
		return NAN_DEF;
	}
public:
	double Range;
	double Offset;
	template<class TYPE,class NUM> DuffingSolv(const TYPE dt, const TYPE TRange, const NUM NumSeg) {
		Output.resize((int)(TRange / dt + 1));
		Analyze.resize((int)(TRange / dt + 1));
		Solver.Init();
		SetCoffdfdt = FALSE;
		Range = InitRange;
		Offset = InitOffset;
		Loop = 0;
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
		Analyze[0] = InitStatus;
		if (SetCoffdfdt) {
			for (auto Itr = Analyze.begin(); Itr != Analyze.end(); Itr++) {
				Solver.CalcRungeSolveOnce(&(*Itr));
				CoffList.UpdateCoff();
				Solver.SetCoff(CoffList);
			}
		}
		else {
			Solver.CalcRungeSolve(&Analyze);
		}

	}
	inline void OutputCSV() {
		std::ofstream OutputFile("RungKutta001.csv",std::ofstream::trunc|std::ofstream::out);
		OutputFile << "t,x,p" << std::endl;
		//for (auto Rungitr = Output.begin(); Rungitr != Output.end(); Rungitr++) {
		//	OutputFile << (*Rungitr).t << "," << (*Rungitr).x << "," << (*Rungitr).p << std::endl;
		//}
		for (auto Rungitr =0; Rungitr < Analyze.size(); Rungitr++) {
			OutputFile << Analyze[Rungitr].t << "," << Analyze[Rungitr].x << "," << Analyze[Rungitr].p << std::endl;
		}


	}
	inline void MakeGraphData() {
		SearchMax();
		Quantum();
	}
	inline void Quantum() {
		auto Data = Analyze;
		auto Outitr = Output.begin();
		for (auto itr = Data.begin(); itr !=Data.end(); itr++) {
			(*Outitr).x = (*itr).x / (MaxPara.x / Range);
			(*Outitr).p = (*itr).p / (MaxPara.p / Range);
			(*Outitr).t = (*itr).t / (Data.size()*CoffList.Getdt() / 2.0) - 1.0;
			Outitr++;
		}

	}
	inline void SearchMax() {
		MaxPara.x = Analyze[0].x;
		MaxPara.p = Analyze[0].p;

		for (auto i = 1; i <  Analyze.size(); i++) {
			if (MaxPara.x <  Analyze[i].x)
			{
				MaxPara.x = Analyze[i].x;
			}
			if (MaxPara.p < Analyze[i].p)
			{
				MaxPara.p = Analyze[i].p;
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
	inline void DrawBall() {
		glPushMatrix();
		glMassColor;
		glTranslated(0.0, Output[Loop%Output.size()].x,0.3);
		glutSolidSphere(0.1, 40, 40);
		glPopMatrix();
	}
	inline void DrawTime() {
		std::string TimeDisp;
		glScaleColor;
		TimeDisp = std::to_string(Analyze[Loop%Output.size()].t);
		glRasterPos2d(-0.70,0.70);
		glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, reinterpret_cast<const unsigned char*>(TimeDisp.c_str()));

	}
	inline void DrawOutputModel() {
			//Å@ÇÃï`âÊ
		DrawBall();
		DrawTime();
		Loop++;
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
	inline static void dispModel() {

		glClear(GL_COLOR_BUFFER_BIT);
		Rung_global_pointer->DrawOutputModel();
		Rung_global_pointer->DrawAxis(&DuffingSolv::NotScale, &DuffingSolv::GetScaleX);
		glutSwapBuffers();
		glutPostRedisplay();

	}

	inline static void MouseWheel(int wheel_number, int direction, int x, int y) {
		if (direction == 1) { Rung_global_pointer->Range += 0.01; }
		else { Rung_global_pointer->Range -= 0.01; }
		Rung_global_pointer->Quantum();
		glutPostRedisplay();
	}
	inline static void MouseWheelX(int wheel_number, int direction, int x, int y) {
		if (direction == 1) { Rung_global_pointer->Offset += 0.1; }
		else { Rung_global_pointer->Offset -= 0.1; }
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
		glutDisplayFunc(dispTX);
		glutMouseWheelFunc(MouseWheelX);
		glutInitWindowPosition(600, 100);
		glutInitWindowSize(400, 300);
		winID[1] = glutCreateWindow("Duffing T-P");
		glutDisplayFunc(dispTP);
		glutMouseWheelFunc(MouseWheelX);
		glutInitWindowPosition(200, 500);
		glutInitWindowSize(400, 300);
		winID[2] = glutCreateWindow("Duffing X-P");
		glutDisplayFunc(dispXP);
		glutMouseWheelFunc(MouseWheel);
		glutInitWindowPosition(600, 500);
		glutInitWindowSize(400, 300);
		winID[3] = glutCreateWindow("Duffing Physhics Model");
		glutDisplayFunc(dispModel);
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