#include <GL/glew.h>
#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "tinyxml.h"
#include "tinystr.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

// Este desenha todo o array por cada frame, por isso é que os modelos são desenhados a mexer;
// faz o desenho de tudo uma vez por frame.

// Tem o timer agora para tratar das rotaçoes e translaçoes com time.
// 1 segundo são 1000 milisegundos e o gluGet obtém o tempo em milisegundos.
// Quanto maior for o time, mais lento são as rotações/translações, daí (timer/(1000*tempo))

// Se o tempo = 1, ou seja se quiser fazer uma rotação completa/andar a orbita inteira em 1 segundo, entao o l = (timer/1000) começa em 0 e demora 1 segundo a chegar a l=1.
// Se tempo = 10, entao l = (timer/1000*10) entao para chegar de 0 a 1 demora 10 segundos (por causa do /10).

// Rotação -> mesma fórmula só que é *360 para dar o angulo de rotação, ou seja, se tempo = 1, entao l = (timer/100)*360 demora 1 segundo a ter de l=0 a l=360.


using namespace std;
#define POINT_COUNT 5
float alfa = M_PI/2, beta = M_PI/6, raio = 200.0;
vector< std::vector<long double> > vert;
int i = 0;
float x,y,z;
vector<std::string> modelos;
vector<float> pontosT;
int nvertices;
int inicioVertices = 0;
int preenchido = 0;

float timer;

float temp = 0.0;
float tempAux = 0.0;

vector<GLfloat> vertices;

GLuint buffer;

int cont = 0;
int contSoma = 0;
float tempo = 1.0;

float pointx = 0.0;
float pointy = 0.0;
float pointz = 0.0;

void normalize(float *a) { 
	float l = sqrt(a[0]*a[0] + a[1] * a[1] + a[2] * a[2]); 
	a[0] = a[0]/l; 
	a[1] = a[1]/l; 
	a[2] = a[2]/l;
}

void cross(float *a, float *b, float *res) { 
	res[0] = a[1]*b[2] - a[2]*b[1]; 
	res[1] = a[2]*b[0] - a[0]*b[2]; 
	res[2] = a[0]*b[1] - a[1]*b[0];
}

void buildRotMatrix(float *x, float *y, float *z, float *m) {
	m[0] = x[0];
	m[1] = x[1]; 
	m[2] = x[2]; 
	m[3] = 0; 
	m[4] = y[0]; 
	m[5] = y[1];
	m[6] = y[2];
	m[7] = 0;
	m[8] = z[0]; 
	m[9] = z[1];
	m[10] = z[2];
	m[11] = 0; 
	m[12] = 0;
	m[13] = 0; 
	m[14] = 0; 
	m[15] = 1;
}

void getCatmullRomPoint(float t, int *indices, float *res,float *deriv) {

	// catmull-rom matrix
	float m[4][4] = {	{-0.5f,  1.5f, -1.5f,  0.5f},
						{ 1.0f, -2.5f,  2.0f, -0.5f},
						{-0.5f,  0.0f,  0.5f,  0.0f},
						{ 0.0f,  1.0f,  0.0f,  0.0f}};
						
	res[0] = 0.0; res[1] = 0.0; res[2] = 0.0;
	deriv[0] = 0.0; deriv[1] = 0.0; deriv[2] = 0.0;	
	float ax[4];
	float ay[4];
	float az[4];
	

	ax[0]=m[0][0]*pontosT[indices[0]*3+0]+m[0][1]*pontosT[indices[1]*3+0]+m[0][2]*pontosT[indices[2]*3+0]+m[0][3]*pontosT[indices[3]*3+0];
	ax[1]=m[1][0]*pontosT[indices[0]*3+0]+m[1][1]*pontosT[indices[1]*3+0]+m[1][2]*pontosT[indices[2]*3+0]+m[1][3]*pontosT[indices[3]*3+0];
	ax[2]=m[2][0]*pontosT[indices[0]*3+0]+m[2][1]*pontosT[indices[1]*3+0]+m[2][2]*pontosT[indices[2]*3+0]+m[2][3]*pontosT[indices[3]*3+0];
	ax[3]=m[3][0]*pontosT[indices[0]*3+0]+m[3][1]*pontosT[indices[1]*3+0]+m[3][2]*pontosT[indices[2]*3+0]+m[3][3]*pontosT[indices[3]*3+0];

	ay[0]=m[0][0]*pontosT[indices[0]*3+1]+m[0][1]*pontosT[indices[1]*3+1]+m[0][2]*pontosT[indices[2]*3+1]+m[0][3]*pontosT[indices[3]*3+1];
	ay[1]=m[1][0]*pontosT[indices[0]*3+1]+m[1][1]*pontosT[indices[1]*3+1]+m[1][2]*pontosT[indices[2]*3+1]+m[1][3]*pontosT[indices[3]*3+1];
	ay[2]=m[2][0]*pontosT[indices[0]*3+1]+m[2][1]*pontosT[indices[1]*3+1]+m[2][2]*pontosT[indices[2]*3+1]+m[2][3]*pontosT[indices[3]*3+1];
	ay[3]=m[3][0]*pontosT[indices[0]*3+1]+m[3][1]*pontosT[indices[1]*3+1]+m[3][2]*pontosT[indices[2]*3+1]+m[3][3]*pontosT[indices[3]*3+1];

	az[0]=m[0][0]*pontosT[indices[0]*3+2]+m[0][1]*pontosT[indices[1]*3+2]+m[0][2]*pontosT[indices[2]*3+2]+m[0][3]*pontosT[indices[3]*3+2];
	az[1]=m[1][0]*pontosT[indices[0]*3+2]+m[1][1]*pontosT[indices[1]*3+2]+m[1][2]*pontosT[indices[2]*3+2]+m[1][3]*pontosT[indices[3]*3+2];
	az[2]=m[2][0]*pontosT[indices[0]*3+2]+m[2][1]*pontosT[indices[1]*3+2]+m[2][2]*pontosT[indices[2]*3+2]+m[2][3]*pontosT[indices[3]*3+2];
	az[3]=m[3][0]*pontosT[indices[0]*3+2]+m[3][1]*pontosT[indices[1]*3+2]+m[3][2]*pontosT[indices[2]*3+2]+m[3][3]*pontosT[indices[3]*3+2];
	

	// resultado final do produto de matrizes do pdf
	res[0]=pow(t,3)*ax[0]+pow(t,2)*ax[1]+t*ax[2]+ax[3];
	res[1]=pow(t,3)*ay[0]+pow(t,2)*ay[1]+t*ay[2]+ay[3];
	res[2]=pow(t,3)*az[0]+pow(t,2)*az[1]+t*az[2]+az[3];

	deriv[0]=3*pow(t,2)*ax[0]+2*t*ax[1]+ax[2];
	deriv[1]=3*pow(t,2)*ay[0]+2*t*ay[1]+ay[2];
	deriv[2]=3*pow(t,2)*az[0]+2*t*az[1]+az[2];

	// Compute point res = T * M * P
	// where Pi = p[indices[i]]
	// ...
}

void getGlobalCatmullRomPoint(float gt, float *res, float *deriv) {

	float t = gt * cont; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment

	// indices store the points
	int indices[4]; 
	indices[0] = ((index + cont-1)%cont) + contSoma;	
	indices[1] = ((indices[0]+1)%cont) + contSoma;
	indices[2] = ((indices[1]+1)%cont) + contSoma; 
	indices[3] = ((indices[2]+1)%cont) + contSoma;

	getCatmullRomPoint(t, indices, res, deriv);
}

void renderCatmullRomCurve() {

// desenhar a curva usando segmentos de reta - GL_LINE_LOOP
	float res[3];
	float deriv[3];
	
	//gera a linha a partir dos pontos gerados acima
	glBegin(GL_LINE_LOOP);
	for(int i=0; i<100; i++){
		getGlobalCatmullRomPoint(i/100.0,res,deriv);
		glVertex3f(res[0],res[1],res[2]);


	}
	
	glEnd();
	

}

void changeSize(int w, int h) {
    
    // Prevent a divide by zero, when window is too short
    // (you cant make a window with zero width).
    if(h == 0)
        h = 1;
    
    // compute window's aspect ratio
    float ratio = w * 1.0 / h;
    
    // Set the projection matrix as current
    glMatrixMode(GL_PROJECTION);
    // Load Identity Matrix
    glLoadIdentity();
    
    // Set the viewport to be the entire window
    glViewport(0, 0, w, h);
    
    // Set perspective
    gluPerspective(45.0f ,ratio, 1.0f ,1000.0f);
    
    // return to the model view matrix mode
    glMatrixMode(GL_MODELVIEW);
}

void processKeys(unsigned char c, int xx, int yy) {
  switch (c){
  case 'q': raio=raio+1.0; break;
  case 'e': raio=raio-1.0; break;
  case 'w': if((beta+0.1<=M_PI/2) && (beta+0.1>=-M_PI/2)) beta+=0.1; break;
  case 's': if((beta-0.1<=M_PI/2) && (beta-0.1>=-M_PI/2)) beta-=0.1; break;
  case 'm': alfa+=0.1; break;
  case 'n': alfa-=0.1; break;
  default: break;
  }
  glutPostRedisplay();
}

// tamM é o indice inicial para desenhar e tamMax é o numero de indices do buffer que vai usar para desenhar
// pelo que parece também dá se usar de 0 a vertices.size()
void desenhaModelo(int tamM, int tamMax){

	glGenBuffers(1,&buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STATIC_DRAW);

	glVertexPointer(3,GL_FLOAT,0,0);
	glDrawArrays(GL_TRIANGLES, tamM, tamMax);



}

void renderScene(void) {

	static float t = 0;
	float res[3];
	float deriv[3];
	float m[16]; 
	float z[3];

	int tamB = 0;

	contSoma = 0;
	cont = 0;



    glEnableClientState(GL_VERTEX_ARRAY);

    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // set the camera
    glLoadIdentity();
    gluLookAt(raio*sin(alfa)*cos(beta),raio*sin(beta),raio*cos(beta)*cos(alfa),
              0.0,0.0,0.0,
              0.0f,1.0f,0.0f);

    for(int i=0;i<modelos.size();){
		if(strcmp("group",modelos.at(i).c_str())==0){
			glPushMatrix();
			i++;
		}
		
		else if(strcmp("/group",modelos.at(i).c_str())==0){
			glPopMatrix();
			i++;
		}

		else if(strcmp("translate",modelos.at(i).c_str())==0){
			contSoma = contSoma + cont;
			timer = glutGet(GLUT_ELAPSED_TIME);
			tempo = atof(modelos.at(i+1).c_str());
			cont = atoi(modelos.at(i+2).c_str());
			renderCatmullRomCurve();
			getGlobalCatmullRomPoint(timer/(tempo*1000),res,deriv);
			glTranslatef(res[0],res[1],res[2]);
			i+=3;
		}

		else if(strcmp("rotate",modelos.at(i).c_str())==0){
			timer = glutGet(GLUT_ELAPSED_TIME);
			float ang;
			temp = atof(modelos.at(i+1).c_str());
			/*printf("%f %f %f\n",timer,timer/(temp*1000)*360, timer/(temp*1000));
			if(tempAux==0.0) tempAux = temp;
			if(tempAux!=0.0) ang = (360/temp)*(temp - tempAux);*/
			//float ang = atof(modelos.at(i+1).c_str());
			float x = atof(modelos.at(i+2).c_str());
			float y = atof(modelos.at(i+3).c_str());
			float z = atof(modelos.at(i+4).c_str());
			printf("%f\n",timer/(temp*1000)*360);
			glRotatef(timer/(temp*1000)*360,x,y,z);
			i+=5;
		}

		else if(strcmp("scale",modelos.at(i).c_str())==0){
			float sx = atof(modelos.at(i+1).c_str());
			float sy = atof(modelos.at(i+2).c_str());
			float sz = atof(modelos.at(i+3).c_str());
			glScalef(sx,sy,sz);
			i+=4;
		}

		else if(strcmp("ring",modelos.at(i).c_str())==0){
            		float inr = atof(modelos.at(i+1).c_str());
            		float outr = atof(modelos.at(i+2).c_str());
            		int nside = atoi(modelos.at(i+3).c_str());
            		int ring = atoi(modelos.at(i+4).c_str());
            		glutSolidTorus(inr,outr,nside,ring);
            		i+=5;
        	}

		else if(strcmp("color",modelos.at(i).c_str())==0){
            		float r = atof(modelos.at(i+1).c_str());
           		float g = atof(modelos.at(i+2).c_str());
            		float b = atof(modelos.at(i+3).c_str());
            		glColor3f(r/255.0,g/255.0,b/255.0);
            		i+=4;
		}

		else if(strcmp("model",modelos.at(i).c_str())==0){
			long double p1,p2,p3,p4,p5,p6,p7,p8,p9;
			ifstream fich;
			fich.open(modelos.at(i+1).c_str());
			if(preenchido==0){
    		if (fich.is_open()){
			fich >> nvertices;
        		for(int k = 0; k < nvertices; k=k+3){
				fich >> p1;
				fich >> p2;
				fich >> p3;

				fich >> p4;
				fich >> p5;
				fich >> p6;

				fich >> p7;
				fich >> p8;
				fich >> p9;

				/*glBegin(GL_TRIANGLES);
				glVertex3f(p1, p2, p3);
				glVertex3f(p4, p5, p6);
				glVertex3f(p7, p8, p9);
				glEnd();*/

				vertices.push_back(p1);
				vertices.push_back(p2);
				vertices.push_back(p3);
				vertices.push_back(p4);
				vertices.push_back(p5);
				vertices.push_back(p6);
				vertices.push_back(p7);
				vertices.push_back(p8);
				vertices.push_back(p9);
				}
		     	}
		     }
		 fich.close();
		 desenhaModelo(tamB, nvertices);
		 tamB += nvertices;
		 i=i+2;
	     }
    }
    preenchido = 1;

    // End of frame
    glutSwapBuffers();
    //tempAux --;
    //t += 0.001;
}


void lerXML(TiXmlElement* e){

	const char* pAttrib;
	const char* s[3] = {"X", "Y", "Z"};
	const char* t[4] = {"time", "axisX", "axisY", "axisZ"};
	const char* mod[4] = {"file", "diffR", "diffG", "diffB"};

	while(e){				
			if(strcmp("group",e->Value()) == 0){
				if(e==NULL) printf("Erro no group.\n");
				else{ modelos.push_back(e->Value()); lerXML(e->FirstChildElement()); modelos.push_back("/group"); e=e->NextSiblingElement(); }
			}
				
				
			else if(strcmp("models",e->Value()) == 0){
				if(e==NULL) printf("Erro no models.\n");
				TiXmlElement* m = e->FirstChildElement("model");
				while(m){
					modelos.push_back("color");
					for(int i=1;i<4;i++){
						pAttrib = m->Attribute(mod[i]);
						if(pAttrib)
							modelos.push_back(pAttrib);
						else modelos.push_back("255");
					}
					modelos.push_back(m->Value());
					modelos.push_back(m->Attribute(mod[0]));
					m=m->NextSiblingElement();
				}
				e=e->NextSiblingElement();
			}

			else if(strcmp("ring",e->Value()) == 0){
				if(e==NULL) printf("Erro no translate.\n");
				modelos.push_back(e->Value());
				TiXmlAttribute* pAttrib=e->FirstAttribute();
				while(pAttrib){
					modelos.push_back(pAttrib->Value());
					pAttrib=pAttrib->Next();
				}
				e=e->NextSiblingElement();
			}

			else if(strcmp("translate",e->Value()) == 0){
				char str[15];
				int cont = 0;
				if(e==NULL) printf("Erro no translate.\n");
				TiXmlElement* m = e->FirstChildElement("point");
				if(m==NULL) printf("Erro no translate \n");
				modelos.push_back(e->Value());
				modelos.push_back(e->Attribute("time"));
				while(m){
					cont++;
				for(int i=0;i<3;i++){
					pAttrib = m->Attribute(s[i]);
					if(pAttrib){
						pontosT.push_back(atof(pAttrib));}
					else pontosT.push_back(0.0);
				}
				m=m->NextSiblingElement();
				}
				sprintf(str, "%d", cont);
				modelos.push_back(str);
				e=e->NextSiblingElement();
			}

			else if(strcmp("scale",e->Value()) == 0){
				if(e==NULL) printf("Erro no scale.\n");
				modelos.push_back(e->Value());
				for(int i=0;i<3;i++){
					pAttrib = e->Attribute(s[i]);
					if(pAttrib)
						modelos.push_back(pAttrib);
					else modelos.push_back("0.0");
				}
				e=e->NextSiblingElement();
			}

			else if(strcmp("rotate",e->Value()) == 0){
				if(e==NULL) printf("Erro no rotate.\n");
				modelos.push_back(e->Value());
				for(int i=0;i<4;i++){
					pAttrib = e->Attribute(t[i]);
					if(pAttrib)
						modelos.push_back(pAttrib);
					else modelos.push_back("0.0");
				}
				e=e->NextSiblingElement();
			}

	}

}

	


int main(int argc, char* argv[]) {

	TiXmlDocument doc(argv[1]);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		TiXmlElement* titleElement = doc.FirstChildElement( "scene" )->FirstChildElement( "group" );
		TiXmlElement* e = titleElement->FirstChildElement();
		lerXML(e);
	}
	else printf("Failed to load file\n");

	//for(int i=0;i<modelos.size();i++) cout << modelos.at(i) << "\n";
    
    // init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(1200,800);
    glutCreateWindow("Scene");
    
    // Required callback registry 
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutIdleFunc(renderScene);
    
    // Callback registration for keyboard processing
    //glutSpecialFunc(processSpecialKeys);
    glutKeyboardFunc(processKeys);
    
    //  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glClearColor(0.0f,0.0f,0.0f,0.0f);
    
    if( glewInit() != GLEW_OK )
    {
        printf( "Failed to init glew!\n" );
        return -1;
    }

    // enter GLUT's main cycle
    glutMainLoop();
    
    return 1;
}
