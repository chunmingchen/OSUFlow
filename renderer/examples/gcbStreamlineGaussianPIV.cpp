/*

gcbStreamline:
The is a demo to show how to generate streamlines by OSUFlow and how to render the generated streamlines.

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <list>
#include <iterator>
#include <vector>

#include "gcb.h"
#include "OSUFlow.h"
#include "LineRendererInOpenGL.h"

// Monte Carlo seed samples
const int samples = 300;

// max steps
const int max_steps = 500;

const float pseedx[3] = {30, 85, 140};
const float seedy=86, seedz=1;
const float seedwidth = 8;
const float seedstep = 2;

float alpha = .01;

char *szVecFilePath;	// ADD-BY-LEETEN 09/29/2012
OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListSeedTrace*> sl_list; 
float center[3], len[3]; 
// ADD-BY-LEETEN 07/07/2010-BEGIN
list<VECTOR4> liv4Colors;
// ADD-BY-LEETEN 07/07/2010-END
CLineRendererInOpenGL cLineRenderer;

vector<VECTOR3> genSeeds()
{
    vector<VECTOR3> seeds;
    int x,y, i, count=0;
    for (i=0; i<3; i++)
        for (y=-seedwidth; y<=seedwidth; y+=seedstep)
            for (x=-seedwidth; x<=seedwidth; x+=seedstep)
            {
                seeds.push_back(VECTOR3( pseedx[i]+x, seedy+y, seedz ));
            }
    return seeds;
}

////////////////////////////////////////////////////////////////////////////
void compute_streamlines() 
{
	LOG("");

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 

  //osuflow->SetRegularSeedPoints(from, to, seedsPerDim);

  vector<VECTOR3> seeds;
  seeds = genSeeds();
  int nSeeds = seeds.size();

  osuflow->SetSeedPoints(&seeds[0], nSeeds);

  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(1, 5); 
  osuflow->GenStreamLines(sl_list , FORWARD_DIR, max_steps, 0, samples);
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 

	// ADD-BY-LEETEN 07/07/2010-BEGIN
    int n = sl_list.size(); // Jimmy added: list.size() takes O(n) time
	for(int i = 0; i < n; i++)
	{
		VECTOR4 v4Color;
        switch((i/samples)%9)
		{

        case 0: v4Color = VECTOR4(1.0f, 0.0f, 0.0f, 1.f);	break;
        case 1: v4Color = VECTOR4(1.0f, 0.5f, 0.0f, 1.f);	break;
        case 2: v4Color = VECTOR4(1.0f, 1.0f, 0.0f, 1.f);	break;
        //case 3: v4Color = VECTOR4(0.5f, 1.0f, 0.0f, 1.f);	break;
        case 3: v4Color = VECTOR4(0.0f, 1.0f, 0.0f, 1.f);	break;
        case 4: v4Color = VECTOR4(0.0f, 1.0f, 0.5f, 1.f);	break;
        case 5: v4Color = VECTOR4(0.0f, 1.0f, 1.0f, 1.f);	break;
        //case 7: v4Color = VECTOR4(0.0f, 0.5f, 1.0f, 1.f);	break;
        case 6: v4Color = VECTOR4(0.0f, 0.0f, 1.0f, 1.f);	break;
        case 7: v4Color = VECTOR4(0.5f, 0.0f, 1.0f, 1.f);	break;
        case 8: v4Color = VECTOR4(1.0f, 0.0f, 1.0f, 1.f);	break;
        //case 11: v4Color = VECTOR4(1.0f, 0.0f, 0.5f, 1.f);	break;
#if 0
		case 0: v4Color = VECTOR4(1.0f, 0.0f, 0.0f, 1.0f);	break;
		case 1: v4Color = VECTOR4(0.0f, 1.0f, 0.0f, 1.0f);	break;
		case 2: v4Color = VECTOR4(0.0f, 0.0f, 1.0f, 1.0f);	break;
		case 3: v4Color = VECTOR4(1.0f, 1.0f, 0.0f, 1.0f);	break;
		case 4: v4Color = VECTOR4(1.0f, 0.0f, 1.0f, 1.0f);	break;
		case 5: v4Color = VECTOR4(0.0f, 1.0f, 1.0f, 1.0f);	break;
		case 6: v4Color = VECTOR4(1.0f, 1.0f, 1.0f, 1.0f);	break;
#endif
		}
        liv4Colors.push_back(v4Color*.25);
	}
	// ADD-BY-LEETEN 07/07/2010-END
	cLineRenderer._Update();
}


void draw_cube(float w, float h, float d, float r, float g, float b)
{
  glPushMatrix();
  glScalef(w,h,d);
  glColor3f(r, g, b);
  glutWireCube(1.0);   // draw a wire cube
  glPopMatrix();
}

void draw_streamlines() 
{
    glPushAttrib(
		GL_LIGHTING_BIT |
		0
	);

    cLineRenderer._SetBoundingBox(minLen[0], minLen[1], minLen[2], maxLen[0], maxLen[1], maxLen[2]);


    // draw transparent traces
    glEnable (GL_BLEND);
    glDepthMask(GL_FALSE);
    glBlendFunc(GL_CONSTANT_ALPHA, GL_ONE);
    glBlendColor(1,1,1,alpha);
    //glDisable(GL_DEPTH_TEST);
    {
        cLineRenderer._Draw();
    }
    glDisable (GL_BLEND);
    //glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    if (0)
    {
        // show axes
        glBegin(GL_LINES);
        glColor3f(1,0,0);
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);
        glColor3f(0,1,0);
        glVertex3f(0,0,0);
        glVertex3f(0,1,0);
        glColor3f(0,0,1);
        glVertex3f(0,0,0);
        glVertex3f(0,0,1);
        glEnd();
    }

    float maxlen = max(len[0], max(len[1], len[2])) *.5;
    draw_cube(len[0]/maxlen, len[1]/maxlen, len[2]/maxlen,  .5, .5, .5);


	glPopAttrib();

}

///////////////////////////////////////////////////////////////////////////////
void
_KeyboardFunc(unsigned char ubKey, int iX, int iY)
{
	switch(ubKey)
	{
	case 's':
		compute_streamlines();
		glutPostRedisplay();
		break;

	case 'h':
		{
			int iHalo;
			cLineRenderer._GetInteger(CLineRenderer::ENABLE_HALO, &iHalo);
			iHalo = !iHalo;
			cLineRenderer._SetInteger(CLineRenderer::ENABLE_HALO, iHalo);
		}
		glutPostRedisplay();
		break;

	case 'l':
		{
			int iLighting;
			cLineRenderer._GetInteger(CLineRenderer::ENABLE_LIGHTING, &iLighting);
			iLighting = !iLighting;
			cLineRenderer._SetInteger(CLineRenderer::ENABLE_LIGHTING, iLighting);
		}

		glutPostRedisplay();
		break;

	// ADD-BY-LEETEN 09/29/2012-BEGIN
	case 'S':
		{
			VECTOR3 v3Min, v3Max;
			osuflow->Boundary(v3Min, v3Max);
			float pfDomainMin[4];
			float pfDomainMax[4];
			for(size_t d = 0; d < 3; d++)
			{
				pfDomainMin[d] = v3Min[d];
				pfDomainMax[d] = v3Max[d];
			}
			pfDomainMin[3] = 0.0f;
			pfDomainMax[3] = 0.0f;

			char szFilename[1024];
			strcpy(szFilename, szVecFilePath);
			strcat(szFilename, ".trace");

			OSUFlow::WriteFlowlines(
				pfDomainMin,
				pfDomainMax,
				&sl_list,
				NULL,
				szFilename);
			LOG(printf("Save the streamlines to %s", szFilename));
		}
		break;
	// ADD-BY-LEETEN 09/29/2012-END

    case '+':
        if (alpha<1.f)
            alpha /= .8f;
        printf("alpha=%f\n", alpha);
        glutPostRedisplay();
        break;
    case '-':
        if (alpha>1e-5)
            alpha *= .8f;
        printf("alpha=%f\n", alpha);
        glutPostRedisplay();
        break;

	}
}

void
_DisplayFunc()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render the scene
    draw_streamlines(); 

	// NOTE: Call glutSwapBuffers() at the end of your display function
	glutSwapBuffers();
}

void
init()
{
	LOG(printf("Initialize here."));
	glEnable(GL_DEPTH_TEST);

	// setup light 0
	static GLfloat pfLightAmbient[4] =	{0.0f, 0.0f, 0.0f, 1.0f};
	static GLfloat pfLightColor[4] =	{0.5f, 0.5f, 0.5f, 1.0f};
	glLightfv(GL_LIGHT0, GL_AMBIENT,		pfLightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,		pfLightColor);
	glLightfv(GL_LIGHT0, GL_SPECULAR,		pfLightColor);
	glLightf(GL_LIGHT0, GL_SPOT_EXPONENT,	4.0f);
	cLineRenderer._UpdateLighting();

	// ADD-BY-LEETEN 08/14/2010-BEGIN
	LOG(printf("The vector field is ready. Press key 's' to generate the primtives."));
	// ADD-BY-LEETEN 08/14/2010-END
}

void 
quit()
{
	LOG(printf("Clean up here."));
}

int
main(int argc, char* argv[])
{
	///////////////////////////////////////////////////////////////
	// when use GCB, it is still needed to initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_STENCIL );

	///////////////////////////////////////////////////////////////
	// initialize OSU flow
	osuflow = new OSUFlow(); 

	// load the scalar field
	LOG(printf("read file %s\n", argv[1])); 

    osuflow->LoadUncertainDataInGaussian(argv[1], argv[2], true); //true: a steady flow field
    osuflow->NormalizeField(true);

	szVecFilePath = argv[1];	// ADD-BY-LEETEN 09/29/2012

	// comptue the bounding box of the streamlines 
	VECTOR3 minB, maxB; 
	osuflow->Boundary(minLen, maxLen); // get the boundary 
	minB[0] = minLen[0]; minB[1] = minLen[1];  minB[2] = minLen[2];
	maxB[0] = maxLen[0]; maxB[1] = maxLen[1];  maxB[2] = maxLen[2];
	//  osuflow->SetBoundary(minB, maxB);  // set the boundary. just to test
										 // the subsetting feature of OSUFlow
	printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
								minLen[0], maxLen[0], minLen[1], maxLen[1], 
								minLen[2], maxLen[2]); 

	center[0] = (minLen[0]+maxLen[0])/2.0; 
	center[1] = (minLen[1]+maxLen[1])/2.0; 
	center[2] = (minLen[2]+maxLen[2])/2.0; 
	printf("center is at %f %f %f \n", center[0], center[1], center[2]); 
	len[0] = maxLen[0]-minLen[0]; 
	len[1] = maxLen[1]-minLen[1]; 
	len[2] = maxLen[2]-minLen[2]; 

	///////////////////////////////////////////////////////////////
	cLineRenderer._SetBoundingBox(
		minLen[0], minLen[1], minLen[2], 
		maxLen[0], maxLen[1], maxLen[2]);
	cLineRenderer._SetDataSource(&sl_list);
	// ADD-BY-LEETEN /2010-BEGIN
	cLineRenderer._SetColorSource(&liv4Colors);
    cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_TRACE);
    //cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_ALL_WHITE);
	// ADD-BY-LEETEN 07/07/2010-END

	///////////////////////////////////////////////////////////////
	glutCreateWindow("GCB Streamline");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbStreamline.cpp,v $
Revision 1.5  2010/10/01 20:38:23  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/08/15 04:23:06  leeten

[08/15/2010]
1. [ADD] Change the description of the application.
2. [ADD] Print a message a the end of _init() to indicate the user to press the hotkey 's'.

Revision 1.2  2010/07/07 18:00:33  leeten

[07/07/2010]
1. [ADD] Specify the color for all traces.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/