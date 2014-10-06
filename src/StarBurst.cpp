//
//
//	StarBurst.cpp   :::   This Is the Starburst XBMC Visualization
//	V0.75			Written by Dinomight
//					dylan@castlegate.net
//				
// 
//
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "types.h"
#include "timer.h"
#include <ctype.h>

#include <xbmc/xbmc_vis_dll.h>
#include <GL/gl.h>
#include <GL/glu.h>

typedef enum _WEIGHT {
	WEIGHT_NONE = 0,
	WEIGHT_A    = 1,
	WEIGHT_B    = 2,
	WEIGHT_C		= 3
} WEIGHT;

#define sqr(x) (x*x)


#define	FREQ_DATA_SIZE 1024			// size of frequency data wanted
#define MAX_BARS 720				// number of bars in the Spectrum
#define MIN_PEAK_DECAY_SPEED 0		// decay speed in dB/frame
#define MAX_PEAK_DECAY_SPEED 4
#define MIN_RISE_SPEED 0.01f		// fraction of actual rise to allow
#define MAX_RISE_SPEED 1
#define MIN_FALL_SPEED 0.01f		// fraction of actual fall to allow
#define MAX_FALL_SPEED 1
#define MIN_FREQUENCY 1				// allowable frequency range
#define MAX_FREQUENCY 24000
#define MIN_LEVEL 0					// allowable level range
#define MAX_LEVEL 96
#define TEXTURE_HEIGHT 256
#define TEXTURE_MID 128
#define TEXTURE_WIDTH 1
#define MAX_CHANNELS 2

	float	m_pScreen[MAX_BARS*2];			// Current levels on the screen
	float	m_pPeak[MAX_BARS*2];			// Peak levels
	float	m_pWeight[FREQ_DATA_SIZE/2+1];		// A/B/C weighted levels for speed
	float	m_pFreq[MAX_BARS*2];			// Frequency data

	
	int		m_iSampleRate;
	int		m_width;
	int		m_height;
	float		m_centerx;
	float		m_centery;
	float		m_fRotation=0.0f;
	float		angle=0.0f;

	float		startradius;	//radius at which to start each bar
	float		minbar;			//minimum length of a bar
	float		spinrate;		// rate at witch to spin vis
	double		timepassed;
	double		lasttime = 0.0;
	double		currenttime = 0.0;

	float				r1; //floats used for bar colors;
	float				g1;
	float				b1;
	float				a1;
	float				r2;
	float				g2;
	float				b2;
	float				a2;

	int				m_iBars;				// number of bars to draw
	bool			m_bLogScale;			// true if our frequency is on a log scale
	bool			m_bShowPeaks;			// show peaks?
	bool			m_bAverageLevels;		// show average levels?
	float			m_fPeakDecaySpeed;		// speed of decay (in dB/frame)
	float			m_fRiseSpeed;			// division of rise to actually go up
	float			m_fFallSpeed;			// division of fall to actually go up
	float			m_fMinFreq;				// wanted frequency range
	float			m_fMaxFreq;
	float			m_fMinLevel;			// wanted level range
	float			m_fMaxLevel;
	WEIGHT			m_Weight;				// weighting type to be applied
	bool			m_bMixChannels;			// Mix channels, or stereo?
	bool			m_bSeperateBars;



	#define POLE1 20.598997*20.598997	// for A/B/C weighting
	#define POLE2 12194.217*12194.217	// for A/B/C weighting
	#define POLE3 107.65265*107.65265	// for A weighting
	#define POLE4 737.86223*737.86223	// for A weighting
	#define POLE5 158.5*158.5			// for B weighting


// A structure for our custom vertex type
struct CUSTOMVERTEX
{
    float x, y, z, rhw; // The transformed position for the vertex
    CRGBA color;        // The vertex color
};


CUSTOMVERTEX g_Vertices[MAX_BARS*4];
CTimer gTimer;

int inline sign(float a) {return (a>0) ? 1 : -1;}

int htoi(char *str) /* Convert hex string to integer */
{
  unsigned int digit, number = 0;
  while (*str)
  {
    if (isdigit(*str))
      digit = *str - '0';
    else
      digit = tolower(*str)-'a'+10;
    number<<=4;
    number+=digit;
    str++;
  }
  return number;  
}

void CreateArrays()
{
  CUSTOMVERTEX nullV = { 0.0f, 0.0f, 0.0f, 0.0f, CRGBA(0, 0, 0, 0), }; 

  for (int i=0; i<m_iBars*2; i++)
  {
    m_pScreen[i] = 0.0f;
    m_pPeak[i] = 0.0f;
    m_pFreq[i] = 0.0f;

    g_Vertices[i*2] = nullV;
    g_Vertices[i*2+1] = nullV;

  }
  // and the weight array
  if (m_Weight == WEIGHT_NONE)
    return;
  // calculate the weights (squared)
  float f2;
  for (int i=0; i<FREQ_DATA_SIZE/2+1; i++)
  {
    f2 = (float)sqr((float)i*m_iSampleRate/FREQ_DATA_SIZE);
    if (m_Weight == WEIGHT_A)
      m_pWeight[i] = (float)sqr(POLE2*sqr(f2)/(f2+POLE1)/(f2+POLE2)/sqrt(f2+POLE3)/sqrt(f2+POLE4));
    else if (m_Weight == WEIGHT_B)
      m_pWeight[i] = (float)sqr(POLE2*f2*sqrt(f2)/(f2+POLE1)/(f2+POLE2)/sqrt(f2+POLE5));
    else  // m_Weight == WEIGHT_C
      m_pWeight[i] = (float)sqr(POLE2*f2/(f2+POLE1)/(f2+POLE2));
  }
}

void SetupCamera()
{
  //Here we will setup the camera.
  //The camera has three settings: "Camera Position", "Look at Position" and "Up Direction"
  //We have set the following:
  //Camera Position: (0, 0, -30)
  //Look at Position: (0, 0, 0)
  //Up direction: Y-Axis.

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 0.0, -50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

void SetupPerspective()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, 1.0, 1.0, 100.0);
}

void SetupRotation(float x, float y, float z)
{
  ////Here we will rotate our view around the x, y and z axis.
  glMatrixMode(GL_MODELVIEW);
  glRotatef(x/M_PI*180, 1.0, 0.0, 0.0);
  glRotatef(y/M_PI*180, 0.0, 1.0, 0.0);
  glRotatef(z/M_PI*180, 0.0, 0.0, 1.0);
}

extern "C" void SetDefaults()
{
  m_iBars = 50;
  m_bLogScale=false;
  m_fPeakDecaySpeed=0.5f;
  m_fRiseSpeed=0.5f;
  m_fFallSpeed=0.5f;
  m_Weight = WEIGHT_NONE;
  m_bMixChannels = true;
  m_fMinFreq = 80;
  m_fMaxFreq = 16000;
  m_fMinLevel = 1;
  m_fMaxLevel = 80;
  m_bShowPeaks = true;
  m_bAverageLevels = false;
  spinrate = 1.0f/30.0f;
  startradius = 0.0f;
  minbar = 0.0f;
  //inital color
  r2 = 255;
  g2 = 200;
  b2 = 0;
  a2 = 255;
  //finalColor
  r1 = 163;
  g1 = 192;
  b1 = 255;
  a1 = 255;
  // color Diff
  r2 -= r1;
  g2 -= g1;
  b2 -= b1;
  a2 -= a1;
}

bool InitGeometry()
{
    // Initialize three vertices for rendering a triangle
    CUSTOMVERTEX g_Vertices[] =
    {
        { 200.0f, 200.0f, 0.5f, 1.0f, CRGBA(0, 255, 0, 255) }, // x, y, z, rhw, color
        { 300.0f, 200.0f, 0.5f, 1.0f, CRGBA(0, 255, 0, 255) },
		{ 300.0f, 300.0f, 0.5f, 1.0f, CRGBA(0, 255, 0, 255) },
		{ 200.0f, 300.0f, 0.5f, 1.0f, CRGBA(0, 255, 0, 255) },
		{ 200.0f, 300.0f, 0.5f, 1.0f, CRGBA(0, 255, 0, 255) },

		};

    return true;
}

ADDON_STATUS ADDON_Create(void* hdl, void* props)
{
  if (!props)
    return ADDON_STATUS_UNKNOWN;

  VIS_PROPS* visProps = (VIS_PROPS*)props;

  m_width = visProps->width;
  m_height = visProps->height;
  m_centerx = m_width/2.0f + visProps->x;
  m_centery = m_height/2.0f + visProps->y;
  SetDefaults();
  CreateArrays();

  return ADDON_STATUS_OK;
}

extern "C" void Start(int iChannels, int iSamplesPerSec, int iBitsPerSample, const char* szSongName)
{
  m_iSampleRate = iSamplesPerSec;
  CreateArrays();

  InitGeometry();
  gTimer.Init();
}

extern "C" void AudioData(const float* pAudioData, int iAudioDataLength, float *pFreqData, int iFreqDataLength)
{
  if (iFreqDataLength>FREQ_DATA_SIZE)
    iFreqDataLength = FREQ_DATA_SIZE;
  // weight the data using A,B or C-weighting
  if (m_Weight != WEIGHT_NONE)
  {
    for (int i=0; i<iFreqDataLength+2; i+=2)
    {
      pFreqData[i]*=m_pWeight[i>>1];
      pFreqData[i+1]*=m_pWeight[i>>1];
    }
  }
  // Group data into frequency bins by averaging (Ignore the constant term)
  int jmin=2;
  int jmax;
  // FIXME:  Roll conditionals out of loop
  for (int i=0, iBin=0; i < m_iBars; i++, iBin+=2)
  {
    m_pFreq[iBin]=0.000001f;	// almost zero to avoid taking log of zero later
    m_pFreq[iBin+1]=0.000001f;
    if (m_bLogScale)
      jmax = (int) (m_fMinFreq*pow(m_fMaxFreq/m_fMinFreq,(float)i/m_iBars)/m_iSampleRate*iFreqDataLength + 0.5f);
    else
      jmax = (int) ((m_fMinFreq + (m_fMaxFreq-m_fMinFreq)*i/m_iBars)/m_iSampleRate*iFreqDataLength + 0.5f);
    // Round up to nearest multiple of 2 and check that jmin is not jmax
    jmax<<=1;
    if (jmax > iFreqDataLength) jmax = iFreqDataLength;
    if (jmax==jmin) jmin-=2;
    for (int j=jmin; j<jmax; j+=2)
    {
      if (m_bMixChannels)
      {
        if (m_bAverageLevels)
          m_pFreq[iBin]+=pFreqData[j]+pFreqData[j+1];
        else 
        {
          if (pFreqData[j]>m_pFreq[iBin])
            m_pFreq[iBin]=pFreqData[j];
          if (pFreqData[j+1]>m_pFreq[iBin])
            m_pFreq[iBin]=pFreqData[j+1];
        }
      }
      else
      {
        if (m_bAverageLevels)
        {
          m_pFreq[iBin]+=pFreqData[j];
          m_pFreq[iBin+1]+=pFreqData[j+1];
        }
        else
        {
          if (pFreqData[j]>m_pFreq[iBin])
            m_pFreq[iBin]=pFreqData[j];
          if (pFreqData[j+1]>m_pFreq[iBin+1])
            m_pFreq[iBin+1]=pFreqData[j+1];
        }
      }
    }
    if (m_bAverageLevels)
    {
      if (m_bMixChannels)
        m_pFreq[iBin] /=(jmax-jmin);
      else
      {
        m_pFreq[iBin] /= (jmax-jmin)/2;
        m_pFreq[iBin+1] /= (jmax-jmin)/2;
      }
    }
    jmin = jmax;
  }
  // Transform data to dB scale, 0 (Quietest possible) to 96 (Loudest)
  for (int i=0; i < m_iBars*2; i++)
  {
    m_pFreq[i] = 10*log10(m_pFreq[i]);
    if (m_pFreq[i] > MAX_LEVEL)
      m_pFreq[i] = MAX_LEVEL;
    if (m_pFreq[i] < MIN_LEVEL)
      m_pFreq[i] = MIN_LEVEL;
  }
}

extern "C" void Render()
{
  gTimer.Update();

  timepassed = gTimer.GetDeltaTime();

  float PI = 3.141592653589793f;
  float devisions = (2.0f*PI)/(m_iBars);
  float dwidth = devisions/2.3f;

  angle += (2.0f*PI)/(spinrate*1000)*(timepassed/250000.0);

  for (int i=0; i < m_iBars*2; i++)
  {
    // truncate data to the users range
    if (m_pFreq[i] > m_fMaxLevel)
      m_pFreq[i] = m_fMaxLevel;
    m_pFreq[i]-=m_fMinLevel;
    if (m_pFreq[i] < 0)
      m_pFreq[i] = 0;
    // Smooth out the movement
    if (m_pFreq[i] > m_pScreen[i])
      m_pScreen[i] += (m_pFreq[i]-m_pScreen[i])*m_fRiseSpeed;
    else
      m_pScreen[i] -= (m_pScreen[i]-m_pFreq[i])*m_fFallSpeed;
    // Work out the peaks
    if (m_pScreen[i] >= m_pPeak[i])
    {
      m_pPeak[i] = m_pScreen[i];
    }
    else
    {
      m_pPeak[i]-=m_fPeakDecaySpeed;
      if (m_pPeak[i] < 0)
        m_pPeak[i] = 0;
    }
  }

  if (angle >2.0f*PI)
    angle -= 2.0f*PI;
  float x1 = 0;
  float y1 = 0;
  float x2 = 0;
  float y2 = 0;
  float radius=0;
  int iChannels = m_bMixChannels ? 1 : 2;

  //	for (int j=0; j<iChannels; j++){
  int j = 0;

  int points = 4;
  float scaler = (m_height/2 - minbar - startradius)/(m_fMaxLevel - m_fMinLevel);
  CRGBA color1 = CRGBA((int)r1,(int)g1,(int)b1,(int)a1); 
  for (int i=0; i < m_iBars*2; i+=2)
  {
    radius =  m_pScreen[i+j] * scaler + minbar + startradius;

    x1 = sin(angle - dwidth) * radius;		
    y1 = cos(angle - dwidth) * radius;
    x2 = sin(angle + dwidth) * radius;		
    y2 = cos(angle + dwidth) * radius;		
    float x3 = sin(angle) * startradius;
    float y3 = cos(angle) * startradius;

    float colorscaler = ((m_pScreen[i+j])/(m_fMaxLevel - m_fMinLevel));

    CRGBA color2 = CRGBA((int)((colorscaler*r2)+r1),(int)((colorscaler*g2)+g1),(int)((colorscaler*b2)+b1),(int) ((colorscaler*a2)+a1));
    //color1 = color2;
    CUSTOMVERTEX b = { m_centerx + x3, m_centery + y3, 0.5f, 1.0f, color2, };
    CUSTOMVERTEX a1 = { m_centerx + x1, m_centery + y1, 0.5f, 1.0f, color2, }; 
    CUSTOMVERTEX a2 = { m_centerx + x2, m_centery + y2, 0.5f, 1.0f, color2, }; 
    g_Vertices[(((i+2)/2 -1)*points)] = b;
    g_Vertices[(((i+2)/2 -1)*points)+1] = a1;
    g_Vertices[(((i+2)/2 -1)*points)+2] = a2;
    g_Vertices[(((i+2)/2 -1)*points)+3] = b;

    angle += devisions;
  }

  glBegin(GL_TRIANGLE_STRIP);
  for (size_t i=0;i<MAX_BARS*4;++i)
  {
    glColor3f(g_Vertices[i].color.r, g_Vertices[i].color.g, g_Vertices[i].color.b);
    glVertex3f(g_Vertices[i].x, g_Vertices[i].y, g_Vertices[i].z);
  }
  glEnd();
}

extern "C" void ADDON_Stop()
{
}

//-- GetSubModules ------------------------------------------------------------
// Return any sub modules supported by this vis
//-----------------------------------------------------------------------------
extern "C" unsigned int GetSubModules(char ***names)
{
  return 0; // this vis supports 0 sub modules
}

//-- OnAction -----------------------------------------------------------------
// Handle XBMC actions such as next preset, lock preset, album art changed etc
//-----------------------------------------------------------------------------
extern "C" bool OnAction(long flags, const void *param)
{
  bool ret = false;
  return ret;
}

//-- GetPresets ---------------------------------------------------------------
// Return a list of presets to XBMC for display
//-----------------------------------------------------------------------------
extern "C" unsigned int GetPresets(char ***presets)
{
  return 0;
}

//-- GetPreset ----------------------------------------------------------------
// Return the index of the current playing preset
//-----------------------------------------------------------------------------
extern "C" unsigned GetPreset()
{
  return 0;
}

//-- IsLocked -----------------------------------------------------------------
// Returns true if this add-on use settings
//-----------------------------------------------------------------------------
extern "C" bool IsLocked()
{
  return false;
}

//-- GetInfo ------------------------------------------------------------------
// Tell XBMC our requirements
//-----------------------------------------------------------------------------
extern "C" void GetInfo(VIS_INFO* pInfo)
{
  pInfo->bWantsFreq = true;
  pInfo->iSyncDelay = 16;
}

//-- Destroy ------------------------------------------------------------------
// Do everything before unload of this add-on
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------
extern "C" void ADDON_Destroy()
{
}

//-- HasSettings --------------------------------------------------------------
// Returns true if this add-on use settings
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------
extern "C" bool ADDON_HasSettings()
{
  return false;
}

//-- GetStatus ---------------------------------------------------------------
// Returns the current Status of this visualisation
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------
extern "C" ADDON_STATUS ADDON_GetStatus()
{
  return ADDON_STATUS_OK;
}

//-- GetSettings --------------------------------------------------------------
// Return the settings for XBMC to display
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------
extern "C" unsigned int ADDON_GetSettings(ADDON_StructSetting ***sSet)
{
  return 0;
}

//-- FreeSettings --------------------------------------------------------------
// Free the settings struct passed from XBMC
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------

extern "C" void ADDON_FreeSettings()
{
}

//-- SetSetting ---------------------------------------------------------------
// Set a specific Setting value (called from XBMC)
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------
extern "C" ADDON_STATUS ADDON_SetSetting(const char *strSetting, const void* value)
{
  return ADDON_STATUS_OK;
}

//-- Announce -----------------------------------------------------------------
// Receive announcements from XBMC
// !!! Add-on master function !!!
//-----------------------------------------------------------------------------
extern "C" void ADDON_Announce(const char *flag, const char *sender, const char *message, const void *data)
{
}
