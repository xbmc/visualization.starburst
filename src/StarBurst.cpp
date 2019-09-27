//
//
//  StarBurst.cpp   :::   This Is the Starburst XBMC Visualization
//  V0.75      Written by Dinomight
//          dylan@castlegate.net
//
// 
//
//////////////////////////////////////////////////////////////////////

#include "StarBurst.h"

CVisualizationStarBurst::CVisualizationStarBurst()
{
  m_width = Width();
  m_height = Height();
  m_centerx = m_width/2.0f + X();
  m_centery = m_height/2.0f + Y();
  SetDefaults();
  CreateArrays();
}

bool CVisualizationStarBurst::Start(int channels, int samplesPerSec, int bitsPerSample, std::string songName)
{
  (void)channels;
  (void)bitsPerSample;
  (void)songName;

  std::string fraqShader = kodi::GetAddonPath("resources/shaders/" GL_TYPE_STRING "/frag.glsl");
  std::string vertShader = kodi::GetAddonPath("resources/shaders/" GL_TYPE_STRING "/vert.glsl");
  if (!LoadShaderFiles(vertShader, fraqShader) || !CompileAndLink())
  {
    kodi::Log(ADDON_LOG_ERROR, "Failed to create or compile shader");
    return false;
  }

  m_iSampleRate = samplesPerSec;
  CreateArrays();

  InitGeometry();

#ifdef HAS_GL
  glGenBuffers(2, m_vertexVBO);
#endif

  m_modelProjMat = glm::ortho(0.0f, float(Width()), float(Height()), 0.0f);

  m_oldTime = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
  m_startOK = true;
  return true;
}

void CVisualizationStarBurst::Stop()
{
  if (!m_startOK)
    return;

  m_startOK = false;

#ifdef HAS_GL
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glDeleteBuffers(2, m_vertexVBO);
  m_vertexVBO[0] = 0;
  m_vertexVBO[1]= 0;
#endif
}

void CVisualizationStarBurst::Render()
{
  if (!m_startOK)
    return;

  double currentTime = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
  double timepassed = currentTime - m_oldTime;
  m_oldTime = currentTime;

  float PI = 3.141592653589793f;
  float devisions = (2.0f*PI)/(m_iBars);
  float dwidth = devisions/2.3f;

  m_angle += (2.0f*PI)/(spinrate)*(timepassed/25.0);

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

  if (m_angle >2.0f*PI)
    m_angle -= 2.0f*PI;
  float x1 = 0;
  float y1 = 0;
  float x2 = 0;
  float y2 = 0;
  float radius=0;
  int iChannels = m_bMixChannels ? 1 : 2;

  //  for (int j=0; j<iChannels; j++){
  int j = 0;

  int points = 4;
  float scaler = (m_height/2 - minbar - startradius)/(m_fMaxLevel - m_fMinLevel);
  glm::vec4 color1 = glm::vec4(m_r1, m_g1, m_b1, m_a1);
  for (int i=0; i < m_iBars*2; i+=2)
  {
    radius =  m_pScreen[i+j] * scaler + minbar + startradius;

    x1 = sin(m_angle - dwidth) * radius;
    y1 = cos(m_angle - dwidth) * radius;
    x2 = sin(m_angle + dwidth) * radius;
    y2 = cos(m_angle + dwidth) * radius;
    float x3 = sin(m_angle) * startradius;
    float y3 = cos(m_angle) * startradius;

    float colorscaler = ((m_pScreen[i+j])/(m_fMaxLevel - m_fMinLevel));

    glm::vec4 color2 = glm::vec4(((colorscaler*m_r2)+m_r1),
                                 ((colorscaler*m_g2)+m_g1),
                                 ((colorscaler*m_b2)+m_b1),
                                 ((colorscaler*m_a2)+m_a1));
    //color1 = color2;
    glm::vec4 b = glm::vec4(m_centerx + x3, m_centery + y3, 0.5f, 1.0f);
    glm::vec4 a1 = glm::vec4(m_centerx + x1, m_centery + y1, 0.5f, 1.0f);
    glm::vec4 a2 = glm::vec4(m_centerx + x2, m_centery + y2, 0.5f, 1.0f);
    m_positions[(((i+2)/2 -1)*points)] = b;
    m_colors[(((i+2)/2 -1)*points)] = color2;
    m_positions[(((i+2)/2 -1)*points)+1] = a1;
    m_colors[(((i+2)/2 -1)*points)+1] = color2;
    m_positions[(((i+2)/2 -1)*points)+2] = a2;
    m_colors[(((i+2)/2 -1)*points)+2] = color2;
    m_positions[(((i+2)/2 -1)*points)+3] = b;
    m_colors[(((i+2)/2 -1)*points)+3] = color2;

    m_angle += devisions;
  }

#ifdef HAS_GL
  glBindBuffer(GL_ARRAY_BUFFER, m_vertexVBO[0]);
  glVertexAttribPointer(m_aPosition, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), BUFFER_OFFSET(offsetof(glm::vec4, x)));
  glEnableVertexAttribArray(m_aPosition);
  glBufferData(GL_ARRAY_BUFFER, sizeof(m_positions), m_positions, GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER, m_vertexVBO[1]);
  glVertexAttribPointer(m_aColor, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), BUFFER_OFFSET(offsetof(glm::vec4, r)));
  glEnableVertexAttribArray(m_aColor);
  glBufferData(GL_ARRAY_BUFFER, sizeof(m_colors), m_colors, GL_STATIC_DRAW);
#else
  glVertexAttribPointer(m_aPosition, 4, GL_FLOAT, GL_FALSE, 0, m_positions);
  glEnableVertexAttribArray(m_aPosition);

  glVertexAttribPointer(m_aColor, 4, GL_FLOAT, GL_FALSE, 0, m_colors);
  glEnableVertexAttribArray(m_aColor);
#endif

  glDisable(GL_BLEND);
  EnableShader();

  glDrawArrays(GL_TRIANGLE_STRIP, 0, m_iBars*4-2);

  DisableShader();
  glEnable(GL_BLEND);

  glDisableVertexAttribArray(m_aPosition);
  glDisableVertexAttribArray(m_aColor);
}

void CVisualizationStarBurst::AudioData(const float* pAudioData, int iAudioDataLength, float* pFreqData, int iFreqDataLength)
{
  if (iFreqDataLength>FREQ_DATA_SIZE)
    iFreqDataLength = FREQ_DATA_SIZE;
  // weight the data using A,B or C-weighting
  if (m_Weight != WEIGHT_NONE)
  {
    for (int i=0; i<iFreqDataLength+2; i+=2)
    {
      pFreqData[i] *= m_pWeight[i>>1];
      pFreqData[i+1] *= m_pWeight[i>>1];
    }
  }
  // Group data into frequency bins by averaging (Ignore the constant term)
  int jmin=2;
  int jmax;
  // FIXME:  Roll conditionals out of loop
  for (int i=0, iBin=0; i < m_iBars; i++, iBin+=2)
  {
    m_pFreq[iBin]=0.000001f;  // almost zero to avoid taking log of zero later
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
}

void CVisualizationStarBurst::OnCompiledAndLinked()
{
  m_uModelProjMatrix = glGetUniformLocation(ProgramHandle(), "u_modelViewProjectionMatrix");

  m_aPosition = glGetAttribLocation(ProgramHandle(), "a_position");
  m_aColor = glGetAttribLocation(ProgramHandle(), "a_color");
}

bool CVisualizationStarBurst::OnEnabled()
{
  glUniformMatrix4fv(m_uModelProjMatrix, 1, GL_FALSE, glm::value_ptr(m_modelProjMat));
  return true;
}

bool CVisualizationStarBurst::InitGeometry()
{
  // Initialize vertices for rendering a triangle

  glm::vec4 positions[] =
  {
    { glm::vec4(200.0f, 200.0f, 0.5f, 1.0f) }, // x, y, z, rhw, color
    { glm::vec4(300.0f, 200.0f, 0.5f, 1.0f) },
    { glm::vec4(300.0f, 300.0f, 0.5f, 1.0f) },
    { glm::vec4(200.0f, 300.0f, 0.5f, 1.0f) },
    { glm::vec4(200.0f, 300.0f, 0.5f, 1.0f) },
  };

  memcpy(m_positions, positions, sizeof(positions));

  glm::vec4 colors[] =
  {
    { glm::vec4(0.0f, 1.0f, 0.0f, 1.0f) }, // x, y, z, rhw, color
    { glm::vec4(0.0f, 1.0f, 0.0f, 1.0f) },
    { glm::vec4(0.0f, 1.0f, 0.0f, 1.0f) },
    { glm::vec4(0.0f, 1.0f, 0.0f, 1.0f) },
    { glm::vec4(0.0f, 1.0f, 0.0f, 1.0f) },
  };

  memcpy(m_colors, colors, sizeof(colors));

  return true;
}

void CVisualizationStarBurst::CreateArrays()
{
  for (int i=0; i<m_iBars*2; i++)
  {
    m_pScreen[i] = 0.0f;
    m_pPeak[i] = 0.0f;
    m_pFreq[i] = 0.0f;

    m_positions[i*2] = glm::vec4(0.0f);
    m_positions[i*2+1] = glm::vec4(0.0f);
    m_colors[i*2] = glm::vec4(0.0f);
    m_colors[i*2+1] = glm::vec4(0.0f);
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

void CVisualizationStarBurst::SetDefaults()
{
  m_iBars = 40;
  m_bLogScale=false;
  m_fPeakDecaySpeed = 0.5f;
  m_fRiseSpeed = 0.5f;
  m_fFallSpeed = 0.5f;
  m_Weight = WEIGHT_NONE;
  m_bMixChannels = true;
  m_fMinFreq = 80;
  m_fMaxFreq = 16000;
  m_fMinLevel = 0;
  m_fMaxLevel = 0.2;
  m_bShowPeaks = true;
  m_bAverageLevels = false;
  spinrate = 1.0/3.0;
  startradius = 0.0f;
  minbar = 200.0f;
  //inital color
  m_r2 = 1.0f;
  m_g2 = 0.785f;
  m_b2 = 0.0f;
  m_a2 = 1.0f;
  //finalColor
  m_r1 = 0.64f;
  m_g1 = 0.75f;
  m_b1 = 1.0f;
  m_a1 = 1.0f;
  // color Diff
  m_r2 -= m_r1;
  m_g2 -= m_g1;
  m_b2 -= m_b1;
  m_a2 -= m_a1;
}

ADDONCREATOR(CVisualizationStarBurst)
