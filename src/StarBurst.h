#pragma once

#include <chrono>
#include <kodi/addon-instance/Visualization.h>
#include <kodi/gui/gl/Shader.h>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef enum _WEIGHT {
  WEIGHT_NONE = 0,
  WEIGHT_A    = 1,
  WEIGHT_B    = 2,
  WEIGHT_C    = 3
} WEIGHT;

#define sqr(x) (x*x)

#define  FREQ_DATA_SIZE 512 // size of frequency data wanted
#define MAX_BARS 256 // number of bars in the Spectrum
#define MIN_PEAK_DECAY_SPEED 0    // decay speed in dB/frame
#define MAX_PEAK_DECAY_SPEED 4
#define MIN_RISE_SPEED 0.01f    // fraction of actual rise to allow
#define MAX_RISE_SPEED 1
#define MIN_FALL_SPEED 0.01f    // fraction of actual fall to allow
#define MAX_FALL_SPEED 1
#define MIN_FREQUENCY 1   // allowable frequency range
#define MAX_FREQUENCY 24000
#define MIN_LEVEL 0     // allowable level range
#define MAX_LEVEL 96
#define TEXTURE_HEIGHT 256
#define TEXTURE_MID 128
#define TEXTURE_WIDTH 1
#define MAX_CHANNELS 2

#define POLE1 20.598997*20.598997  // for A/B/C weighting
#define POLE2 12194.217*12194.217  // for A/B/C weighting
#define POLE3 107.65265*107.65265  // for A weighting
#define POLE4 737.86223*737.86223  // for A weighting
#define POLE5 158.5*158.5 // for B weighting

class ATTRIBUTE_HIDDEN CVisualizationStarBurst
  : public kodi::addon::CAddonBase,
    public kodi::addon::CInstanceVisualization,
    public kodi::gui::gl::CShaderProgram
{
public:
  CVisualizationStarBurst();
  ~CVisualizationStarBurst() override = default;

  bool Start(int channels, int samplesPerSec, int bitsPerSample, std::string songName) override;
  void Stop() override;
  void Render() override;
  void AudioData(const float* audioData, int audioDataLength, float* freqData, int freqDataLength) override;
  void GetInfo(bool& wantsFreq, int& syncDelay) override { wantsFreq = true; syncDelay = 16; }

  void OnCompiledAndLinked() override;
  bool OnEnabled() override;

private:
  bool InitGeometry();
  void CreateArrays();
  void SetDefaults();

  float m_fWaveform[2][512];

  glm::mat4 m_modelProjMat;

#ifdef HAS_GL
  GLuint m_vertexVBO[2] = { 0 };
#endif
  GLint m_uModelProjMatrix = -1;
  GLint m_aPosition = -1;
  GLint m_aColor = -1;

  bool m_startOK = false;

  float m_pScreen[MAX_BARS*2]; // Current levels on the screen
  float m_pPeak[MAX_BARS*2]; // Peak levels
  float m_pWeight[FREQ_DATA_SIZE/2+1]; // A/B/C weighted levels for speed
  float m_pFreq[MAX_BARS*2]; // Frequency data

  int m_iSampleRate;
  int m_width;
  int m_height;
  float m_centerx;
  float m_centery;
  float m_fRotation = 0.0f;
  float m_angle = 0.0f;

  float startradius; //radius at which to start each bar
  float minbar; //minimum length of a bar
  float spinrate; // rate at witch to spin vis

  float m_r1; //floats used for bar colors;
  float m_g1;
  float m_b1;
  float m_a1;
  float m_r2;
  float m_g2;
  float m_b2;
  float m_a2;

  int m_iBars; // number of bars to draw
  bool m_bLogScale; // true if our frequency is on a log scale
  bool m_bShowPeaks; // show peaks?
  bool m_bAverageLevels; // show average levels?
  float m_fPeakDecaySpeed; // speed of decay (in dB/frame)
  float m_fRiseSpeed; // division of rise to actually go up
  float m_fFallSpeed; // division of fall to actually go up
  float m_fMinFreq; // wanted frequency range
  float m_fMaxFreq;
  float m_fMinLevel; // wanted level range
  float m_fMaxLevel;
  WEIGHT m_Weight; // weighting type to be applied
  bool m_bMixChannels; // Mix channels, or stereo?
  bool m_bSeperateBars;

  glm::vec4 m_positions[MAX_BARS*4]; // The transformed position for the vertex
  glm::vec4 m_colors[MAX_BARS*4]; // The vertex color

  double m_oldTime;
};
