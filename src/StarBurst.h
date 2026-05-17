/*
 *  Copyright (C) 2005-2026 Team Kodi (https://kodi.tv)
 *  Copyright (C) Dinomight (dylan@castlegate.net)
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 *  See LICENSE.md for more information.
 */

#pragma once

#include <chrono>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <kodi/addon-instance/Visualization.h>
#include <kodi/gui/gl/Shader.h>

class MRFFT;

typedef enum _WEIGHT
{
  WEIGHT_NONE = 0,
  WEIGHT_A = 1,
  WEIGHT_B = 2,
  WEIGHT_C = 3
} WEIGHT;

#define sqr(x) (x * x)

#define FREQ_DATA_SIZE 512 // size of frequency data wanted
#define MAX_BARS 256 // number of bars in the Spectrum
#define MIN_PEAK_DECAY_SPEED 0 // decay speed in dB/frame
#define MAX_PEAK_DECAY_SPEED 4
#define MIN_RISE_SPEED 0.01f // fraction of actual rise to allow
#define MAX_RISE_SPEED 1
#define MIN_FALL_SPEED 0.01f // fraction of actual fall to allow
#define MAX_FALL_SPEED 1
#define MIN_FREQUENCY 1 // allowable frequency range
#define MAX_FREQUENCY 24000
#define MIN_LEVEL 0 // allowable level range
#define MAX_LEVEL 96
#define TEXTURE_HEIGHT 256
#define TEXTURE_MID 128
#define TEXTURE_WIDTH 1
#define MAX_CHANNELS 2

#define POLE1 20.598997 * 20.598997 // for A/B/C weighting
#define POLE2 12194.217 * 12194.217 // for A/B/C weighting
#define POLE3 107.65265 * 107.65265 // for A weighting
#define POLE4 737.86223 * 737.86223 // for A weighting
#define POLE5 158.5 * 158.5 // for B weighting

class ATTR_DLL_LOCAL CVisualizationStarBurst : public kodi::addon::CAddonBase,
                                               public kodi::addon::CInstanceVisualization,
                                               public kodi::gui::gl::CShaderProgram
{
public:
  CVisualizationStarBurst();
  ~CVisualizationStarBurst() override = default;

  bool Start(int channels,
             int samplesPerSec,
             int bitsPerSample,
             const std::string& songName) override;
  void Stop() override;
  void Render() override;
  void AudioData(const float* audioData, size_t audioDataLength) override;
  int GetSyncDelay() override { return 16; }

  void OnCompiledAndLinked() override;
  bool OnEnabled() override;

private:
  bool InitGeometry();
  void CreateArrays();

  std::unique_ptr<MRFFT> m_transform;
  std::unique_ptr<float[]> m_freqData;
  size_t m_freqDataLength{0};
  size_t m_prevFreqDataLength{0};

  glm::mat4 m_modelProjMat;

#ifdef HAS_GL
  GLuint m_vertexVBO[2] = {0};
#endif
  GLint m_uModelProjMatrix{-1};
  GLint m_aPosition{-1};
  GLint m_aColor{-1};

  bool m_startOK = false;

  float m_pScreen[MAX_BARS * 2] = {0.0f}; // Current levels on the screen
  float m_pPeak[MAX_BARS * 2] = {0.0f}; // Peak levels
  float m_pWeight[FREQ_DATA_SIZE / 2 + 1] = {0.0f}; // A/B/C weighted levels for speed
  float m_pFreq[MAX_BARS * 2] = {0.0f}; // Frequency data

  int m_iSampleRate;
  int m_width;
  int m_height;
  float m_centerx;
  float m_centery;

  float m_fRotation{0.0f};
  float m_angle{0.0f};
  float startradius{0.0f}; //radius at which to start each bar
  float minbar{200.0f}; //minimum length of a bar
  float spinrate{1.0f / 3.0f}; // rate at witch to spin vis

  float m_r1{0.64f}; //floats used for bar colors;
  float m_g1{0.75f};
  float m_b1{1.0f};
  float m_a1{1.0f};
  float m_r2{1.0f - m_r1};
  float m_g2{0.785f - m_g1};
  float m_b2{0.0f - 1.0f};
  float m_a2{1.0f - m_a1};

  int m_iBars{40}; // number of bars to draw
  bool m_bLogScale{false}; // true if our frequency is on a log scale
  bool m_bShowPeaks{false}; // show peaks?
  bool m_bAverageLevels{false}; // show average levels?
  float m_fPeakDecaySpeed{0.5f}; // speed of decay (in dB/frame)
  float m_fRiseSpeed{0.5f}; // division of rise to actually go up
  float m_fFallSpeed{0.5f}; // division of fall to actually go up
  float m_fMinFreq{80}; // wanted frequency range
  float m_fMaxFreq{16000};
  float m_fMinLevel{0}; // wanted level range
  float m_fMaxLevel{0.2f};
  WEIGHT m_Weight{WEIGHT_NONE}; // weighting type to be applied
  bool m_bMixChannels{true}; // Mix channels, or stereo?

   // The transformed position for the vertex
  glm::vec4 m_positions[MAX_BARS * 4] = {glm::vec4(0.0f)};
  // The vertex color
  glm::vec4 m_colors[MAX_BARS * 4] = {glm::vec4(0.0f)};

  double m_oldTime{0.0};
};
